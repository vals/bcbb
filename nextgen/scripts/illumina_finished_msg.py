#!/usr/bin/env python
"""Script to check for finalized illumina runs and report to messaging server.

Run this script with an hourly cron job; it looks for newly finished output
directories for processing.

Usage:
    illumina_finished_msg.py <YAML local config>
                             [<post-processing config file>]

Supplying a post-processing configuration file skips the messaging step and
we moves directly into analysis processing on the current machine. Use
this if there is no RabbitMQ messaging server and your dump machine is directly
connected to the analysis machine. You will also want to set postprocess_dir in
the YAML local config to the directory to write fastq and analysis files.

The Galaxy config needs to have information on the messaging server and queues.

The local config should have the following information:

    dump_directories: directories to check for machine output
    msg_db: flat file of reported output directories
"""
import os
import operator
import socket
import glob
import getpass
import subprocess
import time
from optparse import OptionParser
from xml.etree.ElementTree import ElementTree
import errno

import logbook

from bcbio.solexa import samplesheet
from bcbio.log import create_log_handler, logger2
from bcbio import utils
from bcbio.distributed import messaging
from bcbio.solexa.flowcell import (get_flowcell_info, get_fastq_dir, get_qseq_dir)
from bcbio.pipeline.config_loader import load_config

LOG_NAME = os.path.splitext(os.path.basename(__file__))[0]
log = logbook.Logger(LOG_NAME)


def main(local_config, post_config_file=None,
         fetch_msg=True, process_msg=True, store_msg=True,
         backup_msg=False, qseq=True, fastq=True, remove_qseq=False,
         compress_fastq=False, casava=False):
    """Main function.
    """
    config = load_config(local_config)
    log_handler = create_log_handler(config, True)

    with log_handler.applicationbound():
        search_for_new(config, local_config, post_config_file,
                       fetch_msg, process_msg, store_msg, backup_msg,
                       qseq, fastq, remove_qseq, compress_fastq, casava)


def search_for_new(config, config_file, post_config_file,
                   fetch_msg, process_msg, store_msg, backup_msg, qseq, fastq,
                   remove_qseq, compress_fastq, casava):
    """Search for any new unreported directories.
    """
    reported = _read_reported(config["msg_db"])
    for dir_name in _get_directories(config["dump_directories"]):
        dir_name_is_reported = any(r_dir.startswith(dir_name) for r_dir in reported)

        if not os.path.isdir(dir_name) or dir_name_is_reported:
            continue  # Try the next dir_name

        # Injects run_name on logging calls.
        # Convenient for run_name on "Subject" for email notifications
        def inject_run_name(record):
            record.extra["run"] = os.path.basename(dir_name)

        with logbook.Processor(inject_run_name):
            if casava and _is_finished_dumping_read_1(dir_name):
                logger2.info("Generating fastq.qz files for reads 1 {}".format(dir_name))
                fastq_dir = None
                _generate_fastq_with_casava(dir_name, config, r1=True)

                post_process_arguments = { \
                "dname": dir_name, \
                "config": config, \
                "config_file": config_file, \
                "fastq_dir": fastq_dir, \
                "post_config_file": post_config_file, \
                "fetch_msg": True, \
                "process_msg": False, \
                "store_msg": store_msg, \
                "backup_msg": False \
                }

                _post_process_run(**post_process_arguments)

                continue  # Just move on to the next dir_name

            if not _is_finished_dumping(dir_name):
                continue

            logger2.info("The instrument has finished dumping on directory {}" + \
                         "".format(dir_name))
            _update_reported(config["msg_db"], dir_name)
            _process_samplesheets(dir_name, config)
            if qseq:
                logger2.info("Generating qseq files for {}".format(dir_name))
                _generate_qseq(get_qseq_dir(dir_name), config)

            fastq_dir = None
            if fastq:
                logger2.info("Generating fastq files for {}".format(dir_name))
                fastq_dir = _generate_fastq(dir_name, config)
                _calculate_md5(fastq_dir)
                if remove_qseq:
                    _clean_qseq(get_qseq_dir(dir_name), fastq_dir)

                if compress_fastq:
                    _compress_fastq(fastq_dir, config)

            if casava:
                logger2.info("Generating fastq.gz files for %s" % dir_name)
                _generate_fastq_with_casava(dir_name, config)

            post_process_arguments = { \
                "dname": dir_name, \
                "config": config, \
                "config_file": config_file, \
                "fastq_dir": fastq_dir, \
                "post_config_file": post_config_file, \
                "fetch_msg": fetch_msg, \
                "process_msg": process_msg, \
                "store_msg": store_msg, \
                "backup_msg": backup_msg \
                }

            _post_process_run(**post_process_arguments)

            # Update the reported database after successful processing
            _update_reported(config["msg_db"], dir_name)

        # Re-read the reported database to make sure it hasn't
        # changed while processing
        reported = _read_reported(config["msg_db"])


def _post_process_run(dname, config, config_file, fastq_dir, post_config_file,
                      fetch_msg, process_msg, store_msg, backup_msg):
    """For a finished directory, send out message or process directly.
    """
    run_module = "bcbio.distributed.tasks"
    # without a configuration file, send out message for processing
    if post_config_file is None:
        store_files, process_files, backup_files = _files_to_copy(dname)
        if process_msg:
            finished_message("analyze_and_upload", run_module, dname,
                             process_files, config, config_file)
        elif fetch_msg:
            finished_message("fetch_data", run_module, dname,
                             process_files, config, config_file)
        if store_msg:
            raise NotImplementedError("Storage server needs update.")
            finished_message("long_term_storage", run_module, dname,
                             store_files, config, config_file)
        if backup_msg:
            finished_message("backup_data", run_module, dname,
                             backup_files, config, config_file)
    # otherwise process locally
    else:
        analyze_locally(dname, post_config_file, fastq_dir)


def analyze_locally(dname, post_config_file, fastq_dir):
    """Run analysis directly on the local machine.
    """
    assert fastq_dir is not None
    post_config = load_config(post_config_file)
    analysis_dir = os.path.join(fastq_dir, os.pardir, "analysis")
    utils.safe_makedir(analysis_dir)
    with utils.chdir(analysis_dir):
        if post_config["algorithm"]["num_cores"] == "messaging":
            prog = post_config["analysis"]["distributed_process_program"]

        else:
            prog = post_config["analysis"]["process_program"]

        cl = [prog, post_config_file, dname]
        run_yaml = os.path.join(dname, "run_info.yaml")
        if os.path.exists(run_yaml):
            cl.append(run_yaml)

        subprocess.check_call(cl)


def _process_samplesheets(dname, config):
    """Process Illumina samplesheets into YAML files for post-processing.
    """
    ss_file = samplesheet.run_has_samplesheet(dname, config)
    if ss_file:
        out_file = os.path.join(dname, "run_info.yaml")
        logger2.info("CSV Samplesheet {0} found, converting to {1}" + \
        "".format(ss_file, out_file))
        samplesheet.csv2yaml(ss_file, out_file)


def _generate_fastq_with_casava(fc_dir, config, r1=False):
    """Perform demultiplexing and generate fastq.gz files for the current
    flowecell using CASAVA (>1.8).
    """
    basecall_dir = os.path.join(fc_dir, "Data", "Intensities", "BaseCalls")

    try:
        casava_dir = config["program"]["casava"]
    except KeyError as e:
        logger2.error("CASAVA directory not configured")
        raise e

    unaligned_dir = os.path.join(fc_dir, "Unaligned")
    samplesheet_file = samplesheet.run_has_samplesheet(fc_dir, config)
    num_mismatches = config["algorithm"].get("mismatches", 1)
    num_cores = config["algorithm"].get("num_cores", 1)

    cl = [os.path.join(casava_dir, "configureBclToFastq.pl")]
    cl.extend(["--input-dir", basecall_dir])
    cl.extend(["--output-dir", unaligned_dir])
    cl.extend(["--sample-sheet", samplesheet_file])
    cl.extend(["--mismatches", str(num_mismatches)])

    options = ["--fastq-cluster-count", "0", \
               "--ignore-missing-stats", \
               "--ignore-missing-bcl", \
               "--ignore-missing-control"]

    cl.extend(options)

    if r1:
        # Run configuration script
        logger2.info("Configuring BCL to Fastq conversion")
        logger2.debug(cl)
        subprocess.check_call(cl)

    # Go to <Unaligned> folder
    with utils.chdir(unaligned_dir):
        # Perform make

        cl = ["nohup", "make", "-j", str(num_cores)]
        if r1:
            cl.append("r1")

        logger2.info("Demultiplexing and converting bcl to fastq.gz")
        logger2.debug(cl)
        subprocess.check_call(cl)

    logger2.debug("Done")


def _generate_fastq(fc_dir, config):
    """Generate fastq files for the current flowcell.
    """
    fc_name, fc_date = get_flowcell_info(fc_dir)
    short_fc_name = "{}_{}".format(fc_date, fc_name)
    fastq_dir = get_fastq_dir(fc_dir)
    basecall_dir = os.path.split(fastq_dir)[0]
    postprocess_dir = config.get("postprocess_dir", "")
    if postprocess_dir:
        long_fc_name = os.path.basename(fc_dir)
        fastq_dir = os.path.join(postprocess_dir, long_fc_name, "fastq")

    if not fastq_dir == fc_dir:
        with utils.chdir(basecall_dir):
            qseq_files = glob.glob("*qseq.txt")
            lanes = (f.split("_")[1] for f in qseq_files)
            unique_lanes = sorted(set(lanes))
            lane_list = ",".join(unique_lanes)
            cl = ["solexa_qseq_to_fastq.py", short_fc_name, lane_list]

            if postprocess_dir:
                cl += ["-o", fastq_dir]

            logger2.debug("Converting qseq to fastq on all lanes.")

            cl = filter(lambda l: l != "", cl)
            subprocess.check_call(cl)

    return fastq_dir


def _calculate_md5(fastq_dir):
    """Calculate the md5sum for the fastq files
    """
    if not os.path.exists(fastq_dir):
        log.warn("There are no fastq files to hash")
        return

    glob_str = "*_fastq.txt"
    fastq_files = glob.glob(os.path.join(fastq_dir, glob_str))

    md5sum_file = os.path.join(fastq_dir, "md5sums.txt")

    with open(md5sum_file, 'w') as fh:
        for fastq_file in fastq_files:
            logger2.debug("Calculating md5 for %s using md5sum" % fastq_file)
            cl = ["md5sum", fastq_file]
            try:
                fh.write(subprocess.check_output(cl))
            except OSError as e:
                if e.errno == errno.ENOENT:
                    log.warn("md5sum not available")
                else:
                    raise e


def _compress_fastq(fastq_dir, config):
    """Compress the fastq files using gzip
    """
    glob_str = "*_fastq.txt"
    fastq_files = glob.glob(os.path.join(fastq_dir, glob_str))
    num_cores = config["algorithm"].get("num_cores", 1)
    active_procs = []
    sleeptime = 60
    # For testing without waiting too long.
    if hasattr(_compress_fastq, "sleeptime"):
        sleeptime = _compress_fastq.sleeptime

    for fastq_file in fastq_files:
        # Sleep for one minute while waiting for an open slot
        while len(active_procs) >= num_cores:
            time.sleep(sleeptime)
            active_procs, _ = _process_status(active_procs)

        logger2.debug("Compressing {} using gzip".format(fastq_file))
        cl = ["gzip", fastq_file]
        active_procs.append(subprocess.Popen(cl))

    # Wait for the last processes to finish
    while len(active_procs) > 0:
        time.sleep(sleeptime)
        active_procs, _ = _process_status(active_procs)


def _process_status(processes):
    """Return a list of the processes that are still active and
       a list of the returncodes of the processes that have finished
    """
    active = []
    retcodes = []
    for p in processes:
        c = p.poll()
        if c is None:
            active.append(p)
        else:
            retcodes.append(c)
    return active, retcodes


def _clean_qseq(bc_dir, fastq_dir):
    """Remove the temporary qseq files if the corresponding fastq file
       has been created
    """
    glob_str = "*_1_fastq.txt"
    fastq_files = glob.glob(os.path.join(fastq_dir, glob_str))

    for fastq_file in fastq_files:
        try:
            lane = int(os.path.basename(fastq_file)[0])
        except ValueError:
            continue

        logger2.debug("Removing qseq files for lane %d" % lane)
        glob_str = "s_%d_*qseq.txt" % lane

        for qseq_file in glob.glob(os.path.join(bc_dir, glob_str)):
            try:
                os.unlink(qseq_file)
            except:
                logger2.debug("Could not remove %s" % qseq_file)


def _generate_qseq(bc_dir, config):
    """Generate qseq files from illumina bcl files if not present.

    More recent Illumina updates do not produce qseq files. Illumina's
    offline base caller (OLB) generates these starting with bcl,
    intensity and filter files.
    """
    if not os.path.exists(os.path.join(bc_dir, "finished.txt")):
        bcl2qseq_log = os.path.join(config["log_dir"], "setupBclToQseq.log")
        cmd = os.path.join(config["program"]["olb"], "bin", "setupBclToQseq.py")
        cl = [cmd, "-L", bcl2qseq_log, "-o", bc_dir, "--in-place", "--overwrite",
              "--ignore-missing-stats", "--ignore-missing-control"]
        # in OLB version 1.9, the -i flag changed to intensities instead of input
        version_cl = [cmd, "-v"]
        p = subprocess.Popen(version_cl, stdout=subprocess.PIPE)
        (out, _) = p.communicate()
        olb_version = float(out.strip().split()[-1].rsplit(".", 1)[0])
        if olb_version > 1.8:
            cl += ["-P", ".clocs"]
            cl += ["-b", bc_dir]
        else:
            cl += ["-i", bc_dir, "-p", os.path.split(bc_dir)[0]]
        subprocess.check_call(cl)
        with utils.chdir(bc_dir):
            processors = config["algorithm"].get("num_cores", 8)
            cl = config["program"].get("olb_make", "make").split() + ["-j", str(processors)]
            subprocess.check_call(cl)


def _is_finished_dumping(directory):
    """Determine if the sequencing directory has all files.

    The final checkpoint file will differ depending if we are a
    single or paired end run.
    """
    # Check final output files; handles both HiSeq, MiSeq and GAII
    run_info = os.path.join(directory, "RunInfo.xml")
    hi_seq_checkpoint = \
    "Basecalling_Netcopy_complete_Read{}.txt".format(_expected_reads(run_info))

    to_check = ["Basecalling_Netcopy_complete_SINGLEREAD.txt",
                "Basecalling_Netcopy_complete_READ2.txt",
                hi_seq_checkpoint]

    # os.path.exists follows case sensitivity of the file system,
    # this method is always case sensitive
    dir_files = os.listdir(directory)
    dump_done = len(set(to_check) & set(dir_files)) > 0

    # Include a check to wait for any ongoing MiSeq analysis
    is_queued_for_analysis = os.path.exists(os.path.join(directory, "QueuedForAnalysis.txt"))
    job_is_completed = os.path.exists(os.path.join(directory, "CompletedJobInfo.xml"))

    miseq_analysis_checkpoint = not is_queued_for_analysis or job_is_completed

    return dump_done and miseq_analysis_checkpoint


def _is_finished_dumping_read_1(directory):
    """Determine if the sequencing directory has all files from read 1, as
    well as the indexed read (read 2).

    This lets CASAVA 1.8  and above start demultiplexing and converting to
    fastq.gz while the last read is still being processed.
    """
    indexed_read_checkpoint = os.path.join(directory, "Basecalling_Netcopy_complete_Read2.txt")
    read_1_finished_checkpoint = os.path.join(directory, "Demultiplexing_done_Read1.txt")

    if os.path.exists(indexed_read_checkpoint) \
    and not os.path.exists(read_1_finished_checkpoint):
        open(read_1_finished_checkpoint, 'w').close()
        return True

    else:
        return False


def _expected_reads(run_info_file):
    """Parse the number of expected reads from the RunInfo.xml file.
    """
    reads = []
    if os.path.exists(run_info_file):
        tree = ElementTree()
        tree.parse(run_info_file)
        read_elem = tree.find("Run/Reads")
        reads = read_elem.findall("Read")

    return len(reads)


def _is_finished_dumping_checkpoint(directory):
    """Recent versions of RTA (1.10 or better), write the complete file.

    This is the most straightforward source but as of 1.10 still does not
    work correctly as the file will be created at the end of Read 1 even
    if there are multiple reads.
    """
    check_file = os.path.join(directory, "Basecalling_Netcopy_complete.txt")
    check_v1, check_v2 = (1, 10)
    if os.path.exists(check_file):
        with open(check_file) as in_handle:
            line = in_handle.readline().strip()

        if line:
            version = line.split()[-1]
            v1, v2 = [float(v) for v in version.split(".")[:2]]
            if ((v1 > check_v1) or (v1 == check_v1 and v2 >= check_v2)):
                return True


def _files_to_copy(directory):
    """Retrieve files that should be remotely copied.
    """
    with utils.chdir(directory):
        image_redo_files = reduce(operator.add,
                                  [glob.glob("*.params"),
                                   glob.glob("Images/L*/C*"),
                                   ["RunInfo.xml", "runParameters.xml"]])

        qseqs = reduce(operator.add,
                    [glob.glob("Data/Intensities/*.xml"),
                     glob.glob("Data/Intensities/BaseCalls/*qseq.txt")
                    ])

        reports = reduce(operator.add,
                        [glob.glob("*.xml"),
                         glob.glob("Data/Intensities/BaseCalls/*.xml"),
                         glob.glob("Data/Intensities/BaseCalls/*.xsl"),
                         glob.glob("Data/Intensities/BaseCalls/*.htm"),
                         glob.glob("Unaligned/Basecall_Stats_*/*"),
                         glob.glob("Unalgiend/Basecall_Stats_*/**/*"),
                         ["Data/Intensities/BaseCalls/Plots", "Data/reports",
                          "Data/Status.htm", "Data/Status_Files", "InterOp"]
                        ])

        run_info = reduce(operator.add,
                        [glob.glob("run_info.yaml"),
                         glob.glob("*.csv"),
                         glob.glob("Unaligned/Project_*/**/*.csv"),
                         glob.glob("Unaligned/Undetermined_indices/**/*.csv"),
                         glob.glob("*.txt")
                        ])

        logs = reduce(operator.add, [["Logs", "Recipe", "Diag", "Data/RTALogs", "Data/Log.txt"]])

        fastq = reduce(operator.add,
                        [glob.glob("Data/Intensities/BaseCalls/*fastq.gz"),
                         glob.glob("Unaligned/Project_*/**/*.fastq.gz"),
                         glob.glob("Unaligned/Undetermined_indices/**/*.fastq.gz"),
                         ["Data/Intensities/BaseCalls/fastq"]
                        ])

        analysis = reduce(operator.add, [glob.glob("Data/Intensities/BaseCalls/Alignment")])

    return (sorted(image_redo_files + logs + reports + run_info + qseqs),
            sorted(reports + fastq + run_info + analysis),
            ["*"])


def _read_reported(msg_db):
    """Retrieve a list of directories previous reported.
    """
    reported = []
    if not os.path.exists(msg_db):
        return reported

    with open(msg_db) as in_handle:
        for line in in_handle:
            reported.append(line.strip())

    return reported


def _get_directories(dump_directories):
    if isinstance(dump_directories, basestring):
        dump_directories = [dump_directories]

    for directory in dump_directories:
        globs = []
        # Glob to capture general flowcells
        globs.extend(glob.glob(os.path.join(directory, "*[Aa]*[Xx][Xx]")))
        # Glob to capture MiSeq flowcells
        globs.extend(glob.glob(os.path.join(directory, "*_M*AMS*")))
        for dname in sorted(globs):
            if os.path.isdir(dname):
                yield dname


def _update_reported(msg_db, new_dname):
    """Add a new directory to the database of reported messages.
    """
    reported = _read_reported(msg_db)
    startswith = lambda dname: dname.startswith(new_dname)
    for r_dir in filter(startswith, reported):
        new_dname = r_dir
        reported.remove(r_dir)

    reported.append("{}\t{}".format(new_dname, time.strftime("%x-%X")))

    with open(msg_db, "w") as out_handle:
        for r_dir in reported:
            out_handle.write("{}\n".format(r_dir))


def finished_message(fn_name, run_module, directory, files_to_copy,
                     config, config_file):
    """Wait for messages with the give tag, passing on to the supplied handler.
    """
    logger2.debug("Calling remote function: {}".format(fn_name))
    user = getpass.getuser()
    hostname = socket.gethostbyaddr(socket.gethostname())[0]
    data = dict(
            machine_type='illumina',
            hostname=hostname,
            user=user,
            directory=directory,
            to_copy=files_to_copy
            )

    dirs = {"work": os.getcwd()}
    try:
        dirs["config"] = os.path.dirname(config_file)
    except AttributeError:
        log.warn("No config file was specified")

    runner = messaging.runner(run_module, dirs, config, config_file, wait=False)
    runner(fn_name, [[data]])


# Testing code: run with 'nosetests -v -s illumina_finished_message.py'
import unittest
import random
import string
import shutil


class IFMTestCase(unittest.TestCase):
    """Class with helper functions for testing illumina_finished_message.
    """
    def make_run_info_xml(self):
        test_xml_file = os.path.join(self.test_dir, "RunInfo.xml")

        test_xml = \
        """<?xml version="1.0"?>
            <RunInfo>
                <Run>
                    <Reads>
                        <Read Number="1" />
                        <Read Number="2" />
                    </Reads>
                </Run>
            </RunInfo>"""

        with open(test_xml_file, "w") as txf:
            txf.write(test_xml)


class FinishedDumpingTest(unittest.TestCase):
    """Tests for _is_finished_dumping
    """
    def setUp(self):
        self.test_dir = "".join(random.choice(string.ascii_uppercase) for i in xrange(9))
        os.mkdir(self.test_dir)

    def test__is_finished_dumping(self):
        """Test determining if a sequencer is finished dumping.
        """
        assert _is_finished_dumping(self.test_dir) == False, \
        "Empty directory was interpreted as finished"

        test_xml_file = os.path.join(self.test_dir, "RunInfo.xml")

        test_xml = \
        """<?xml version="1.0"?>
            <RunInfo>
                <Run>
                    <Reads>
                        <Read Number="1" />
                        <Read Number="2" />
                    </Reads>
                </Run>
            </RunInfo>"""

        with open(test_xml_file, "w") as txf:
            txf.write(test_xml)

        assert _is_finished_dumping(self.test_dir) == False, \
        "Existance of RunInfo.xml should not signal finished"

        with open(os.path.join(self.test_dir, "Basecalling_Netcopy_complete_Read2.txt"), "w") as tfh:
            tfh.write("")

        assert _is_finished_dumping(self.test_dir) == True, \
        """Existance of Basecalling_Netcopy_complete_Read2.txt should signal that
        the dump is finished."""

        miseq_queued = os.path.join(self.test_dir, "QueuedForAnalysis.txt")
        with open(miseq_queued, "w") as mqfh:
            mqfh.write("")

        assert _is_finished_dumping(self.test_dir) == False, \
        """Existance of QueuedForAnalysis.txt should prevent signalling the dump
        being finished."""

        job_completed = os.path.join(self.test_dir, "CompletedJobInfo.xml")
        with open(job_completed, "w") as jcfh:
            jcfh.write("")

        assert _is_finished_dumping(self.test_dir) == True, \
        """Existance of CompletedJobInfo.xml should counteract QueuedForAnalysis.txt,
        and signal that the dump is finished."""

    def tearDown(self):
        shutil.rmtree(self.test_dir)


class ReportedTest(unittest.TestCase):
    """Tests for _read_reported()
    """
    def setUp(self):
        self.temp_files = []

    def test__read_reported(self):
        """Test _read_reported()
        """
        test_file = "".join(random.choice(string.ascii_uppercase) for i in xrange(9))

        assert _read_reported(test_file) == [], \
        "Inputting nonexistant file does not return empty list"

        with open(test_file, "w") as test_handle:
            test_handle.write("TEST")

        self.temp_files.append(test_file)

        test_read = _read_reported(test_file)

        assert len(test_read) == 1, \
        "Unexpected number of lines in test file"

        assert test_read[0] == "TEST", \
        "Unexpected contents of the test file"

    def test__update_reported(self):
        """Test _update_reported()
        """
        test_file = "".join(random.choice(string.ascii_uppercase) for i in xrange(9))

        _update_reported(test_file, "Line 1")

        self.temp_files.append(test_file)

        reported = _read_reported(test_file)

        assert reported[0].startswith("Line 1"), \
        "Nonexistant report was not created and updated"

        reported_time = time.strptime(reported[0].split("\t")[-1], "%x-%X")
        assert isinstance(reported_time, time.struct_time), \
        "Could not parse report time."

        _update_reported(test_file, "Line 2")
        _update_reported(test_file, "Line 3")

        reported = _read_reported(test_file)

        n = len(reported)
        assert n == 3, \
        "Unexpected number of lines in test_file: {}".format(n)

        # Update already existing line "Line 1"
        _update_reported(test_file, "Line 1")

        reported = _read_reported(test_file)

        assert len(reported) == 3, \
        "Unexpected number of lines in test_file"

        for exp_line, line in zip(["Line 2", "Line 3", "Line 1"], reported):
            assert line.startswith(exp_line), \
            "Unmatching lines: " \
            "Expected '{}*', got '{}'".format(exp_line, line)

        _update_reported(test_file, "Line")

        reported = _read_reported(test_file)

        assert len(reported) == 1, \
        "File was not reduced to a single line."

        assert reported[0].startswith("Line"), \
        "Unmatching lines: Expected '{}*', got '{}'".format(exp_line, line)

        reported = reported[0].split("\t")
        assert len(reported) == 4, \
        "Unexpected number of reported dates: Expected 3, got {}".format(len(reported) - 1)

    def tearDown(self):
        for temp_file in self.temp_files:
            os.remove(temp_file)


class MainTest(IFMTestCase):
    """Tests for the main functionality of the script.
    """
    def setUp(self):
        self.test_dir = "test_data"
        self.bc_dir = "test_data/111009_SN1_0002_AB0CDDECXX/Data/Intensities/BaseCalls/"
        try:
            os.makedirs(self.bc_dir)

        except OSError as e:
            if e.errno == errno.EEXIST:
                pass
            else:
                raise e

        self.msg_db = "test_data/transferred.db"
        open(self.msg_db, "w").close()

        self.make_run_info_xml()
        open(os.path.join(self.test_dir, "111009_SN1_0002_AB0CDDECXX/Basecalling_Netcopy_complete_Read2.txt"), "w").close()

        self.kwords = {
            "config": {
                "msg_db": self.msg_db,
                "dump_directories": "test_data",
                "algorithm": {}
                },
            "config_file": None,
            "post_config_file": None,

            "fetch_msg": False,
            "process_msg": False,
            "store_msg": False,
            "backup_msg": False,

            "qseq": False,
            "fastq": False,
            "remove_qseq": False,
            "compress_fastq": False,
            "casava": False
            }

    def test_search_for_new_fastq_empty(self):
        self.kwords["fastq"] = True

        search_for_new(**self.kwords)

        faulty_resulting_fastq = os.path.join(self.bc_dir, "fastq/_111009_AB0CDDECXX_fastq.txt")
        assert not os.path.exists(faulty_resulting_fastq), \
        "Unexpected fastq was created"

    def test_search_for_new_fastq_single(self):
        self.kwords["fastq"] = True
        qseq_str = "0\t" * 8 + ("A" * 4 + "\t") * 2 + "1\t"
        with open(os.path.join(self.bc_dir, "s_3_1_1108_qseq.txt"), "w") as h:
            h.write(qseq_str)

        search_for_new(**self.kwords)

        resulting_fastq = os.path.join(self.bc_dir, "fastq/3_111009_AB0CDDECXX_fastq.txt")
        assert os.path.exists(resulting_fastq), \
        "Single read fastq was not created"

        given_qseq = os.path.join(self.bc_dir, "s_3_1_1108_qseq.txt")
        assert os.path.exists(given_qseq), \
        "Qseq file 1 was removed"

    def test_search_for_new_fastq_paired(self):
        self.kwords["fastq"] = True

        qseq_str = "0\t" * 8 + ("A" * 4 + "\t") * 2 + "1\t"
        with open(os.path.join(self.bc_dir, "s_3_1_1108_qseq.txt"), "w") as h:
            h.write(qseq_str)
        qseq_str = "0\t" * 8 + ("A" + "\t") * 2 + "1\t"
        with open(os.path.join(self.bc_dir, "s_3_2_1108_qseq.txt"), "w") as h:
            h.write(qseq_str)
        qseq_str = "0\t" * 8 + ("A" * 4 + "\t") * 2 + "1\t"
        with open(os.path.join(self.bc_dir, "s_3_3_1108_qseq.txt"), "w") as h:
            h.write(qseq_str)

        search_for_new(**self.kwords)

        resulting_fastq = os.path.join(self.bc_dir, "fastq/3_111009_AB0CDDECXX_1_fastq.txt")
        assert os.path.exists(resulting_fastq), \
        "Fastq file 1 not created"

        resulting_fastq = os.path.join(self.bc_dir, "fastq/3_111009_AB0CDDECXX_2_fastq.txt")
        assert os.path.exists(resulting_fastq), \
        "Fastq file 2 not created"

        given_qseq = os.path.join(self.bc_dir, "s_3_1_1108_qseq.txt")
        assert os.path.exists(given_qseq), \
        "Qseq file 1 was removed"

        given_qseq = os.path.join(self.bc_dir, "s_3_2_1108_qseq.txt")
        assert os.path.exists(given_qseq), \
        "Qseq file 2 was removed"

        given_qseq = os.path.join(self.bc_dir, "s_3_3_1108_qseq.txt")
        assert os.path.exists(given_qseq), \
        "Qseq file 3 was removed"

    def test_search_for_new_fastq_paired_clean_qseq(self):
        self.kwords["fastq"] = True
        self.kwords["remove_qseq"] = True

        qseq_str = "0\t" * 8 + ("A" * 4 + "\t") * 2 + "1\t"
        with open(os.path.join(self.bc_dir, "s_3_1_1108_qseq.txt"), "w") as h:
            h.write(qseq_str)
        qseq_str = "0\t" * 8 + ("A" + "\t") * 2 + "1\t"
        with open(os.path.join(self.bc_dir, "s_3_2_1108_qseq.txt"), "w") as h:
            h.write(qseq_str)
        qseq_str = "0\t" * 8 + ("A" * 4 + "\t") * 2 + "1\t"
        with open(os.path.join(self.bc_dir, "s_3_3_1108_qseq.txt"), "w") as h:
            h.write(qseq_str)

        search_for_new(**self.kwords)

        resulting_fastq = os.path.join(self.bc_dir, "fastq/3_111009_AB0CDDECXX_1_fastq.txt")
        assert os.path.exists(resulting_fastq), \
        "Fastq file 1 not created"

        resulting_fastq = os.path.join(self.bc_dir, "fastq/3_111009_AB0CDDECXX_2_fastq.txt")
        assert os.path.exists(resulting_fastq), \
        "Fastq file 2 not created"

        given_qseq = os.path.join(self.bc_dir, "s_3_1_1108_qseq.txt")
        assert not os.path.exists(given_qseq), \
        "Qseq file 1 was not removed"

        given_qseq = os.path.join(self.bc_dir, "s_3_2_1108_qseq.txt")
        assert not os.path.exists(given_qseq), \
        "Qseq file 2 was not removed"

        given_qseq = os.path.join(self.bc_dir, "s_3_3_1108_qseq.txt")
        assert not os.path.exists(given_qseq), \
        "Qseq file 3 was not removed"

    def test_search_for_new_fastq_paired_compress(self):
        self.kwords["fastq"] = True
        self.kwords["compress_fastq"] = True

        qseq_str = "0\t" * 8 + ("A" * 4 + "\t") * 2 + "1\t"
        with open(os.path.join(self.bc_dir, "s_3_1_1108_qseq.txt"), "w") as h:
            h.write(qseq_str)
        qseq_str = "0\t" * 8 + ("A" + "\t") * 2 + "1\t"
        with open(os.path.join(self.bc_dir, "s_3_2_1108_qseq.txt"), "w") as h:
            h.write(qseq_str)
        qseq_str = "0\t" * 8 + ("A" * 4 + "\t") * 2 + "1\t"
        with open(os.path.join(self.bc_dir, "s_3_3_1108_qseq.txt"), "w") as h:
            h.write(qseq_str)

        _compress_fastq.sleeptime = 1
        search_for_new(**self.kwords)

        resulting_fastq = os.path.join(self.bc_dir, "fastq/3_111009_AB0CDDECXX_1_fastq.txt")
        assert not os.path.exists(resulting_fastq), \
        "Fastq file 1 was not removed during compression"

        resulting_fastq = os.path.join(self.bc_dir, "fastq/3_111009_AB0CDDECXX_2_fastq.txt")
        assert not os.path.exists(resulting_fastq), \
        "Fastq file 2 was not removed during compression"

        resulting_fastq = os.path.join(self.bc_dir, "fastq/3_111009_AB0CDDECXX_1_fastq.txt.gz")
        assert os.path.exists(resulting_fastq), \
        "Fastq file 1 was not compressed"

        resulting_fastq = os.path.join(self.bc_dir, \
            "fastq/3_111009_AB0CDDECXX_2_fastq.txt.gz")
        assert os.path.exists(resulting_fastq), \
        "Fastq file 2 was not compressed"

        given_qseq = os.path.join(self.bc_dir, "s_3_1_1108_qseq.txt")
        assert os.path.exists(given_qseq), \
        "Qseq file 1 was removed"

        given_qseq = os.path.join(self.bc_dir, "s_3_2_1108_qseq.txt")
        assert os.path.exists(given_qseq), \
        "Qseq file 2 was removed"

        given_qseq = os.path.join(self.bc_dir, "s_3_3_1108_qseq.txt")
        assert os.path.exists(given_qseq), \
        "Qseq file 3 was removed"

    def tearDown(self):
        shutil.rmtree("test_data")

from mock import MagicMock


class CasavaTest(IFMTestCase):
    """Tests for using Casava 1.8 to make fastq files
    """
    def setUp(self):
        self.test_dir = "tst_" + "".join( \
            random.choice(string.ascii_uppercase) for i in xrange(3))
        self.bc_dir = os.path.join(self.test_dir, \
            "111009_SN1_0002_AB0CDDECXX/Data/Intensities/BaseCalls/")
        os.makedirs(self.bc_dir)

        self.msg_db = os.path.join(self.test_dir, "transferred.db")
        open(self.msg_db, "w").close()

        self.kwords = {
            "config": {
                "msg_db": self.msg_db,
                "dump_directories": self.test_dir,
                "algorithm": {},
                "program": {}
                },
            "config_file": None,
            "post_config_file": None,

            "fetch_msg": False,
            "process_msg": False,
            "store_msg": False,
            "backup_msg": False,

            "qseq": False,
            "fastq": False,
            "remove_qseq": False,
            "compress_fastq": False,
            "casava": True
            }

    def test_search_for_new_read_1_no_casava(self):
        open(os.path.join(self.test_dir, "111009_SN1_0002_AB0CDDECXX", \
            "Basecalling_Netcopy_complete_Read2.txt"), "w").close()

        with self.assertRaises(KeyError) as ke:
            search_for_new(**self.kwords)

        self.assertEqual(ke.exception.message, "casava")

    def test_search_for_new_read_1(self):
        self.kwords["config"]["program"]["casava"] = "casava_dir"

        subprocess.check_call = MagicMock()

        search_for_new(**self.kwords)

        self.assertFalse(subprocess.check_call.called)

        open(os.path.join(self.test_dir, "111009_SN1_0002_AB0CDDECXX", \
            "Basecalling_Netcopy_complete_Read2.txt"), "w").close()

        # Make the messaging.runner functional return a mock function
        runner = MagicMock()
        messaging.runner = MagicMock(side_effect=lambda *args, **kwargs: runner)

        search_for_new(**self.kwords)

        self.assertEqual(runner.call_args[0][0], "fetch_data")

    def test_search_for_new_read_2(self):
        self.kwords["config"]["program"]["casava"] = "casava_dir"

        open(os.path.join(self.test_dir, "111009_SN1_0002_AB0CDDECXX", \
            "Basecalling_Netcopy_complete_Read2.txt"), "w").close()
        open(os.path.join(self.test_dir, "111009_SN1_0002_AB0CDDECXX", \
            "Demultiplexing_done_Read1.txt"), "w").close()

        subprocess.check_call = MagicMock()
        search_for_new(**self.kwords)
        
        self.assertFalse(subprocess.check_call.called)


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-b", "--backup", dest="backup_msg",
            action="store_true", default=False)
    parser.add_option("-d", "--nofetch", dest="fetch_msg",
            action="store_false", default=True)
    parser.add_option("-p", "--noprocess", dest="process_msg",
            action="store_false", default=True)
    parser.add_option("-s", "--nostore", dest="store_msg",
            action="store_false", default=True)
    parser.add_option("-f", "--nofastq", dest="fastq",
            action="store_false", default=True)
    parser.add_option("-z", "--compress-fastq", dest="compress_fastq",
            action="store_true", default=False)
    parser.add_option("-q", "--noqseq", dest="qseq",
            action="store_false", default=True)
    parser.add_option("-c", "--casava", dest="casava",
            action="store_true", default=False)
    parser.add_option("-r", "--remove-qseq", dest="remove_qseq",
            action="store_true", default=False)
    parser.add_option("-m", "--miseq", dest="miseq",
            action="store_true", default=False)

    (options, args) = parser.parse_args()

    # Option --miseq implies --noprocess, --nostore, --nofastq, --noqseq
    if options.miseq:
        options.fetch_msg = False
        options.process_msg = False
        options.store_msg = False
        options.backup_msg = True
        options.fastq = False
        options.qseq = False
        options.casava = False

    kwargs = dict(fetch_msg=options.fetch_msg, process_msg=options.process_msg,
                  store_msg=options.store_msg, backup_msg=options.backup_msg,
                  fastq=options.fastq, qseq=options.qseq, remove_qseq=options.remove_qseq,
                  compress_fastq=options.compress_fastq, casava=options.casava)

    main(*args, **kwargs)
