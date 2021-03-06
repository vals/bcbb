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

import logbook

from bcbio.solexa import samplesheet
from bcbio.log import create_log_handler, logger2
from bcbio import utils
from bcbio.distributed import messaging
from bcbio.solexa.flowcell import (get_flowcell_info, get_fastq_dir, get_qseq_dir)
from bcbio.pipeline.config_loader import load_config

LOG_NAME = os.path.splitext(os.path.basename(__file__))[0]
log = logbook.Logger(LOG_NAME)

def main(local_config, post_config_file=None, fetch_msg=True, process_msg=True, store_msg=True, 
         backup_msg=False, qseq=True, fastq=True, remove_qseq=False, compress_fastq=False, casava=False):
    config = load_config(local_config)
    log_handler = create_log_handler(config, True)

    with log_handler.applicationbound():
        search_for_new(config, local_config, post_config_file, fetch_msg, \
            process_msg, store_msg, backup_msg, qseq, fastq, remove_qseq, compress_fastq, casava)


def search_for_new(config, config_file, post_config_file, fetch_msg, \
        process_msg, store_msg, backup_msg, qseq, fastq, remove_qseq, compress_fastq, casava):
    """Search for any new unreported directories.
    """
    reported = _read_reported(config["msg_db"])
    for dname in _get_directories(config):
        if os.path.isdir(dname) and not any(dir.startswith(dname) for dir in reported):
            # Injects run_name on logging calls.
            # Convenient for run_name on "Subject" for email notifications
            run_setter = lambda record: record.extra.__setitem__('run', os.path.basename(dname))
            with logbook.Processor(run_setter):
                if casava and _is_finished_dumping_read_1(dname):
                    logger2.info("Generating fastq.gz files for read 1 of %s" % dname)
                    fastq_dir = None
                    _generate_fastq_with_casava(dname, config, r1=True)
                    _post_process_run(dname, config, config_file,
                                      fastq_dir, post_config_file,
                                      fetch_msg=True, process_msg=False,
                                      store_msg=store_msg,backup_msg=False)
    
                if _is_finished_dumping(dname):
                    logger2.info("The instrument has finished dumping on directory %s" % dname)
                    _update_reported(config["msg_db"], dname)
                    _process_samplesheets(dname, config)
                    if qseq:
                        logger2.info("Generating qseq files for %s" % dname)
                        _generate_qseq(get_qseq_dir(dname), config)

                    fastq_dir = None
                    if fastq:
                        logger2.info("Generating fastq files for %s" % dname)
                        fastq_dir = _generate_fastq(dname, config)
                        if remove_qseq: _clean_qseq(get_qseq_dir(dname), fastq_dir)
                        _calculate_md5(fastq_dir)
                        if compress_fastq: _compress_fastq(fastq_dir, config)
                    if casava:
                        logger2.info("Generating fastq.gz files for %s" % dname)
                        _generate_fastq_with_casava(dname, config)

                    _post_process_run(dname, config, config_file,
                                      fastq_dir, post_config_file,
                                      fetch_msg, process_msg, store_msg, backup_msg)

                    # Update the reported database after successful processing
                    _update_reported(config["msg_db"], dname)

                # Re-read the reported database to make sure it hasn't
                # changed while processing.
                reported = _read_reported(config["msg_db"])


def _post_process_run(dname, config, config_file, fastq_dir, post_config_file,
                      fetch_msg, process_msg, store_msg, backup_msg):
    """With a finished directory, send out message or process directly.
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
        logger2.info("CSV Samplesheet %s found, converting to %s" %
                 (ss_file, out_file))
        samplesheet.csv2yaml(ss_file, out_file)


def _generate_fastq_with_casava(fc_dir, config, r1=False):
    """Perform demultiplexing and generate fastq.gz files for the current
    flowecell using CASAVA (>1.8).
    """
    basecall_dir = os.path.join(fc_dir, "Data", "Intensities", "BaseCalls")
    casava_dir = config["program"].get("casava")
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
    short_fc_name = "%s_%s" % (fc_date, fc_name)
    fastq_dir = get_fastq_dir(fc_dir)
    basecall_dir = os.path.split(fastq_dir)[0]
    postprocess_dir = config.get("postprocess_dir", "")
    if postprocess_dir:
        fastq_dir = os.path.join(postprocess_dir, os.path.basename(fc_dir), "fastq")

    if not fastq_dir == fc_dir:# and not os.path.exists(fastq_dir):

        with utils.chdir(basecall_dir):
            lanes = sorted(list(set([f.split("_")[1] for f in
                glob.glob("*qseq.txt")])))
            cl = ["solexa_qseq_to_fastq.py", short_fc_name,
                  ",".join(lanes)]
            if postprocess_dir:
                cl += ["-o", fastq_dir]

            logger2.debug("Converting qseq to fastq on all lanes.")
            subprocess.check_call(cl)

    return fastq_dir

def _calculate_md5(fastq_dir):
    """Calculate the md5sum for the fastq files
    """
    glob_str = "*_fastq.txt"
    fastq_files = glob.glob(os.path.join(fastq_dir,glob_str))
    
    md5sum_file = os.path.join(fastq_dir,"md5sums.txt")
    with open(md5sum_file,'w') as fh:
        for fastq_file in fastq_files:
            logger2.debug("Calculating md5 for %s using md5sum" % fastq_file)
            cl = ["md5sum",fastq_file]
            fh.write(subprocess.check_output(cl))

def _compress_fastq(fastq_dir, config):
    """Compress the fastq files using gzip
    """
    glob_str = "*_fastq.txt"
    fastq_files = glob.glob(os.path.join(fastq_dir,glob_str))
    num_cores = config["algorithm"].get("num_cores",1)
    active_procs = []
    for fastq_file in fastq_files:
        # Sleep for one minute while waiting for an open slot
        while len(active_procs) >= num_cores:
            time.sleep(60)
            active_procs, _ = _process_status(active_procs)
            
        logger2.debug("Compressing %s using gzip" % fastq_file)
        cl = ["gzip",fastq_file]
        active_procs.append(subprocess.Popen(cl))
    
    # Wait for the last processes to finish
    while len(active_procs) > 0:
        time.sleep(60)
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
    fastq_files = glob.glob(os.path.join(fastq_dir,glob_str))
    
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
    hi_seq_checkpoint = "Basecalling_Netcopy_complete_Read%s.txt" % \
                        _expected_reads(run_info)
    # include a check to wait for any ongoing MiSeq analysis
    miseq_analysis_checkpoint = (not os.path.exists(os.path.join(directory, "QueuedForAnalysis.txt")) or os.path.exists(os.path.join(directory, "CompletedJobInfo.xml")))
    
    to_check = ["Basecalling_Netcopy_complete_SINGLEREAD.txt",
                "Basecalling_Netcopy_complete_READ2.txt",
                hi_seq_checkpoint]
    return (reduce(operator.or_,
            [os.path.exists(os.path.join(directory, f)) for f in to_check]) and miseq_analysis_checkpoint)


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
    if os.path.exists(msg_db):
        with open(msg_db) as in_handle:
            for line in in_handle:
                reported.append(line.strip())
    return reported


def _get_directories(config):
    for directory in config["dump_directories"]:
        globs = []
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
    for d in [dir for dir in reported if dir.startswith(new_dname)]:
        new_dname = d
        reported.remove(d)
    reported.append("%s\t%s" % (new_dname, time.strftime("%x-%X")))

    with open(msg_db, "w") as out_handle:
        for dir in reported:
            out_handle.write("%s\n" % dir)


def finished_message(fn_name, run_module, directory, files_to_copy,
                     config, config_file):
    """Wait for messages with the give tag, passing on to the supplied handler.
    """
    logger2.debug("Calling remote function: %s" % fn_name)
    user = getpass.getuser()
    hostname = socket.gethostbyaddr(socket.gethostname())[0]
    data = dict(
            machine_type='illumina',
            hostname=hostname,
            user=user,
            directory=directory,
            to_copy=files_to_copy
            )
    dirs = {"work": os.getcwd(),
            "config": os.path.dirname(config_file)}
    runner = messaging.runner(run_module, dirs, config, config_file, wait=False)
    runner(fn_name, [[data]])

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
                  store_msg=options.store_msg, backup_msg=options.backup_msg, fastq=options.fastq,
                  qseq=options.qseq, remove_qseq=options.remove_qseq, compress_fastq=options.compress_fastq, casava=options.casava)

    main(*args, **kwargs)
