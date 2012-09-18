
# Testing code: run with 'nosetests -v -s illumina_finished_message.py'
import unittest
import random
import string
import shutil
import os
import time
import subprocess
import errno

from bcbio.distributed import messaging


# import sys
# sys.path.insert(0, os.path.abspath('..'))
# import illumina_finished_msg

# from ..scripts.illumina_finished_message import _compress_fastq

from illumina_finished_msg import _compress_fastq
from illumina_finished_msg import _is_finished_dumping
from illumina_finished_msg import _read_reported
from illumina_finished_msg import _update_reported
from illumina_finished_msg import search_for_new

# import illumina_finished_message
# _compress_fastq = illumina_finished_message._compress_fastq


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

        # # Make the messaging.runner functional return a mock function
        # runner = MagicMock()
        # messaging.runner = MagicMock(side_effect=lambda *args, **kwargs: runner)

        # search_for_new(**self.kwords)

        # self.assertEqual(runner.call_args[0][0], "fetch_data")
