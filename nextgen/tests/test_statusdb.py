"""Gather output metrics and upload to statusdb"""

import os
import unittest
import shutil
import contextlib
import subprocess

@contextlib.contextmanager
def make_workdir(link=True):
    dirname = os.path.join(os.path.dirname(__file__), "110106_FC70BUKAAXX")
    if os.path.exists(dirname):
        if os.path.islink(dirname):
            os.remove(dirname)
        else:
            shutil.rmtree(dirname)
    src = os.path.join(os.path.dirname(__file__), "test_automated_output")
    orig_dir = os.getcwd()
    if link:
        os.symlink(src, dirname)
    else:
        shutil.copytree(src, dirname)
    try:
        os.chdir(dirname)
        yield
    finally:
        os.chdir(orig_dir)

class StatusDBTest(unittest.TestCase):
    """Setup test case for statusdb"""
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data", "automated")

    def test_1_statusdb(self):
        """Run full automated analysis pipeline with multiplexing, adding statusdb
        """
        self.setUp()
        config_file = os.path.join(self.data_dir, "post_process-statusdb.yaml")

        with make_workdir():
            cl = ["automated_initial_analysis.py",
                  config_file,
                  os.path.join(self.data_dir, os.pardir, "110106_FC70BUKAAXX"),
                  os.path.join(self.data_dir, "run_info.yaml")]
            subprocess.check_call(cl)

