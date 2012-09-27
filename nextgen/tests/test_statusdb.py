"""Gather output metrics and upload to statusdb"""

import os
import unittest
import shutil
import contextlib
import subprocess
import couchdb

from bcbio.qc import FlowcellQCMetrics, QCMetrics, SampleQCMetrics, LaneQCMetrics
from bcbio.log import logger, setup_logging, version
from bcbio.pipeline.config_loader import load_config


data_dir = os.path.join(os.path.dirname(__file__), "test_automated_output")
sample_kw = dict(path=os.path.abspath(data_dir), flowcell="FC70BUKAAXX", date="110106", lane="1", barcode_name="test", barcode_id="1", sample_prj=None)
lane_kw = dict(path=os.path.abspath(data_dir), flowcell="FC70BUKAAXX", date="110106", lane="1")

data_dir = os.path.join(os.path.dirname(__file__), "data", "automated")
runinfo_yaml = os.path.join(data_dir, "run_info.yaml")
fc_kw = dict(path=os.path.abspath(data_dir), fc_date="110106", fc_name="FC70BUKAAXX", run_info_yaml=runinfo_yaml)

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
        config_file = os.path.join(self.data_dir, "post_process-statusdb.yaml")
        config = load_config(config_file)
        setup_logging(config)
        fc_date = "110106"
        fc_name = "FC70BUKAAXX"
        run_info_yaml = os.path.join(self.data_dir, "run_info.yaml")
        workdir = os.path.join(os.path.dirname(__file__), "110106_FC70BUKAAXX")
        fc_dir = os.path.join(self.data_dir, os.pardir, "110106_FC70BUKAAXX")

    def test_1_objects(self):
        obj = QCMetrics()
        obj = SampleQCMetrics(**sample_kw)
        obj.read_picard_metrics()
        obj.parse_fastq_screen()
        obj.read_fastqc_metrics()
        print obj 
        obj = LaneQCMetrics(**lane_kw)
        obj.parse_filter_metrics()
        obj.parse_bc_metrics()
        print obj
        obj = FlowcellQCMetrics(**fc_kw)
        print obj

