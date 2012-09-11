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

        

    # def test_1_statusdb(self):
    #     """Run full automated analysis pipeline with multiplexing, adding statusdb
    #     """
    #     self.setUp()
    #     config_file = os.path.join(self.data_dir, "post_process-statusdb.yaml")

    #     with make_workdir():
    #         cl = ["automated_initial_analysis.py",
    #               config_file,
    #               os.path.join(self.data_dir, os.pardir, "110106_FC70BUKAAXX"),
    #               os.path.join(self.data_dir, "run_info.yaml")]
    #         subprocess.check_call(cl)

    # def test_2_object_equality(self):
    #     def _couchdoc_to_dict(obj):
    #         d = {}
    #         for i in obj.iterkeys():
    #             d[i] = obj[i]
    #         return d
    #     self.setUp()
    #     from bcbio.qc import FlowcellQCMetrics
    #     from bcbio.log import logger, setup_logging, version
    #     from bcbio.pipeline.config_loader import load_config
    #     config_file = os.path.join(self.data_dir, "post_process-statusdb.yaml")
    #     config = load_config(config_file)
    #     setup_logging(config)
    #     fc_date = "110106"
    #     fc_name = "FC70BUKAAXX"
    #     run_info_yaml = os.path.join(self.data_dir, "run_info.yaml")
    #     workdir = os.path.join(os.path.dirname(__file__), "110106_FC70BUKAAXX")
    #     fc_dir = os.path.join(self.data_dir, os.pardir, "110106_FC70BUKAAXX")
    #     qc_obj = FlowcellQCMetrics(fc_date, fc_name, run_info_yaml, workdir, fc_dir)
    #     couch = couchdb.Server(url="http://%s" % statusdb_url)
    #     db=couch['qc']
    #     dbobj = db.get(qc_obj.get_db_id())
    #     qc_obj["creation_time"] = dbobj["creation_time"]
    #     qc_obj["modification_time"] = dbobj["modification_time"]
    #     if dbobj == qc_obj:
    #         print "Identical document and dict"
    #     else:
    #         print "Different document and dict"
    #     dbobj_dict = _couchdoc_to_dict(dbobj)
    #     print dbobj_dict
    #     print dbobj_dict["creation_time"]
    #     print qc_obj
    #     if dbobj_dict == qc_obj:
    #         print "Identical dbobj_dict"
    #     else:
    #         print "Different dbdict and dict"

        #print qc_obj.__class__
        #print dir(dbobj)
