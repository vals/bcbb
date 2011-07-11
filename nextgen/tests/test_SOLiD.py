"""
Test SOLiD related stuff
"""
import os
import sys
import unittest
import shutil

from bcbio.solid import *


class test_Targeted_Frag(unittest.TestCase):
    """Test the Targeted Frag pipeline"""
    def setUp(self):
        self.file_dir = os.path.join(os.path.dirname(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects", "j_doe_00_02")
        


class testSOLiDProjects(unittest.TestCase):
    def setUp(self):
        testdir = os.path.join(os.getcwd(), "projects", "j_doe_00_02")
        if os.path.exists(testdir):
            shutil.rmtree(testdir)
        #self.wtse = WT_SingleRead("Test", "Test", "ref", testdir, "csfasta", "filterref", "exons_gtf", "junction_ref", None)
        tfkw = {'runname':'test',
                'samplename':'testsample',
                'reference':'reference',
                'basedir':testdir,
                'targetfile':'target',
                'annotation_gtf_file':'annotation'
                }
        self.tf = TargetedFrag(**tfkw)

    # def testWT_SingleRead(self):
    #     print self.wtse.global_ini()
    #     print self.wtse.wt_single_read_ini()

    # def testTargetedFrag(self):
    #     print self.tf.global_ini()
    #     print self.tf.saet_ini()
    #     print self.tf.small_indel_frag_ini()
    #     print self.tf.enrichment_ini()
    #     print self.tf.targeted_frag_workflow_ini()

    def testTargetedFragPrimer(self):
        print self.tf.primersets['F3']
        self.tf.init_project()

