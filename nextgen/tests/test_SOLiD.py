"""
Test SOLiD related stuff

TODO: I'm using Applied's test data. The examples have been performed on chr20, a 61M reference file, not suitably 
small for distribution. The reference files should be downloaded and installed by the test script. Currently, the
reference I'm using resides at SciLife.

TODO: ditto for the target files
"""
import os
import sys
import unittest
import shutil

from bcbio.solid.solid import *

kw = {

    'reference' : "/bubo/home/h1/perun/opt/bcbb.git/nextgen/tests/data/solid/reference/chr20.fasta",
    'targetfile' : "/bubo/proj/b2010037/private/projects/perun/RARS.exome.kalkyl.git/data/SureSelect_All_Exon_G3362.bed",
    'saettargetfile' : "/bubo/proj/b2010037/private/projects/perun/RARS.exome.kalkyl.git/data/SureSelect_All_Exon_G3362.saet.bed",
    'cmap': os.path.join(os.path.dirname(__file__), "data", "solid", "reference", "GRCh37.dna.cmap")
    }
kw.update(
    runname = 'test_SOLiD',
    samplename = '20090114_2_HuRef_LMP_4CP1',
    basedir = os.path.join(os.path.dirname(__file__), "projects", "j_doe_00_02"),
    file_base = '20090114_2_HuRef_LMP_4CP1',
    )

# class test_Targeted_Frag(unittest.TestCase):
#     """Test the Targeted Frag pipeline"""
#     def setUp(self):
#         self.file_dir = os.path.join(os.path.dirname(__file__))
#         self.proj_dir = os.path.join(self.file_dir, "projects", "j_doe_00_02")
#         if os.path.exists(self.proj_dir):
#             shutil.rmtree(self.proj_dir)
#         self.sample = TargetedFrag(**kw)

#     def testTargetedFrag(self):
#         """Test installation of targeted frag config files"""
#         self.sample.init_project()

class test_Targeted_PE(unittest.TestCase):
    """Test the Targeted PE pipeline"""
    def setUp(self):
        self.file_dir = os.path.join(os.path.dirname(__file__))
        self.proj_dir = os.path.join(self.file_dir, "projects", "j_doe_00_03")
        if os.path.exists(self.proj_dir):
            shutil.rmtree(self.proj_dir)
        kw.update(
            runname = "test_Targeted_PE",
            basedir = os.path.join(os.path.dirname(__file__), "projects", "j_doe_00_03")
            )
        self.sample = TargetedPE(**kw)
    def testTargetedPE(self):
        """Test installation of targeted PE config files"""
        self.sample.init_project()

