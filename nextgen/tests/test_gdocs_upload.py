"""A test that will attempt to generate demultiplex reports and upload to Google Docs
"""
import os
import subprocess
import unittest
import shutil
import contextlib
import collections
import yaml
import random
from test_automated_analysis import make_workdir
from bcbio.solexa.flowcell import get_flowcell_info
from bcbio.utils import UnicodeWriter
from bcbio.google.sequencing_report import create_report_on_gdocs
from bcbio.pipeline.config_loader import load_config
import xml.etree.ElementTree

read_qc = ('<!--Illumina RTA Data--><Summary Read="1" densityRatio="0.173619792"><Lane key="1" TileCount="48" ClustersRaw="1860069" ClustersRawSD="361361.1" ClustersPF="1831969" ClustersPFSD="352563.3" PrcPFClusters="98.5" PrcPFClustersSD="0.22" Phasing="0.205" Prephasing="0.457" CalledCyclesMin="100" CalledCyclesMax="100" PrcAlign="1.28" PrcAlignSD="0.054" ErrRatePhiX="0.31" ErrRatePhiXSD="0.091" ErrRate35="0.10" ErrRate35SD="0.050" ErrRate75="0.17" ErrRate75SD="0.045" ErrRate100="0.31" ErrRate100SD="0.091" FirstCycleIntPF="6608" FirstCycleIntPFSD="229.6" PrcIntensityAfter20CyclesPF="79.1" PrcIntensityAfter20CyclesPFSD="1.51" /><Lane key="2" TileCount="48" ClustersRaw="1892872" ClustersRawSD="357670.3" ClustersPF="1863913" ClustersPFSD="348240.3" PrcPFClusters="98.5" PrcPFClustersSD="0.20" Phasing="0.197" Prephasing="0.462" CalledCyclesMin="100" CalledCyclesMax="100" PrcAlign="1.28" PrcAlignSD="0.066" ErrRatePhiX="0.44" ErrRatePhiXSD="0.272" ErrRate35="0.09" ErrRate35SD="0.008" ErrRate75="0.31" ErrRate75SD="0.353" ErrRate100="0.44" ErrRate100SD="0.272" FirstCycleIntPF="6577" FirstCycleIntPFSD="224.7" PrcIntensityAfter20CyclesPF="79.3" PrcIntensityAfter20CyclesPFSD="1.70" /><Lane key="3" TileCount="48" ClustersRaw="5882812" ClustersRawSD="859392.4" ClustersPF="4526246" ClustersPFSD="1773558.1" PrcPFClusters="77.1" PrcPFClustersSD="29.36" Phasing="0.249" Prephasing="0.396" CalledCyclesMin="0" CalledCyclesMax="100" PrcAlign="0.23" PrcAlignSD="0.106" ErrRatePhiX="0.53" ErrRatePhiXSD="0.228" ErrRate35="0.18" ErrRate35SD="0.173" ErrRate75="0.27" ErrRate75SD="0.102" ErrRate100="0.53" ErrRate100SD="0.228" FirstCycleIntPF="2600" FirstCycleIntPFSD="350.7" PrcIntensityAfter20CyclesPF="171.9" PrcIntensityAfter20CyclesPFSD="27.66" /><Lane key="4" TileCount="48" ClustersRaw="5139881" ClustersRawSD="867157.5" ClustersPF="4029611" ClustersPFSD="1747588.3" PrcPFClusters="79.2" PrcPFClustersSD="31.65" Phasing="0.249" Prephasing="0.395" CalledCyclesMin="0" CalledCyclesMax="100" PrcAlign="0.26" PrcAlignSD="0.120" ErrRatePhiX="0.47" ErrRatePhiXSD="0.222" ErrRate35="0.21" ErrRate35SD="0.265" ErrRate75="0.32" ErrRate75SD="0.220" ErrRate100="0.47" ErrRate100SD="0.222" FirstCycleIntPF="2342" FirstCycleIntPFSD="326.4" PrcIntensityAfter20CyclesPF="158.6" PrcIntensityAfter20CyclesPFSD="71.21" /><Lane key="5" TileCount="48" ClustersRaw="5744301" ClustersRawSD="538240.6" ClustersPF="4696004" ClustersPFSD="453621.1" PrcPFClusters="82.0" PrcPFClustersSD="6.67" Phasing="0.180" Prephasing="0.404" CalledCyclesMin="100" CalledCyclesMax="100" PrcAlign="0.28" PrcAlignSD="0.031" ErrRatePhiX="0.41" ErrRatePhiXSD="0.139" ErrRate35="0.18" ErrRate35SD="0.154" ErrRate75="0.26" ErrRate75SD="0.105" ErrRate100="0.41" ErrRate100SD="0.139" FirstCycleIntPF="5644" FirstCycleIntPFSD="160.2" PrcIntensityAfter20CyclesPF="72.9" PrcIntensityAfter20CyclesPFSD="2.56" /><Lane key="6" TileCount="48" ClustersRaw="6001645" ClustersRawSD="457375.5" ClustersPF="4866307" ClustersPFSD="293214.2" PrcPFClusters="81.2" PrcPFClustersSD="1.65" Phasing="0.175" Prephasing="0.391" CalledCyclesMin="100" CalledCyclesMax="100" PrcAlign="0.24" PrcAlignSD="0.026" ErrRatePhiX="0.46" ErrRatePhiXSD="0.101" ErrRate35="0.19" ErrRate35SD="0.090" ErrRate75="0.29" ErrRate75SD="0.074" ErrRate100="0.46" ErrRate100SD="0.101" FirstCycleIntPF="5655" FirstCycleIntPFSD="178.7" PrcIntensityAfter20CyclesPF="69.9" PrcIntensityAfter20CyclesPFSD="1.62" /><Lane key="7" TileCount="48" ClustersRaw="5006987" ClustersRawSD="977030.1" ClustersPF="24481" ClustersPFSD="167834.1" PrcPFClusters="0.4" PrcPFClustersSD="2.94" Phasing="0.211" Prephasing="0.000" CalledCyclesMin="0" CalledCyclesMax="0" PrcAlign="0.00" PrcAlignSD="0.002" ErrRatePhiX="0.00" ErrRatePhiXSD="0.000" ErrRate35="0.00" ErrRate35SD="0.000" ErrRate75="0.00" ErrRate75SD="0.000" ErrRate100="0.00" ErrRate100SD="0.000" FirstCycleIntPF="4714" FirstCycleIntPFSD="2287.2" PrcIntensityAfter20CyclesPF="8.9" PrcIntensityAfter20CyclesPFSD="20.87" /><Lane key="8" TileCount="48" ClustersRaw="6154610" ClustersRawSD="193792.1" ClustersPF="2011986" ClustersPFSD="1202392.4" PrcPFClusters="32.6" PrcPFClustersSD="19.34" Phasing="0.085" Prephasing="0.110" CalledCyclesMin="0" CalledCyclesMax="100" PrcAlign="0.08" PrcAlignSD="0.045" ErrRatePhiX="1.73" ErrRatePhiXSD="0.302" ErrRate35="1.69" ErrRate35SD="0.401" ErrRate75="1.48" ErrRate75SD="0.303" ErrRate100="1.73" ErrRate100SD="0.302" FirstCycleIntPF="5821" FirstCycleIntPFSD="153.4" PrcIntensityAfter20CyclesPF="62.2" PrcIntensityAfter20CyclesPFSD="9.29" /></Summary>',
           '<!--Illumina RTA Data--><Summary Read="2" ReadType=" (Index)" densityRatio="0.173619792"><Lane key="1" TileCount="48" ClustersRaw="1860069" ClustersRawSD="361361.1" ClustersPF="1831969" ClustersPFSD="352563.3" PrcPFClusters="98.5" PrcPFClustersSD="0.22" Phasing="0.000" Prephasing="0.000" CalledCyclesMin="0" CalledCyclesMax="0" PrcAlign="0.00" PrcAlignSD="0.000" ErrRatePhiX="0.00" ErrRatePhiXSD="0.000" ErrRate35="0.00" ErrRate35SD="0.000" ErrRate75="0.00" ErrRate75SD="0.000" ErrRate100="0.00" ErrRate100SD="0.000" FirstCycleIntPF="7045" FirstCycleIntPFSD="185.0" PrcIntensityAfter20CyclesPF="0.0" PrcIntensityAfter20CyclesPFSD="0.00" /><Lane key="2" TileCount="48" ClustersRaw="1892872" ClustersRawSD="357670.3" ClustersPF="1863913" ClustersPFSD="348240.3" PrcPFClusters="98.5" PrcPFClustersSD="0.20" Phasing="0.000" Prephasing="0.000" CalledCyclesMin="0" CalledCyclesMax="0" PrcAlign="0.00" PrcAlignSD="0.000" ErrRatePhiX="0.00" ErrRatePhiXSD="0.000" ErrRate35="0.00" ErrRate35SD="0.000" ErrRate75="0.00" ErrRate75SD="0.000" ErrRate100="0.00" ErrRate100SD="0.000" FirstCycleIntPF="6989" FirstCycleIntPFSD="230.2" PrcIntensityAfter20CyclesPF="0.0" PrcIntensityAfter20CyclesPFSD="0.00" /><Lane key="3" TileCount="48" ClustersRaw="5882812" ClustersRawSD="859392.4" ClustersPF="4526246" ClustersPFSD="1773558.1" PrcPFClusters="77.1" PrcPFClustersSD="29.36" Phasing="0.000" Prephasing="0.000" CalledCyclesMin="0" CalledCyclesMax="0" PrcAlign="0.00" PrcAlignSD="0.000" ErrRatePhiX="0.00" ErrRatePhiXSD="0.000" ErrRate35="0.00" ErrRate35SD="0.000" ErrRate75="0.00" ErrRate75SD="0.000" ErrRate100="0.00" ErrRate100SD="0.000" FirstCycleIntPF="2457" FirstCycleIntPFSD="198.5" PrcIntensityAfter20CyclesPF="0.0" PrcIntensityAfter20CyclesPFSD="0.00" /><Lane key="4" TileCount="48" ClustersRaw="5139881" ClustersRawSD="867157.5" ClustersPF="4029611" ClustersPFSD="1747588.3" PrcPFClusters="79.2" PrcPFClustersSD="31.65" Phasing="0.000" Prephasing="0.000" CalledCyclesMin="0" CalledCyclesMax="0" PrcAlign="0.00" PrcAlignSD="0.000" ErrRatePhiX="0.00" ErrRatePhiXSD="0.000" ErrRate35="0.00" ErrRate35SD="0.000" ErrRate75="0.00" ErrRate75SD="0.000" ErrRate100="0.00" ErrRate100SD="0.000" FirstCycleIntPF="7045" FirstCycleIntPFSD="386.9" PrcIntensityAfter20CyclesPF="0.0" PrcIntensityAfter20CyclesPFSD="0.00" /><Lane key="5" TileCount="48" ClustersRaw="5744301" ClustersRawSD="538240.6" ClustersPF="4696004" ClustersPFSD="453621.1" PrcPFClusters="82.0" PrcPFClustersSD="6.67" Phasing="0.000" Prephasing="0.000" CalledCyclesMin="0" CalledCyclesMax="0" PrcAlign="0.00" PrcAlignSD="0.000" ErrRatePhiX="0.00" ErrRatePhiXSD="0.000" ErrRate35="0.00" ErrRate35SD="0.000" ErrRate75="0.00" ErrRate75SD="0.000" ErrRate100="0.00" ErrRate100SD="0.000" FirstCycleIntPF="5827" FirstCycleIntPFSD="214.7" PrcIntensityAfter20CyclesPF="0.0" PrcIntensityAfter20CyclesPFSD="0.00" /><Lane key="6" TileCount="48" ClustersRaw="6001645" ClustersRawSD="457375.5" ClustersPF="4866307" ClustersPFSD="293214.2" PrcPFClusters="81.2" PrcPFClustersSD="1.65" Phasing="0.000" Prephasing="0.000" CalledCyclesMin="0" CalledCyclesMax="0" PrcAlign="0.00" PrcAlignSD="0.000" ErrRatePhiX="0.00" ErrRatePhiXSD="0.000" ErrRate35="0.00" ErrRate35SD="0.000" ErrRate75="0.00" ErrRate75SD="0.000" ErrRate100="0.00" ErrRate100SD="0.000" FirstCycleIntPF="5945" FirstCycleIntPFSD="230.8" PrcIntensityAfter20CyclesPF="0.0" PrcIntensityAfter20CyclesPFSD="0.00" /><Lane key="7" TileCount="48" ClustersRaw="5006987" ClustersRawSD="977030.1" ClustersPF="24481" ClustersPFSD="167834.1" PrcPFClusters="0.4" PrcPFClustersSD="2.94" Phasing="0.000" Prephasing="0.000" CalledCyclesMin="0" CalledCyclesMax="0" PrcAlign="0.00" PrcAlignSD="0.000" ErrRatePhiX="0.00" ErrRatePhiXSD="0.000" ErrRate35="0.00" ErrRate35SD="0.000" ErrRate75="0.00" ErrRate75SD="0.000" ErrRate100="0.00" ErrRate100SD="0.000" FirstCycleIntPF="13" FirstCycleIntPFSD="64.0" PrcIntensityAfter20CyclesPF="0.0" PrcIntensityAfter20CyclesPFSD="0.00" /><Lane key="8" TileCount="48" ClustersRaw="6154610" ClustersRawSD="193792.1" ClustersPF="2011986" ClustersPFSD="1202392.4" PrcPFClusters="32.6" PrcPFClustersSD="19.34" Phasing="0.000" Prephasing="0.000" CalledCyclesMin="0" CalledCyclesMax="0" PrcAlign="0.00" PrcAlignSD="0.000" ErrRatePhiX="0.00" ErrRatePhiXSD="0.000" ErrRate35="0.00" ErrRate35SD="0.000" ErrRate75="0.00" ErrRate75SD="0.000" ErrRate100="0.00" ErrRate100SD="0.000" FirstCycleIntPF="1111" FirstCycleIntPFSD="49.8" PrcIntensityAfter20CyclesPF="0.0" PrcIntensityAfter20CyclesPFSD="0.00" /></Summary>',
           '<!--Illumina RTA Data--><Summary Read="3" densityRatio="0.173619792"><Lane key="1" TileCount="48" ClustersRaw="1860069" ClustersRawSD="361361.1" ClustersPF="1831969" ClustersPFSD="352563.3" PrcPFClusters="98.5" PrcPFClustersSD="0.22" Phasing="0.190" Prephasing="0.557" CalledCyclesMin="208" CalledCyclesMax="208" PrcAlign="1.28" PrcAlignSD="0.056" ErrRatePhiX="0.44" ErrRatePhiXSD="0.160" ErrRate35="0.11" ErrRate35SD="0.039" ErrRate75="0.25" ErrRate75SD="0.115" ErrRate100="0.44" ErrRate100SD="0.160" FirstCycleIntPF="6093" FirstCycleIntPFSD="242.3" PrcIntensityAfter20CyclesPF="76.0" PrcIntensityAfter20CyclesPFSD="1.05" /><Lane key="2" TileCount="48" ClustersRaw="1892872" ClustersRawSD="357670.3" ClustersPF="1863913" ClustersPFSD="348240.3" PrcPFClusters="98.5" PrcPFClustersSD="0.20" Phasing="0.190" Prephasing="0.556" CalledCyclesMin="208" CalledCyclesMax="208" PrcAlign="1.27" PrcAlignSD="0.068" ErrRatePhiX="0.47" ErrRatePhiXSD="0.146" ErrRate35="0.11" ErrRate35SD="0.012" ErrRate75="0.24" ErrRate75SD="0.042" ErrRate100="0.47" ErrRate100SD="0.146" FirstCycleIntPF="6095" FirstCycleIntPFSD="260.7" PrcIntensityAfter20CyclesPF="75.5" PrcIntensityAfter20CyclesPFSD="1.21" /><Lane key="3" TileCount="48" ClustersRaw="5882812" ClustersRawSD="859392.4" ClustersPF="4526246" ClustersPFSD="1773558.1" PrcPFClusters="77.1" PrcPFClustersSD="29.36" Phasing="0.276" Prephasing="0.477" CalledCyclesMin="0" CalledCyclesMax="208" PrcAlign="0.25" PrcAlignSD="0.076" ErrRatePhiX="3.51" ErrRatePhiXSD="7.658" ErrRate35="0.63" ErrRate35SD="1.157" ErrRate75="2.35" ErrRate75SD="5.385" ErrRate100="3.51" ErrRate100SD="7.658" FirstCycleIntPF="1875" FirstCycleIntPFSD="693.7" PrcIntensityAfter20CyclesPF="340.7" PrcIntensityAfter20CyclesPFSD="378.91" /><Lane key="4" TileCount="48" ClustersRaw="5139881" ClustersRawSD="867157.5" ClustersPF="4029611" ClustersPFSD="1747588.3" PrcPFClusters="79.2" PrcPFClustersSD="31.65" Phasing="0.265" Prephasing="0.445" CalledCyclesMin="0" CalledCyclesMax="208" PrcAlign="0.27" PrcAlignSD="0.111" ErrRatePhiX="0.64" ErrRatePhiXSD="0.181" ErrRate35="0.24" ErrRate35SD="0.099" ErrRate75="0.40" ErrRate75SD="0.148" ErrRate100="0.64" ErrRate100SD="0.181" FirstCycleIntPF="2084" FirstCycleIntPFSD="82.7" PrcIntensityAfter20CyclesPF="192.7" PrcIntensityAfter20CyclesPFSD="13.38" /><Lane key="5" TileCount="48" ClustersRaw="5744301" ClustersRawSD="538240.6" ClustersPF="4696004" ClustersPFSD="453621.1" PrcPFClusters="82.0" PrcPFClustersSD="6.67" Phasing="0.188" Prephasing="0.478" CalledCyclesMin="208" CalledCyclesMax="208" PrcAlign="0.27" PrcAlignSD="0.027" ErrRatePhiX="0.92" ErrRatePhiXSD="0.212" ErrRate35="0.38" ErrRate35SD="0.276" ErrRate75="0.56" ErrRate75SD="0.164" ErrRate100="0.92" ErrRate100SD="0.212" FirstCycleIntPF="5332" FirstCycleIntPFSD="240.1" PrcIntensityAfter20CyclesPF="70.1" PrcIntensityAfter20CyclesPFSD="1.53" /><Lane key="6" TileCount="48" ClustersRaw="6001645" ClustersRawSD="457375.5" ClustersPF="4866307" ClustersPFSD="293214.2" PrcPFClusters="81.2" PrcPFClustersSD="1.65" Phasing="0.183" Prephasing="0.433" CalledCyclesMin="208" CalledCyclesMax="208" PrcAlign="0.23" PrcAlignSD="0.029" ErrRatePhiX="1.01" ErrRatePhiXSD="0.157" ErrRate35="0.36" ErrRate35SD="0.090" ErrRate75="0.63" ErrRate75SD="0.133" ErrRate100="1.01" ErrRate100SD="0.157" FirstCycleIntPF="5351" FirstCycleIntPFSD="234.9" PrcIntensityAfter20CyclesPF="67.5" PrcIntensityAfter20CyclesPFSD="0.72" /><Lane key="7" TileCount="48" ClustersRaw="5006987" ClustersRawSD="977030.1" ClustersPF="24481" ClustersPFSD="167834.1" PrcPFClusters="0.4" PrcPFClustersSD="2.94" Phasing="0.000" Prephasing="0.000" CalledCyclesMin="0" CalledCyclesMax="0" PrcAlign="0.00" PrcAlignSD="0.000" ErrRatePhiX="0.00" ErrRatePhiXSD="0.000" ErrRate35="0.00" ErrRate35SD="0.000" ErrRate75="0.00" ErrRate75SD="0.000" ErrRate100="0.00" ErrRate100SD="0.000" FirstCycleIntPF="4804" FirstCycleIntPFSD="1839.4" PrcIntensityAfter20CyclesPF="1.3" PrcIntensityAfter20CyclesPFSD="8.55" /><Lane key="8" TileCount="48" ClustersRaw="6154610" ClustersRawSD="193792.1" ClustersPF="2011986" ClustersPFSD="1202392.4" PrcPFClusters="32.6" PrcPFClustersSD="19.34" Phasing="0.442" Prephasing="0.931" CalledCyclesMin="0" CalledCyclesMax="208" PrcAlign="0.03" PrcAlignSD="0.040" ErrRatePhiX="4.10" ErrRatePhiXSD="0.540" ErrRate35="1.89" ErrRate35SD="0.216" ErrRate75="2.54" ErrRate75SD="0.606" ErrRate100="4.10" ErrRate100SD="0.540" FirstCycleIntPF="4737" FirstCycleIntPFSD="1966.5" PrcIntensityAfter20CyclesPF="47.3" PrcIntensityAfter20CyclesPFSD="25.19" /></Summary>')

class GDocsUploadTest(unittest.TestCase):
    """Setup a full automated analysis and run the pipeline.
    """
    def setUp(self):
        make_workdir()
        # Make up some barcode numbers
        self.workdir = os.path.join(os.path.dirname(__file__), "test_automated_output")
        self.data_dir = os.path.join(os.path.dirname(__file__), "data", "automated")
        
        # Parse the run_info
        run_info_file = os.path.join(self.data_dir, "run_info-gdocs.yaml")
        with open(run_info_file) as fh:
            self.run_info = yaml.load(fh)
            
        # Make up bogus run names
        self.runname = ("111014_SN0000_0001_AB0AAAACXX","111014_SN0000_0002_BB0AAAACXX")
        
        # Create the run directories (create them if necessary)
        for name in self.runname:
            analysisdir = os.path.join(self.workdir, name)
            if os.path.exists(analysisdir):
                shutil.rmtree(analysisdir)
            os.makedirs(analysisdir)
            self._make_bc_metrics(name,analysisdir)
            self._make_qc_metrics(name,analysisdir)
        

    def _make_qc_metrics(self, runname, analysisdir):
        """Writes RTA quality data for each read"""
        
        fc_name, fc_date = get_flowcell_info(runname)
        run_info_file = os.path.join(analysisdir,"RunInfo.xml")
        run_info_xml = "<RunInfo><Run Id=\"%s\" Number=\"%s\"><Flowcell>%s</Flowcell><Instrument>SN0000</Instrument><Date>%s</Date><Reads><Read Number=\"1\" NumCycles=\"101\" IsIndexedRead=\"N\" /><Read Number=\"2\" NumCycles=\"7\" IsIndexedRead=\"Y\" /><Read Number=\"3\" NumCycles=\"101\" IsIndexedRead=\"N\" /></Reads><FlowcellLayout LaneCount=\"8\" SurfaceCount=\"2\" SwathCount=\"3\" TileCount=\"8\" /><AlignToPhiX><Lane>1</Lane><Lane>2</Lane><Lane>3</Lane><Lane>4</Lane><Lane>5</Lane><Lane>6</Lane><Lane>7</Lane><Lane>8</Lane></AlignToPhiX></Run></RunInfo>" % (runname,1,fc_name,fc_date)
        xmlobj = xml.etree.ElementTree.fromstring(run_info_xml)
        xml.etree.ElementTree.ElementTree(xmlobj).write(run_info_file,"utf-8",True)
        
        qc_dir = os.path.join(analysisdir,"Data","reports","Summary")
        # Create the directory if it doesn't exist
        if not os.path.exists(qc_dir):      
            os.makedirs(qc_dir)
            
        for read in (1,2,3):
            xmlfile = os.path.join(qc_dir,"read%s.xml" % read)
            xmlobj = xml.etree.ElementTree.fromstring(read_qc[read-1])
            xml.etree.ElementTree.ElementTree(xmlobj).write(xmlfile,"utf-8",True)
            
    def _make_bc_metrics(self, runname, analysisdir):
        """Parses the run_info and generates lane folders and barcode metrics corresponding to the lanes and barcodes used"""
        fc_name, fc_date = get_flowcell_info(runname)
        barcode_dir_suffix = "_%s_%s_barcode" % (fc_date,fc_name)
        
        for lane in self.run_info:
            lane_name = str(lane['lane'])
            bc_dir = os.path.join(analysisdir,"%s%s" % (lane_name,barcode_dir_suffix))
            
            # Create the directory if it doesn't exist
            if not os.path.exists(bc_dir):      
                os.makedirs(bc_dir)
            
            # Create, or if it exists, append to the bc_metrics file
            bc_file = os.path.join(bc_dir,"%s_%s_%s_bc.metrics" % (lane_name,fc_date,fc_name))
            with open(bc_file,"a") as fh:
                bcw = UnicodeWriter(fh,dialect='excel-tab')
                
                # Loop over the barcodes and generate random read counts
                bcs = lane.get("multiplex",[])
                for bc in bcs:
                    bc_id = str(bc['barcode_id'])
                    bc_count = random.randint(1,10000000)
                    bcw.writerow([bc_id,bc_count])
                # Lastly write some unmatched counts, or in case no multiplex data was given, a 'trim' entry
                if len(bcs):
                    bcw.writerow(['unmatched',random.randint(1,10000000)])
                else:
                    bcw.writerow(['trim',random.randint(1,100000000)])

    def test_create_bc_report(self):
        """Create a demultiplex report and upload it to gdocs
        """
        # Parse the config
        config_file = os.path.join(self.data_dir, "post_process.yaml")
        self.config = load_config(config_file)

        # Loop over the runs
        for name in self.runname:
            print "\nProcessing %s" % name
            fc_name, fc_date = get_flowcell_info(name)
            analysisdir = os.path.join(self.workdir, name)
            assert create_report_on_gdocs(fc_date, fc_name, {'details': self.run_info}, {"work": analysisdir, "flowcell": analysisdir}, self.config), "Report creation failed"
