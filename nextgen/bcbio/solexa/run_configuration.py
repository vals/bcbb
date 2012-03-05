"""A helper class for parsing the illumina machine configuration for a run
"""

import os
import xml.etree.ElementTree as ET

class IlluminaConfiguration:

    def __init__(self, base_dir,config_xml="RunInfo.xml"):
        self._dir = base_dir
        self._run_info = os.path.join(base_dir,config_xml)
        assert os.path.exists(self._run_info), "The illumina configuration xml file %s does not exist" % self._run_info
        self._parse_configuration()
        
    def _parse_configuration(self):
        """Parse the run info xml configuration
        """
        tree = ET.parse(self._run_info)
        self._parse_run(tree.getroot())
    
    def _parse_run(self,elem):
        run = elem.find("Run")
        
        if run is not None:
            self._run_id = run.get("Id")
            self._run_number = run.get("Number")
            
            self._parse_flowcell(run)
            self._parse_instrument(run)
            self._parse_date(run)
            self._parse_reads(run)
            
    def _parse_flowcell(self,elem):
        node = elem.find("Flowcell")
        if node is not None:
            self._flowcell = node.text
        else:
            self._flowcell = None
        node = elem.find("FlowcellLayout")
        if node is not None:
            self._lanecount = node.get("LaneCount")
            self._surfacecount = node.get("SurfaceCount")
            self._swathcount = node.get("SwathCount")
            self._tilecount = node.get("TileCount")
        else:
            self._lanecount = None
            self._surfacecount = None
            self._swathcount = None
            self._tilecount = None
            
    def _parse_instrument(self,elem):
        node = elem.find("Instrument")
        if node is not None:
            self._instrument = node.text
        else:
            self._instrument = None
            
    def _parse_date(self,elem):
        node = elem.find("Date")
        if node is not None:
            self._date = node.text
        else:
            self._date = None
        
    def _parse_reads(self,elem):
        nodes = elem.findall("Reads/Read") or []
        self._readcount = len(nodes)
        self._reads = {}
        for read in nodes:
            num = read.get("Number")
            cycles = read.get("NumCycles")
            is_index = read.get("IsIndexedRead") == "Y"
            self._reads[num] = {'cycles': cycles, 'index': is_index}
        
    def reads(self):
        return self._reads
    def readcount(self):
        return self._readcount
    def indexread(self):
        indexread = []
        for num, read in self.reads().items():
            if read['index']:
                indexread.append(num)
        return indexread
    
    def date(self):
        return self._date
    def instrument(self):
        return self._instrument
    def flowcell(self):
        return self._flowcell
    def lanecount(self):
        return self._lanecount
    def tilecount(self):
        return self._tilecount
    def run_id(self):
        return self._run_id
    def run_number(self):
        return self._run_number
    
    def to_string(self):
        str = "Run id: %s\n\tStarted on %s, on instrument %s (run number %s)\n\tFlowcell: %s, %s lanes, %s tiles\n\t%s reads, read %s is an index read" % (self.run_id(),self.date(),self.instrument(),self.run_number(),self.flowcell(),self.lanecount(),self.tilecount(),self.readcount(),",".join(self.indexread()))
        for num, read in self.reads().items():
            str = ''.join((str,"\n\t\tRead %s is %s cycles" % (num,read['cycles'])))
        return str

        
    