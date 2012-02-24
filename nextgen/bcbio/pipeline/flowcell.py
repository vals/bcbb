"""A module for managing metadata and layout, as specified in run_info.yaml, as well as read counts for the samples on a flowcell"""

import copy
import re
import os
import glob

from bcbio.pipeline.run_info import get_run_info
from bcbio.google import ( _to_unicode, _from_unicode )
from bcbio.utils import UnicodeReader

def format_project_name(unformated_name):
    """Make the project name adhere to a stricter formatting convention"""
    regexp = r'^(.+?)_(\d{2})_(\d{2})(.*)$'
    m = re.match(regexp,unformated_name)
    if not m or len(m.groups()) < 3:
        return unformated_name
    
    name = m.group(1).strip()
    year = m.group(2).strip()
    month = m.group(3).strip()
    suffix = m.group(4).strip()
   
    # Replace any non-period delimiters
    delimiter = "_"
    p = re.compile('(_)')
    name = p.sub('.',name)
    
    # Format the name
    project_name = "%s_%s_%s%s" % (name,year,month,suffix)
    return project_name

def get_barcode_metrics(workdir):
    """Parse the [lane]_*_bc.metrics files in the *_barcode directories into a dictionary"""
    
    bc_files = []
    if workdir is not None:
        bc_files = glob.glob(os.path.join(workdir,"*_barcode","*_bc.metrics"))
    if not len(bc_files) > 0:
        return None
    
    bc_metrics = {}
    for bc_file in bc_files:
        m = re.match(r'^(\d+)\_',os.path.basename(bc_file))
        if not m or len(m.groups()) != 1:
            continue
        lane = str(m.group(1))
        bc_metrics[lane] = {}
        with open(bc_file) as bcfh:
            csvr = UnicodeReader(bcfh,dialect='excel-tab')
            for row in csvr:
                bc_metrics[lane][str(row[0])] = int(row[1])
            
    return bc_metrics

def get_flowcell(fc_dir, run_info_yaml, config={}):
    # Just get the name of the flowcell directory minus the path
    dirname = os.path.basename(os.path.normpath(fc_dir))
    fc_name, fc_date, run_info = get_run_info(dirname,config,run_info_yaml)
    return Flowcell(fc_name,fc_date,run_info.get("details",[]),fc_dir)

def get_project_name(description):
    """Parse out the project name from the lane description"""
    m = re.match(r'(?:.*\s+)?(\S+)',(description or ""),re.I)
    if m and len(m.groups()) > 0:
        return format_project_name(m.group(1).strip())
    return description
       
def get_sample_name(barcode_name):
    """Extract the sample name by stripping the barcode index part of the sample description""" 
    regexp = r'^(.+?)[\.\-_]?ind?(?:ex)?[\.\-_]?\d+$'
    m = re.search(regexp,(barcode_name or ""),re.I)
    if not m or len(m.groups()) == 0:
        return barcode_name
    return m.group(1)
    
class Flowcell:
    """A class for managing information about a flowcell"""
    
    def __init__(self,fc_name,fc_date,data,fc_dir=None):
        self.set_fc_dir(fc_dir)
        self.set_fc_date(fc_date)
        self.set_fc_name(fc_name)
        self.set_lanes(data)
        # Attempts to set the read counts on creation
        self.set_read_counts()
                
    def get_fc_date(self):
        return self.fc_date
    def set_fc_date(self,fc_date):
        self.fc_date = fc_date
        
    def get_fc_dir(self):
        return self.fc_dir
    def set_fc_dir(self,fc_dir):
        self.fc_dir = fc_dir
        
    def get_fc_name(self):
        return self.fc_name
    def set_fc_name(self,fc_name):
        self.fc_name = fc_name
        
    def get_lane_by_name(self,name):
        for lane in self.get_lanes():
            if (lane.get_name() == name):
                return lane
        return None
    
    def add_lane(self,lane):
        self.lanes.append(lane)
    def get_lanes(self):
        return self.lanes
    def set_lanes(self,lanes):
        self.lanes = []
        for lane in lanes:
            self.add_lane(Lane(lane))
        
    def get_project_names(self):
        pnames = {}
        for lane in self.get_lanes():
            for pname in lane.get_project_names():
                if pname is None:
                    continue
                pnames[pname] = 1
        return pnames.keys()
    
    def get_samples(self):
        samples = []
        for lane in self.get_lanes():
            samples.extend(lane.get_samples())
        return samples
            
    def prune_to_project(self,project):
        """Return a new Flowcell object just containing the lanes and samples belonging to a specific project"""
        lanes = []
        fc = None
        for lane in self.get_lanes():
            l = lane.prune_to_project(project)
            if (l):
                lanes.append(l.to_structure())
        if (len(lanes)):
            fc = Flowcell(self.get_fc_name(),self.get_fc_date(),lanes)
        return fc
    
    def set_read_counts(self,read_counts=None):
        """Sets the read counts of the barcoded samples in the lanes of this flowcell.
           Read counts can be supplied in a dictionary with lane number as key to a
           dictionary with barcode indexes (or 'unmatched' or 'trim') as key and counts 
           as values. If no read counts are supplied, attempts to parse the read counts 
           from the flowcell directory, assuming that the read counts can be found using a glob
           like [lane]_*_barcode/[lane]_*_bc.metrics
        """
        if read_counts is None:
            read_counts = get_barcode_metrics(self.get_fc_dir()) or {}
            
        for name in read_counts.keys():
            lane = self.get_lane_by_name(name)
            if not lane:
                lane = Lane({"lane": name, "description": "Unexpected lane"})
                self.add_lane(lane)
            lane.set_read_counts(read_counts[name])
    
    @staticmethod
    def columns():
        cols = []#"fc_date","fc_name"]
        cols.extend(Lane.columns())
        return cols
    
    def to_rows(self):
        rows = []
        for lane in self.get_lanes():
            for row in lane.to_rows():
                #l = [self.get_fc_date(),self.get_fc_name()]
                #l.extend(row)
                rows.append(row)
        return rows
    
    def to_structure(self):
        struct = []
        for lane in self.get_lanes():
            struct.append(lane.to_structure())
        return {"details": struct}
        
class Lane:
    """A class for managing information about a lane"""
    
    def __init__(self,data):
        self.set_data(data)
        self.set_description(data.get("description",None))
        self.set_name(data.get("lane",None))
        self.set_samples(data)
        self.set_files([])

    def get_project_names(self):
        pnames = {}
        for sample in self.get_samples():
            pnames[sample.get_project()] = 1
        return pnames.keys()
        
    def set_read_counts(self,read_counts):
        for barcode_id in read_counts.keys():
            sample = self.get_sample_by_barcode(barcode_id)
            if (sample is None):
                sample = BarcodedSample({"name": "Unexpected barcode", "barcode_id": barcode_id})
                self.add_sample(sample)
            sample.set_read_count(read_counts[barcode_id])
                
    def get_sample_by_barcode(self,barcode_id):
        for sample in self.get_samples():
            if (str(sample.get_barcode_id()) == str(barcode_id)):
                return sample
        return None
    
    def get_samples_by_project(self,project):
        samples = []
        for sample in self.get_samples():
            if (sample.get_project() == project):
                samples.append(sample)
        return samples
    
    def prune_to_project(self,project,exclude_unmatched=False):
        """Return a new Lane object just containing the samples belonging to a specific project"""
        samples = []
        lane = None
        for sample in self.get_samples_by_project(project):
            samples.append(sample.to_structure())
        if (len(samples)):
            struct = self.to_structure()
            
            # Add the unmatched samples unless specifically asked not to
            if (not exclude_unmatched):
                sample = self.get_sample_by_barcode("unmatched")
                if (sample):
                    samples.append(sample.to_structure())
                    
            struct["multiplex"] = samples
            struct["description"] = "%s-pruned" % project
            lane = Lane(struct)
        return lane
    
    @staticmethod
    def columns():
        cols = ["lane","description"]
        cols.extend(BarcodedSample.columns())
        return cols
    
    def to_rows(self):
        rows = []
        for sample in self.get_samples():
            s = sample.to_rows()
            s.insert(0,self.get_name())
            s.insert(1,self.get_description())
            rows.append(s)
        return rows
            
    def to_structure(self):
        struct = {}
        if (self.get_description()):
            struct["description"] = self.get_description()
        if (self.get_name()):
            struct["lane"] = self.get_name()
        if (self.get_files()):
            struct["files"] = self.get_files()
        if (self.get_samples()):
            struct["multiplex"] = []
            for sample in self.get_samples():
                struct["multiplex"].append(sample.to_structure())
        else:
            if (self.get_analysis()):
                struct["analysis"] = self.get_analysis()
            if (self.get_genome_build()):
                struct["genome_build"] = self.get_genome_build()
        return struct

    def get_barcode_ids(self):
        bcids = None
        if (self.get_samples()):
            bcids = []
            for sample in self.get_samples():
                bcids.append(sample.get_barcode_id())
        return bcids

    def __str__(self):
        s = "Lane: %s\n\nbarcode ids: %s" % (self.get_name(), self.get_barcode_ids())
        return s
        
    def set_data(self,data):
        self.data = copy.deepcopy(data)

    def _get_analysis(self):
        return self._analysis
    def _set_analysis(self,analysis):
        self._analysis = analysis

    def _get_genome_build(self):
        return self._genome_build
    def _set_genome_build(self,genome_build):
        self._genome_build = genome_build

    def _get_description(self):
        return _from_unicode(self._description)
    def _set_description(self,description):
        self._description = _to_unicode(description)
    
    def _get_name(self):
        return _from_unicode(self._name)
    def set_name(self,name):
        self._name = _to_unicode(str(name))
 
    def add_sample(self,sample):
        self._multiplex.append(sample)
    def _get_samples(self):
        return self._multiplex
    def _set_samples(self,multiplex):
        self._multiplex = []
        for barcode in multiplex:
            self.add_sample(BarcodedSample(barcode,self))

    def _set_files(self, files):
        self._files = files
    def _get_files(self):
        return self._files
       
    analysis = property(_get_analysis,_set_analysis)
    genome_build = property(_get_genome_build,_set_genome_build)
    description = property(_get_description,_set_description)
    name = property(_get_name,_set_name)
    sample = property(_get_sample,_set_sample)

class Sample:
    """A class for managing information about a sample"""
     
    def __init__(self,data,lane=Lane({}),comment=None):
        
        for (yaml,attribute) in self.yaml2attribute():
            setattr(self,attribute,data.get(yaml,getattr(lane,attribute)))
        setattr(self,"lane",lane.name)
        setattr(self,"comment",comment)
        
    def add_sample(self,other,delim=', '):
        if self.analysis != other.analysis:
            self.analysis("%s%s%s" % (self.analysis,delim,other.analysis))
        if self.genome_build != other.genome_build:
            self.genome_build("%s%s%s" % (self.genome_build,delim,other.genome_build))
        if self.sample_name != other.sample_name:
            self.sample_name("%s%s%s" % (self.sample_name,delim,other.sample_name))
        if self.project_name != other.project_name:
            self.project_name("%s%s%s" % (self.project_name,delim,other.project_name))
        if self.lane != other.lane:
            self.lane("%s%s%s" % (self.lane,delim,other.lane))
        if other.read_count is not None:
            self.read_count((self.read_count or 0) + other.read_count)
  
    @staticmethod
    def yaml2attribute():
        return [["description", "project_name"],
                ["name","sample_name"],
                #["sample_prj", "project_name"],
                ["read_count", "read_count"],
                ["rounded_read_count", "rounded_read_count"],
                ["analysis", "analysis"],
                ["genome_build", "genome_build"]]
          
    @staticmethod
    def columns():
        return [attribute for [yaml, attribute] in self.yaml2attribute()]
        
    def to_rows(self):
        return [getattr(self,attribute,None) for [yaml, attribute] in self.yaml2attribute()]
        
    def to_structure(self):
        struct = {}
        for (yaml,attribute) in self.yaml2attribute():
            val = getattr(self,attribute)
            if val is not None:
                struct[yaml] = val
        return struct
    
    def _get_analysis(self):
        return _from_unicode(self._analysis)
    def _set_analysis(self,analysis):
        self._analysis = _to_unicode(analysis)
        
    def _get_genome_build(self):
        return _from_unicode(self._genome_build)
    def _set_genome_build(self,genome_build):
        self._genome_build = _to_unicode(genome_build)
         
    def _get_sample_name(self):
        return _from_unicode(self._sample_name)
    def _set_sample_name(self,sample_name):
        self._sample_name = get_sample_name(_to_unicode(sample_name))
        
    def _get_comment(self):
        return self._comment
    def _set_comment(self,comment):
        self._comment = comment
        
    def _get_lane(self):
        return self._lane
    def _set_lane(self,lane):
        self._lane = lane
          
    def _get_read_count(self):
        if self._read_count:
            try:
                rc = int(self._read_count)
                return rc
            except ValueError:
                pass
        return None
    def _set_read_count(self,read_count):
        self._read_count = read_count

    def _get_rounded_read_count(self,unit=1000000,decimals=2):
        return round((self.read_count or 0)/float(unit),int(decimals))
    def _set_rounded_read_count(self,rounded_read_count):
        # Do nothing, we will always calculate this on the fly to ensure consistency
        pass
    
    def _get_sample_project(self):
        return _from_unicode(self._sample_project)
    def _set_sample_project(self,sample_project):
        self._sample_project = get_project_name(_to_unicode(sample_project))
        
    analysis = property(_get_analysis,_set_analysis)
    genome_build = property(_get_genome_build,_set_genome_build)
    sample_name = property(_get_sample_name,_set_sample_name)
    comment = property(_get_comment,_set_comment)
    lane = property(_get_lane,_set_lane)
    read_count = property(_get_read_count,_set_read_count)
    rounded_read_count = property(_get_rounded_read_count,_set_rounded_read_count)
    sample_project = property(_get_sample_project,_set_sample_project)
    
class BarcodedSample(Sample):
    """A subclass of Sample for managing information about a barcoded sample"""
    
    def __init__(self,data,lane):
        Sample.__init__(self,data,lane)
        
        for (yaml,attribute) in self.yaml2attribute():
            setattr(self,attribute,data.get(yaml,None))

    @staticmethod
    def yaml2attribute():
        return [["barcode_id","barcode_id"],
                ["name", "barcode_name"],
                ["sequence", "barcode_sequence"],
                ["barcode_type", "barcode_type"]]
    
    @staticmethod
    def columns():
        cols = Sample.columns()
        cols.extend([attribute for [yaml, attribute] in self.yaml2attribute()])
        return cols
        
    def to_rows(self):
        rows = Sample.to_rows(self)
        rows.extend([getattr(self,attribute,None) for [yaml, attribute] in self.yaml2attribute()])
        return rows
        
    def to_structure(self):
        struct = Sample.to_structure(self)
        for (yaml,attribute) in self.yaml2attribute():
            val = getattr(self,attribute)
            if val is not None:
                struct[yaml] = val
        return struct
  
    def _get_barcode_id(self):
        return _from_unicode(self._barcode_id)
    def _set_barcode_id(self,barcode_id):
        self._barcode_id = _to_unicode(barcode_id)
               
    def _get_barcode_name(self):
        return _from_unicode(self._barcode_name)
    def _set_barcode_name(self,barcode_name):
        self._barcode_name = _to_unicode(barcode_name)
               
    def _get_barcode_sequence(self):
        return _from_unicode(self._barcode_sequence)
    def _set_barcode_sequence(self,barcode_sequence):
        self._barcode_sequence = _to_unicode(barcode_sequence)
             
    def _get_barcode_type(self):
        return _from_unicode(self._barcode_type)
    def _set_barcode_type(self,barcode_type):
        self._barcode_type = _to_unicode(barcode_type)
        
    barcode_id = property(_get_barcode_id,_set_barcode_id)
    barcode_name = property(_get_barcode_name,_set_barcode_name)
    barcode_sequence = property(_get_barcode_sequence,_set_barcode_sequence)
    barcode_type = property(_get_barcode_type,_set_barcode_type)
    
