"""A module for managing metadata and layout, as specified in run_info.yaml, as well as read counts for the samples on a flowcell"""

import copy
import re

from bcbio.pipeline.run_info import get_run_info
from bcbio.google import ( _to_unicode, _from_unicode )

def format_project_name(unformated_name):
    """Make the project name adhere to the formatting convention"""
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

def get_flowcell(fc_dir, run_info_yaml, config={}):
    fc_name, fc_date, run_info = get_run_info(fc_dir,config,run_info_yaml)
    return Flowcell(fc_name,fc_date,run_info)

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
    
    def __init__(self,fc_name,fc_date,data):
        self.set_fc_date(fc_date)
        self.set_fc_name(fc_name)
        self.set_lanes(data.get("details",[]))
        
    def collapse_on_sample(self):
        """Merge all data based on sample name, creating summarized read counts and concatenated lanes as necessary"""
        samples = {}
        for lane in self.get_lanes():
            for sample in lane.get_samples():
                if sample.get_name() in samples:
                    samples[sample.get_name()].add_sample(sample)
                else:
                    samples[sample.get_name()] = BarcodedSample(sample.to_structure())
        
    def get_fc_date(self):
        return self.fc_date
    def set_fc_date(self,fc_date):
        self.fc_date = fc_date
        
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
            fc = Flowcell(self.get_fc_name(),self.get_fc_date(),{"details": lanes})
        return fc
    
    def set_read_counts(self,read_counts):
        for name in read_counts.keys():
            lane = self.get_lane_by_name(name)
            if not lane:
                lane = Lane({"lane": name, "description": "Unexpected lane"})
                self.add_lane(lane)
            lane.set_read_counts(read_counts[name])
    
    def to_rows(self):
        rows = []
        for lane in self.get_lanes():
            for row in lane.to_rows():
                l = [self.get_fc_date(),self.get_fc_name()]
                l.extend(row)
                rows.append(l)
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
        self.set_samples(data.get("multiplex",[]))
        
    def get_analysis(self):
        return self.data.get("analysis",None)
    
    def set_data(self,data):
        self.data = copy.deepcopy(data)
    
    def get_description(self):
        return _from_unicode(self.description)
    def set_description(self,description):
        self.description = _to_unicode(description)
        
    def get_genome_build(self):
        return self.data.get("genome_build",None)
    
    def get_name(self):
        return _from_unicode(self.name)
    def set_name(self,name):
        self.name = _to_unicode(name)
        
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
    
    def add_sample(self,sample):
        self.multiplex.append(sample)
    def get_samples(self):
        return self.multiplex
    def set_samples(self,multiplex):
        self.multiplex = []
        for barcode in multiplex:
            self.add_sample(BarcodedSample(barcode,self))
        
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
        
    def to_rows(self):
        rows = []
        for sample in self.get_samples():
            s = sample.to_rows()
            s.insert(0,self.get_name())
            rows.append(s)
        return rows
            
    def to_structure(self):
        struct = {}
        if (self.get_description()):
            struct["description"] = self.get_description()
        if (self.get_name()):
            struct["lane"] = self.get_name()
        if (self.get_samples()):
            struct["multiplex"] = []
            for sample in self.get_samples():
                struct["multiplex"].append(sample.to_structure())
        return struct
        
class Sample:
    """A class for managing information about a sample"""
     
    def __init__(self,data,lane=Lane({})):
        self.set_analysis(data.get("analysis",lane.get_analysis()))
        self.set_genome_build(data.get("genome_build",lane.get_genome_build()))
        self.set_name(data.get("name",None))
        self.set_project(data.get("description",lane.get_description()))
        self.set_read_count(data.get("read_count",None))
        
    def add_sample(self,other):
        if (self.get_analysis() != other.get_analysis()):
            self.set_analysis("%s, %s" % (self.get_analysis(),other.get_analysis()))
        if (self.get_genome_build() != other.get_genome_build()):
            self.set_genome_build("%s, %s" % (self.get_genome_build(),other.get_genome_build()))
        if (self.get_name() != other.get_name()):
            self.set_name("%s, %s" % (self.get_name(),other.get_name()))
        if (self.get_project() != other.get_project()):
            self.set_project("%s, %s" % (self.get_project(),other.get_project()))
        if (other.get_read_count() is not None):
            self.set_read_count((self.get_read_count() or 0) + other.get_read_count())
            
    def get_analysis(self):
        return _from_unicode(self.analysis)
    def set_analysis(self,analysis):
        self.analysis = _to_unicode(analysis)
        
    def get_genome_build(self):
        return _from_unicode(self.genome_build)
    def set_genome_build(self,genome_build):
        self.genome_build = _to_unicode(genome_build)
        
    def get_name(self):
        return _from_unicode(self.name)
    def set_name(self,name):
        self.name = get_sample_name(_to_unicode(name))
        
    def get_project(self):
        return _from_unicode(self.project)
    def set_project(self,project):
        self.project = get_project_name(_to_unicode(project))
        
    def get_read_count(self):
        return self.read_count
    def get_rounded_read_count(self,unit=1000000,decimals=2):
        return round(int((self.get_read_count() or 0))/float(unit),int(decimals))
    def set_read_count(self,read_count):
        self.read_count = read_count
    
    def to_rows(self):
        rows = [self.get_project(),self.get_name(),self.get_read_count(),self.get_rounded_read_count()]
        return rows
    
    def to_structure(self):
        struct = {}
        if (self.get_analysis()):
            struct["analysis"] = self.get_analysis()
        if (self.get_genome_build()):
            struct["genome_build"] = self.get_genome_build()
        if (self.get_name()):
            struct["name"] = self.get_name()
        if (self.get_project()):
            struct["description"] = self.get_project()
        if (self.get_read_count() is not None):
            struct["read_count"] = self.get_read_count()
        return struct
    
class BarcodedSample(Sample):
    """A subclass of Sample for managing information about a barcoded sample"""
    
    def __init__(self,data,lane=Lane({})):
        Sample.__init__(self,data,lane)
        self.set_barcode_id(data.get("barcode_id",None))
        self.set_barcode_name(data.get("name",None))
        self.set_barcode_sequence(data.get("sequence",None))
        self.set_barcode_type(data.get("barcode_type",None))

    def get_barcode_id(self):
        return _from_unicode(self.barcode_id)
    def set_barcode_id(self,barcode_id):
        self.barcode_id = _to_unicode(barcode_id)
               
    def get_barcode_name(self):
        return _from_unicode(self.barcode_name)
    def set_barcode_name(self,barcode_name):
        self.barcode_name = _to_unicode(barcode_name)
               
    def get_barcode_sequence(self):
        return _from_unicode(self.barcode_sequence)
    def set_barcode_sequence(self,barcode_sequence):
        self.barcode_sequence = _to_unicode(barcode_sequence)
             
    def get_barcode_type(self):
        return _from_unicode(self.barcode_type)
    def set_barcode_type(self,barcode_type):
        self.barcode_type = _to_unicode(barcode_type)
    
    def to_rows(self):
        rows = Sample.to_rows(self)
        rows.extend([self.get_barcode_id(), self.get_barcode_name(), self.get_barcode_sequence(), self.get_barcode_type()])
        return rows
    
    def to_structure(self):
        struct = Sample.to_structure(self)
        if (self.get_barcode_id()):
            struct["barcode_id"] = self.get_barcode_id()
        if (self.get_barcode_sequence()):
            struct["sequence"] = self.get_barcode_sequence()
        if (self.get_barcode_type()):
            struct["barcode_type"] = self.get_barcode_type()
        return struct
  
    