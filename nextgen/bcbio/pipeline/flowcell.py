"""A module for managing metadata and layout, as specified in run_info.yaml,
as well as read counts for the samples on a flowcell.
"""
import copy
import re
import os
import glob
import yaml
import string

from bcbio.solexa.flowcell import get_flowcell_info
from bcbio.google import (_to_unicode, _from_unicode)
from bcbio.utils import UnicodeReader
from bcbio.log import logger2

from bcbio.log import logger2 as log


def format_project_name(unformated_name):
    """Make the project name adhere to a stricter formatting convention.
    """
    regexp = r'^(.+?)_(\d{2})_(\d{2})(.*)$'
    m = re.match(regexp, unformated_name)
    if not m or len(m.groups()) < 3:
        return unformated_name

    name = m.group(1).strip()
    year = m.group(2).strip()
    month = m.group(3).strip()
    suffix = m.group(4).strip()

    # Replace any non-period delimiters
    delimiter = "_"
    p = re.compile('(%s)' % delimiter)
    name = p.sub('.', name)

    # Make sure that the initial letters are in capitals
    name = string.capwords(name, ".")

    # Format the name
    project_name = "%s_%s_%s%s" % (name, year, month, suffix)
    return project_name


def get_barcode_metrics(workdir):
    """Parse the [lane]_*_bc.metrics files in the *_barcode directories into
    a dictionary.

    If the samples have been demultiplexed by CASAVA, parses the
    Unaligned/Basecall_Stats_*/Demultiplex_Stats.htm file for barcode metrics.
    """
    bc_files = []
    if workdir is not None:
        bc_files.extend(glob.glob(os.path.join(workdir, "*_barcode", "*.bc_metrics")))
        bc_files.extend(glob.glob(os.path.join(workdir, "*.bc_metrics")))
    if not len(bc_files) > 0:
        return None

    bc_metrics = {}
    for bc_file in bc_files:

        m = re.match(r'^(\d+)\_', os.path.basename(bc_file))
        if not m or len(m.groups()) != 1:
            continue

        lane = str(m.group(1))
        bc_metrics[lane] = {}
        with open(bc_file) as bcfh:
            csvr = UnicodeReader(bcfh, dialect='excel-tab')
            for row in csvr:
                bc_metrics[lane][str(row[0])] = int(row[1])

    return bc_metrics


def get_flowcell(fc_dir, run_info_yaml, config={}):
    # Just get the name of the flowcell directory minus the path
    fc_name, fc_date = get_flowcell_info(os.path.basename(os.path.normpath(fc_dir)))
    with open(run_info_yaml, "r") as fh:
        run_info = yaml.load(fh)

    return Flowcell(fc_name, fc_date, run_info, fc_dir)


def get_project_name(description):
    """Parse out the project name from the lane description"""
    m = re.match(r'(?:.*\s+)?(\S+)', (description or ""), re.I)
    if m and len(m.groups()) > 0:
        return format_project_name(m.group(1).strip())

    return description

def get_sample_name(barcode_name):
    """Extract the sample name by stripping the barcode index part of the sample description""" 
    name, index = split_sample_name(barcode_name)
    return name
                     
def split_sample_name(sample_name):
    """Split a sample name into parts consisting of 
        - project_name [PNNN]
        - sample_number [NNN]
        - reception_qc [F]
        - prep_version [B]
        - index_id [indexN]
    """
    
    splits = sample_name.split("_")
    prep = ""
    try:
        if len(splits) < 2:
            raise ValueError()
        if splits[0][0] != 'P':
            raise ValueError()
        if type(int(splits[0][1:])) != int:
            raise ValueError()
        while splits[1][-1] in "FB":
            prep = "%c%s" % (splits[1][-1],prep)
            splits[1] = splits[1][0:-1]
        if type(int(splits[1])) != int:
            raise ValueError()
    except:
        logger2.warn("Sample name '%s' does not follow the expected format PXXX_XXX[FB]" % sample_name)
    if len(prep) > 0:
        splits[1] = "%s%s" % (splits[1],prep)
        
    name = []
    index = []
    for s in splits:
        if len(index) == 0 and s.find('index') < 0:
            name.append(s)
        else:
            index.append(s)
    return "_".join(name), "_".join(index)
        
class Flowcell:
    """A class for managing information about a flowcell"""

    def __init__(self, fc_name, fc_date, data, fc_dir=None):
        # Extract the run_items if we are passed a dictionary
        try:
            log.debug("Try making flowcell with this data:")
            log.debug(data)
            d = data.get('details', [])
            data = d
        except AttributeError:
            pass

        self.set_fc_dir(fc_dir)
        self.set_fc_date(fc_date)
        self.set_fc_name(fc_name)
        self.set_lanes(data)
        # Attempts to set the read counts on creation
        self.set_read_counts()

    def get_fc_date(self):
        return self.fc_date

    def set_fc_date(self, fc_date):
        self.fc_date = fc_date

    def get_fc_dir(self):
        return self.fc_dir

    def set_fc_dir(self, fc_dir):
        self.fc_dir = fc_dir

    def get_fc_name(self):
        return self.fc_name

    def set_fc_name(self, fc_name):
        self.fc_name = fc_name

    def get_lane_by_name(self, name):
        for lane in self.get_lanes():
            if (lane.get_name() == name):
                return lane
        return None

    def add_lane(self, lane):
        self.lanes.append(lane)

    def get_lanes(self):
        return self.lanes

    def set_lanes(self, run_items):
        self.lanes = []
        for lane in run_items:
            self.add_lane(Lane(lane))

    def get_project_names(self):
        pnames = set()
        for lane in self.get_lanes():
            for pname in lane.get_project_names():
                if pname is None:
                    continue

                pnames.add(pname)

        return list(pnames)

    def get_samples(self):
        samples = []
        for lane in self.get_lanes():
            samples.extend(lane.get_samples())
        return samples

    def prune_to_project(self, project, exclude_unmatched=False):
        """Return a new Flowcell object just containing the lanes and samples belonging to a specific project"""
        lanes = []
        fc = None
        for lane in self.get_lanes():
            l = lane.prune_to_project(project, exclude_unmatched)
            if (l):
                lanes.append(l.to_structure())
        if (len(lanes)):
            fc = Flowcell(self.get_fc_name(), self.get_fc_date(), lanes)
        return fc

    def set_read_counts(self, read_counts=None):
        """Sets the read counts of the barcoded samples in the lanes of this
        flowcell. Read counts can be supplied in a dictionary with lane number
        as key to a dictionary with barcode indexes (or 'unmatched' or 'trim')
        as key and counts as values. If no read counts are supplied, attempts
        to parse the read counts from the flowcell directory, assuming that the
        read counts can be found using a glob like
        [lane]_*_barcode/[lane]_*.bc_metrics
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
        cols = []  # "fc_date","fc_name"]
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
    """A class for managing information about a lane.
    """

    def __init__(self, data):
        self.set_data(data)
        self.set_description(data.get("description", None))
        self.set_name(data.get("lane", None))
        self.set_samples(data.get("multiplex", []))
        self.set_files([])

    def get_analysis(self):
        return self.data.get("analysis", None)

    def set_data(self, data):
        self.data = copy.deepcopy(data)

    def get_description(self):
        return _from_unicode(self.description)

    def set_description(self, description):
        self.description = _to_unicode(description)

    def get_genome_build(self):
        return self.data.get("genome_build", None)

    def get_name(self):
        return _from_unicode(self.name)

    def set_name(self, name):
        self.name = _to_unicode(str(name))

    def get_project_names(self):
        pnames = {}
        for sample in self.get_samples():
            pnames[sample.get_project()] = 1
        return pnames.keys()

    def set_read_counts(self, read_counts):
        for barcode_id in read_counts.keys():
            sample = self.get_sample_by_barcode(barcode_id)
            if (sample is None):
                sample = BarcodedSample({"name": "Unexpected barcode", "barcode_id": barcode_id})
                self.add_sample(sample)
            sample.set_read_count(read_counts[barcode_id])

    def get_sample_by_barcode(self, barcode_id):
        for sample in self.get_samples():
            if (str(sample.get_barcode_id()) == str(barcode_id)):
                return sample
        return None

    def add_sample(self, sample):
        self.multiplex.append(sample)

    def get_samples(self):
        return self.multiplex

    def set_samples(self, multiplex):
        self.multiplex = []
        for barcode in multiplex:
            self.add_sample(BarcodedSample(barcode, self))

    def set_files(self, files):
        self.files = files

    def get_files(self):
        return self.files

    def get_samples_by_project(self, project):
        samples = []
        for sample in self.get_samples():
            if (sample.get_project() == project):
                samples.append(sample)
        return samples

    def prune_to_project(self, project, exclude_unmatched=False):
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
        cols = ["lane", "description"]
        cols.extend(BarcodedSample.columns())
        return cols

    def to_rows(self):
        rows = []
        for sample in self.get_samples():
            s = sample.to_rows()
            s.insert(0, self.get_name())
            s.insert(1, self.get_description())
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
        s = "Lane: %s\nbarcode ids: %s" % (self.get_name(), self.get_barcode_ids())
        return s


class Sample:
    """A class for managing information about a sample.
    """

    def __init__(self, data, lane=Lane({}), comment=None):
        self.set_analysis(data.get("analysis", lane.get_analysis()))
        self.set_genome_build(data.get("genome_build", lane.get_genome_build()))
        self.set_name(data.get("name", None))
        self.set_full_name(data.get("full_name", data.get("name", None)))
        self.set_project(data.get("sample_prj", \
            data.get("description", lane.get_description())))
        self.set_read_count(data.get("read_count", None))
        self.set_description(data.get("description", None))
        self.set_lane(lane.get_name())
        self.set_comment(comment)

    def add_sample(self, other, delim=', '):
        if (self.get_analysis() != other.get_analysis()):
            self.set_analysis("%s%s%s" \
                % (self.get_analysis(), delim, other.get_analysis()))
        if (self.get_genome_build() != other.get_genome_build()):
            self.set_genome_build("%s%s%s" \
                % (self.get_genome_build(), delim, other.get_genome_build()))
        if (self.get_name() != other.get_name()):
            self.set_name("%s%s%s" \
                % (self.get_name(), delim, other.get_name()))
        if (self.get_project() != other.get_project()):
            self.set_project("%s%s%s" \
                % (self.get_project(), delim, other.get_project()))
        if (self.get_lane() != other.get_lane()):
            self.set_lane("%s%s%s" \
                % (self.get_lane(), delim, other.get_lane()))
        if (other.get_read_count() is not None):
            self.set_read_count((self.get_read_count() or 0) \
                + other.get_read_count())

    def get_analysis(self):
        return _from_unicode(self.analysis)

    def set_analysis(self, analysis):
        self.analysis = _to_unicode(analysis)

    def get_genome_build(self):
        return _from_unicode(self.genome_build)

    def set_genome_build(self, genome_build):
        self.genome_build = _to_unicode(genome_build)

    def get_name(self):
        return _from_unicode(self.name)

    def set_name(self, name):
        self.name = get_sample_name(_to_unicode(name))

    def get_full_name(self):
        return _from_unicode(self.full_name)

    def set_full_name(self, name):
        self.full_name = _to_unicode(name)

    def get_comment(self):
        return self.comment

    def set_comment(self, comment):
        self.comment = comment

    def get_lane(self):
        return self.lane

    def set_lane(self, lane):
        self.lane = lane

    def get_project(self):
        return _from_unicode(self.project)

    def set_project(self, project):
        self.project = get_project_name(_to_unicode(project))

    def get_description(self):
        return _from_unicode(self.description)

    def set_description(self, description):
        self.description = _to_unicode(description)

    def get_read_count(self):
        if self.read_count:
            try:
                rc = int(self.read_count)
                return rc
            except ValueError:
                pass
        return None

    def get_rounded_read_count(self, unit=1000000, decimals=2):
        return round((self.get_read_count() or 0) / float(unit), int(decimals))

    def set_read_count(self, read_count):
        self.read_count = read_count

    @staticmethod
    def columns():
        cols = ["project_name", "sample_name", "read_count", "rounded_read_count"]
        return cols

    def to_rows(self):
        rows = [self.get_project(), self.get_name(), self.get_read_count(), self.get_rounded_read_count()]
        return rows

    def to_structure(self):
        struct = {}
        if (self.get_analysis()):
            struct["analysis"] = self.get_analysis()
        if (self.get_genome_build()):
            struct["genome_build"] = self.get_genome_build()
        if (self.get_name()):
            struct["name"] = self.get_name()
        if (self.get_full_name()):
            struct["full_name"] = self.get_full_name()
        if (self.get_project()):
            struct["sample_prj"] = self.get_project()
        if (self.get_description()):
            struct["description"] = self.get_description()
        if (self.get_read_count() is not None):
            struct["read_count"] = self.get_read_count()

        return struct


class BarcodedSample(Sample):
    """A subclass of Sample for managing information about a barcoded sample.
    """
    def __init__(self, data, lane=Lane({})):
        Sample.__init__(self, data, lane)
        self.set_barcode_id(data.get("barcode_id", None))
        self.set_barcode_name(data.get("name", None))
        self.set_barcode_sequence(data.get("sequence", None))
        self.set_barcode_type(data.get("barcode_type", None))
        self.set_barcode_full_name(data.get("full_name", None))

    def get_barcode_id(self):
        return _from_unicode(self.barcode_id)

    def set_barcode_id(self, barcode_id):
        self.barcode_id = _to_unicode(barcode_id)

    def get_barcode_name(self):
        return _from_unicode(self.barcode_name)

    def set_barcode_name(self, barcode_name):
        self.barcode_name = _to_unicode(barcode_name)

    def get_barcode_full_name(self):
        return _from_unicode(self.barcode_full_name)

    def set_barcode_full_name(self, barcode_full_name):
        self.barcode_full_name = _to_unicode(barcode_full_name)

    def get_barcode_sequence(self):
        return _from_unicode(self.barcode_sequence)

    def set_barcode_sequence(self, barcode_sequence):
        self.barcode_sequence = _to_unicode(barcode_sequence)

    def get_barcode_type(self):
        return _from_unicode(self.barcode_type)

    def set_barcode_type(self, barcode_type):
        self.barcode_type = _to_unicode(barcode_type)

    @staticmethod
    def columns():
        cols = Sample.columns()
        cols.extend(["bcbb_barcode_id", "barcode_name", "barcode_sequence", "barcode_type"])

        return cols

    def to_rows(self):
        rows = Sample.to_rows(self)
        rows.extend([self.get_barcode_id(), \
                     self.get_barcode_name(), \
                     self.get_barcode_sequence(), \
                     self.get_barcode_type()])
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
