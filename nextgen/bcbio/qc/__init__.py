import os
import sys
import re
import yaml
import xml.parsers.expat
import hashlib
import time

import glob
import json
import fnmatch
import numpy as np

from bcbio.log import logger2 as log
from bcbio.broad.metrics import *
from bcbio.pipeline.qcsummary import FastQCParser

class MetricsParser():
    """Basic class for parsing metrics"""
    def __init__(self):
        pass

    def parse_bc_metrics(self, in_handle):
        data = {}
        while 1:
            line = in_handle.readline()
            if not line:
                break
            vals = line.rstrip("\t\n\r").split("\t")
            data[vals[0]] = int(vals[1])
        return data
    
    def parse_filter_metrics(self, in_handle):
        data = {}
        data["reads"] = int(in_handle.readline().rstrip("\n").split(" ")[-1])
        data["reads_aligned"] = int(in_handle.readline().split(" ")[-2])
        data["reads_fail_align"] = int(in_handle.readline().split(" ")[-2])
        return data

    def parse_fastq_screen_metrics(self, in_handle):
        column_names = ["Library", "Unmapped", "Mapped_One_Library", "Mapped_Multiple_Libraries"]
        in_handle.readline()
        data = {}
        while 1:
            line = in_handle.readline()
            if not line:
                break
            vals = line.rstrip("\t\n").split("\t")
            data[vals[0]] = {}
            data[vals[0]]["Unmapped"] = float(vals[1])
            data[vals[0]]["Mapped_One_Library"] = float(vals[2])
            data[vals[0]]["Mapped_Multiple_Libraries"] = float(vals[3])
        return data


class ExtendedPicardMetricsParser(PicardMetricsParser):
    """Extend basic functionality and parse all picard metrics"""

    def __init__(self):
        PicardMetricsParser.__init__(self)

    def _get_command(self, in_handle):
        analysis = None
        while 1:
            line = in_handle.readline()
            if line.startswith("# net.sf.picard.analysis") or line.startswith("# net.sf.picard.sam"):
                break
        return line.rstrip("\n")

    def _read_off_header(self, in_handle):
        while 1:
            line = in_handle.readline()
            if line.startswith("## METRICS"):
                break
        return in_handle.readline().rstrip("\n").split("\t")

    def _read_vals_of_interest(self, want, header, info):
        want_indexes = [header.index(w) for w in header]
        vals = dict()
        for i in want_indexes:
            vals[header[i]] = info[i]
        return vals

    def _parse_align_metrics(self, in_handle):
        command = self._get_command(in_handle)
        header = self._read_off_header(in_handle)
        d = dict([[x, []] for x in header])
        res = dict(command=command, FIRST_OF_PAIR = d, SECOND_OF_PAIR = d, PAIR = d)
        while 1:
            info = in_handle.readline().rstrip("\n").split("\t")
            category = info[0]
            if len(info) <= 1:
                break
            vals = self._read_vals_of_interest(header, header, info)
            res[category] = vals
        return res

    def _parse_dup_metrics(self, in_handle):
        command = self._get_command(in_handle)
        header = self._read_off_header(in_handle)
        info = in_handle.readline().rstrip("\n").split("\t")
        vals = self._read_vals_of_interest(header, header, info)
        histvals = self._read_histogram(in_handle)
        return dict(command=command, metrics = vals, hist = histvals)

    def _parse_insert_metrics(self, in_handle):
        command = self._get_command(in_handle)
        header = self._read_off_header(in_handle)
        info = in_handle.readline().rstrip("\n").split("\t")
        vals = self._read_vals_of_interest(header, header, info)
        histvals = self._read_histogram(in_handle)
        return dict(command=command, metrics = vals, hist = histvals)

    def _parse_hybrid_metrics(self, in_handle):
        command = self._get_command(in_handle)
        header = self._read_off_header(in_handle)
        info = in_handle.readline().rstrip("\n").split("\t")
        vals = self._read_vals_of_interest(header, header, info)
        return dict(command=command, metrics = vals)
    
    def _read_histogram(self, in_handle):
        labels = self._read_to_histogram(in_handle)
        if labels is None:
            return None
        vals = dict([[x, []] for x in labels])
        while 1:
            line = in_handle.readline()
            info = line.rstrip("\n").split("\t")
            if len(info) < len(labels):
                break
            for i in range(0, len(labels)):
                vals[labels[i]].append(info[i])
        return vals

    def _read_to_histogram(self, in_handle):
        while 1:
            line = in_handle.readline()
            if line.startswith("## HISTOGRAM"):
                break
            if not line:
                return None
        return in_handle.readline().rstrip("\n").split("\t")
          
class RunInfoParser():
    """RunInfo parser"""
    def __init__(self):
        self._data = {}
        self._element = None

    def parse(self, fp):
        if os.path.exists(fp):
            self._parse_RunInfo(fp)
        return self._data

    def _start_element(self, name, attrs):
        self._element=name
        if name == "Run":
            self._data["Id"] = attrs["Id"]
            self._data["Number"] = attrs["Number"]
        elif name == "FlowcellLayout":
            self._data["FlowcellLayout"] = attrs
        elif name == "Read":
            self._data["Reads"].append(attrs)
            
    def _end_element(self, name):
        self._element=None

    def _char_data(self, data):
        want_elements = ["Flowcell", "Instrument", "Date"]
        if self._element in want_elements:
            self._data[self._element] = data
        if self._element == "Reads":
            self._data["Reads"] = []

    def _parse_RunInfo(self, fp):
        p = xml.parsers.expat.ParserCreate()
        p.StartElementHandler = self._start_element
        p.EndElementHandler = self._end_element
        p.CharacterDataHandler = self._char_data
        p.ParseFile(fp)


class IlluminaXMLParser():
    """Illumina xml data parser. Parses xml files in flowcell directory."""
    def __init__(self):
        self._data = {}
        self._element = None
        self._tmp = None
        self._header = None

    def _chart_start_element(self, name, attrs):
        self._element = name
        if name == "FlowCellData":
            self._header = attrs
        if name == "Layout":
            n_tiles_per_lane = int(attrs['RowsPerLane']) * int(attrs['ColsPerLane'])
            nrow = int(attrs['NumLanes']) * n_tiles_per_lane
            if self._tmp is None:
                self._tmp = self._header
                self._tmp.update(attrs)
                for i in range(1,int(attrs['NumLanes'])+1):
                    for j in range(1,n_tiles_per_lane+1):
                        key = "%s_%s" % (i,j)
                        self._tmp[key] = {}
                     
        if name == "TL":
            self._tmp[attrs["Key"]][self._index] = {}
            for k in attrs.keys():
                if k == "Key":
                    continue
                if attrs[k] == "NaN":
                    v = None
                else:
                    v = float(attrs[k])
                                        
                self._tmp[attrs["Key"]][self._index] = v

    def _chart_end_element(self, name):
        self._element = None
    def _chart_char_data(self, data):
        pass

    def _parse_charts(self, files):
        for f in files:
            self._index = os.path.basename(f).rstrip(".xml").lstrip("Chart_")
            p = xml.parsers.expat.ParserCreate()
            p.StartElementHandler = self._chart_start_element
            p.EndElementHandler = self._chart_end_element
            p.CharacterDataHandler = self._chart_char_data
            fp = open(f)
            p.ParseFile(fp)
            fp.close()

    def _summary_start_element(self, name, attrs):
        self._element = name
        if name == "Summary":
            self._tmp[self._index] = attrs
        if name == "Lane":
            self._tmp[self._index][attrs['key']] = attrs
    def _summary_end_element(self, name):
        self._element = None
    def _summary_char_data(self, data):
        pass

    def _parse_summary(self, files):
        for f in files:
            self._index = os.path.basename(f).rstrip(".xml")
            p = xml.parsers.expat.ParserCreate()
            p.StartElementHandler = self._summary_start_element
            p.EndElementHandler = self._summary_end_element
            p.CharacterDataHandler = self._summary_char_data
            fp = open(f)
            p.ParseFile(fp)
            fp.close()

    def _clusters_start_element(self, name, attrs):
        self._element = name
        if name == "Data":
            self._tmp[self._index] = attrs
        if name == "Lane":
            self._tmp[self._index][attrs['key']] = attrs
    def _clusters_end_element(self, name):
        self._element = None
    def _clusters_char_data(self, data):
        pass

    def _parse_clusters(self, files):
        for f in files:
            self._index = os.path.basename(f).rstrip(".xml")
            p = xml.parsers.expat.ParserCreate()
            p.StartElementHandler = self._clusters_start_element
            p.EndElementHandler = self._clusters_end_element
            p.CharacterDataHandler = self._clusters_char_data
            fp = open(f)
            p.ParseFile(fp)
            fp.close()

    ## Caution: no assert statements for file existence
    def parse(self, files, fullRTA=False):
        """Full parsing includes all RTA files"""
        if fullRTA:
            error_files = filter(lambda x: os.path.dirname(x).endswith("ErrorRate"), files)
            self._parse_charts(error_files)
            self._data["ErrorRate"] = self._tmp
            self._tmp = None
            FWHM_files = filter(lambda x: os.path.dirname(x).endswith("FWHM"), files)
            self._parse_charts(FWHM_files)
            self._data["FWHM"] = self._tmp
            self._tmp = None
            intensity_files = filter(lambda x: os.path.dirname(x).endswith("Intensity"), files)
            self._parse_charts(intensity_files)
            self._data["Intensity"] = self._tmp
            self._tmp = None
            numgt30_files = filter(lambda x: os.path.dirname(x).endswith("NumGT30"), files)
            self._parse_charts(numgt30_files)
            self._data["NumGT30"] = self._tmp
            self._tmp = None
            chart_files = filter(lambda x: os.path.basename(x).endswith("_Chart.xml"), files)
            self._parse_charts(chart_files)
            self._data["Charts"] = self._tmp
            self._tmp = None

        ## Parse Summary and clusters
        self._tmp = {}
        summary_files = filter(lambda x: os.path.dirname(x).endswith("Summary"), files)
        self._parse_summary(summary_files)
        self._data["Summary"] = self._tmp
        self._tmp = {}
        cluster_files = filter(lambda x: os.path.basename(x).startswith("NumClusters By"), files)
        self._parse_clusters(cluster_files)
        self._data["NumClusters"] = self._tmp

        return self._data

class ExtendedFastQCParser(FastQCParser):
    def __init__(self, base_dir):
        FastQCParser.__init__(self, base_dir)

    def get_fastqc_summary(self):
        metric_labels = ["Per base sequence quality", "Basic Statistics", "Per sequence quality scores",
                         "Per base sequence content", "Per base GC content", "Per sequence GC content",
                         "Per base N content", "Sequence Length Distribution", "Sequence Duplication Levels",
                         "Overrepresented sequences", "Kmer Content"]
        metrics = {x : self._to_dict(self._fastqc_data_section(x)) for x in metric_labels}
        return metrics

    def _to_dict(self, section):
        if len(section) == 0:
            return {}
        header = [x.strip("#") for x in section[0].rstrip("\t").split("\t")]
        d = []
        for l in section[1:]:
            d.append(l.split("\t"))
        data = np.array(d)
        df = {header[i]:data[:,i].tolist() for i in range(0,len(header))}
        return df

##############################
## QCMetrics objects
##############################
class QCMetrics(dict):
    """Generic QCMetrics class"""
    _entity_version = 0.1
    _metrics = []

    def __init__(self):
        self["_id"] = self.get_db_id()
        self["entity_type"] = self.entity_type()
        self["entity_version"] = self.entity_version()
        self["name"] = self.name()
        self["creation_time"] = None
        self["modification_time"] = None
        
    def entity_version(self):
        return str(self._entity_version)

    def entity_type(self):
        return type(self).__name__
    
    # FIXME: should raise error: QCMetrics must be subclassed
    def name(self):
        return "%s" % (self.get("entity_type", None))

    def get_id(self):
        return self.name()
    
    def get_db_id(self):
        return hashlib.md5(self.get_id()).hexdigest()
    
class LaneQCMetrics(QCMetrics):
    """Lane level class for holding qc data"""
    _entity_version = 0.1
    
    def __init__(self, flowcell, date, lane):
        self["lane"] = lane
        self["flowcell"] = flowcell
        self["date"] = date
        self["bc_metrics"] = {}
        self["filter_metrics"] = {}
        QCMetrics.__init__(self)

    def name(self):
        return "%s_%s_%s" % (self.get("lane"), self.get("date"), self.get("flowcell"))

class SampleQCMetrics(QCMetrics):
    """Sample-level class for holding qc metrics data"""
    _entity_version = 0.1
    _metrics = ["picard_metrics","fastqc","fastq_scr"]

    def __init__(self, flowcell, date, lane, barcode_name, barcode_id, sample_prj, sequence=None, barcode_type=None, genomes_filter_out=None, customer_prj=None, customer_sample_name=None):
        self["flowcell"] = flowcell
        self["date"] = date
        self["lane"] = lane
        self["barcode_name"] = barcode_name
        self["barcode_id"] = barcode_id
        self["sample_prj"] = sample_prj
        self["customer_prj"] = customer_prj
        self["customer_sample_name"] = customer_sample_name
        self["sequence"] = sequence
        self["barcode_type"] = barcode_type
        self["genomes_filter_out"] = genomes_filter_out
        self["bc_count"] = None
        self["metrics"] = {}
        QCMetrics.__init__(self)
        self.picard_files = []
        for m in self._metrics:
            self["metrics"][m] = {}

    def name(self):
        return "%s_%s_%s_%s" % (self.get("lane", None), self.get("date", None), self.get("flowcell", None), self.get("barcode_id", None))

    def get_name(self, nophix=False):
        if nophix:
            s = "%s_%s_%s_nophix_%s" % (self["lane"], self["date"], self["flowcell"], self["barcode_id"])
        else:
            s = "%s_%s_%s_%s"  % (self["lane"], self["date"], self["flowcell"], self["barcode_id"])
        return s

class FlowcellQCMetrics(QCMetrics):
    """Flowcell level class for holding qc data."""
    _metrics = ["RunInfo", "run_info_yaml"]
    _entity_version = 0.1

    def __init__(self, fc_date, fc_name, run_info_yaml, flowcell_dir, archive_dir, runinfo="RunInfo.xml", parse=True, fullRTA=False):
        self.flowcell_dir = flowcell_dir
        self.archive_dir = archive_dir
        self.run_id = "%s_%s" % (fc_date, fc_name)
        self.db=None
        self.sample = dict()
        self["lane"] = dict()
        self["metrics"] = dict()
        for m in self._metrics:
            self["metrics"][m] = None
        ## initialize runinfo in case no RunInfo.xml
        self["metrics"]["RunInfo"] = {"Id" : self.run_id, "Flowcell":fc_name, "Date": fc_date, "Instrument": "NA"}
        self._parseRunInfo(runinfo)
        QCMetrics.__init__(self)
        if parse:
            self.parse_run_info_yaml(run_info_yaml)
            self.read_picard_metrics()
            self.read_fastqc_metrics()
            self.parse_filter_metrics()
            self.parse_fastq_screen()
            self.parse_bc_metrics()
            self.parse_illumina_metrics(fullRTA)

    def _parseRunInfo(self, fn="RunInfo.xml"):
        log.info("_parseRunInfo")
        try:
            fp = open(os.path.join(self.archive_dir, fn))
            parser = RunInfoParser()
            data = parser.parse(fp)
            fp.close()
            self["metrics"]["RunInfo"] = data
            self.run_id = data.get("Id")
        except:
            log.warn("No such file %s" % os.path.join(self.archive_dir, fn))


    def parse_run_info_yaml(self, run_info_yaml):
        log.info("parse_run_info_yaml")
        fp = open(run_info_yaml)
        runinfo = yaml.load(fp)
        fp.close()
        for info in runinfo:
            if not self["lane"].has_key(info["lane"]):
                lane = LaneQCMetrics(self.get_full_flowcell(), self.get_date(), info["lane"])
                self["lane"][info["lane"]] = lane
                ## Add sample for unmatched data
                sample = SampleQCMetrics(self.get_full_flowcell(), self.get_date(), info["lane"], "unmatched", "unmatched", "NA", "NA", "NA", "NA")
                bc_index = "%s_%s" % (info["lane"], "unmatched")
                self.sample[bc_index] = sample
            ## Lane could be empty
            try:
                for mp in info["multiplex"]:
                    sample = SampleQCMetrics(self.get_full_flowcell(), self.get_date(), info["lane"], mp["name"], mp["barcode_id"], mp.get("sample_prj", None), mp["sequence"], mp["barcode_type"], mp.get("genomes_filter_out", None))
                    bc_index = "%s_%s" % (info["lane"], mp["barcode_id"])
                    self.sample[bc_index] = sample
            except:
                log.warn("No multiplexing information for lane %s" % info['lane'])
        self["metrics"]["run_info_yaml"] = runinfo

    def get_full_flowcell(self):
        vals = self["metrics"]["RunInfo"]["Id"].split("_")
        return vals[-1]
    def get_flowcell(self):
        return self.get("metrics").get("RunInfo").get("Flowcell")
    def get_date(self):
        return self.get("metrics").get("RunInfo").get("Date")
    def get_run_name(self):
        return "%s_%s" % (self.get_date(), self.get_full_flowcell())
    def name(self):
        return str(self.run_id)

    def to_json(self):
        samples = [self.sample[s] for s in self.sample]
        return json.dumps({'metrics':self["metrics"], 'samples':samples})

    def read_picard_metrics(self):
        log.info("read_picard_metrics")
        picard_parser = ExtendedPicardMetricsParser()
        files = self._get_metrics(self.flowcell_dir)
        # Group files to samples
        for fn in files:
            fn_tgt = os.path.basename(fn).replace("_nophix_", "_")
            re_str = r'(\d+)_(\d+)_([A-Z0-9a-z]+XX)_?([a-zA-Z0-9]+)?-'
            m = re.search(re_str, fn_tgt)
            (lane, date, flowcell, bc) = m.groups()
            bc_index = "%s_%s" % (lane, bc)
            if self.sample.has_key(bc_index):
                self.sample[bc_index].picard_files.append(fn)
                #print >> sys.stderr, "reading metrics %s for sample %s" % (fn, bc_index)
            else:
                log.warn("no sample %s for metrics %s" % (bc_index, fn))

        for s in self.sample:
            metrics = picard_parser.extract_metrics(self.sample[s].picard_files)
            self.sample[s]["metrics"]["picard_metrics"] = metrics

    def _get_metrics(self, indir, re_str='.*.(align|hs|insert|dup)_metrics'):
        matches = []
        for root, dirnames, filenames in os.walk(indir):
            for fn in filenames:
                if re.match(re_str, fn):
                    matches.append(os.path.join(root, fn))
        return matches


    def parse_filter_metrics(self, re_str="*filter[_.]metrics"):
        log.info("parse_filter_metrics")
        for l in self["lane"].keys():
            glob_str = os.path.join(self.flowcell_dir, "nophix", "%s_%s_%s%s" % (l, self.get_date(), self.get_full_flowcell(), re_str))
            f = glob.glob(glob_str)
            self["lane"][l]["filter_metrics"] = {"reads":None, "reads_aligned":None, "reads_fail_align":None}
            try:
                fp = open(f[0])
                parser = MetricsParser()
                data = parser.parse_filter_metrics(fp)
                fp.close()
                self["lane"][l]["filter_metrics"] = data
            except:
                log.warn("No filter nophix metrics for lane %s" % l)

            
    def parse_fastq_screen(self):
        log.info("parse_fastq_screen")
        for s in self.sample:
            g = os.path.join(self.flowcell_dir, "fastq_screen", "%s_%s_[nophix_]*%s_[12]_fastq_screen.txt" % (self.sample[s]["lane"], self.get_run_name(), self.sample[s]["barcode_id"]))
            f = glob.glob(g)
            if len(f) < 1:
                continue
            if not os.path.exists(f[0]):
                continue
            parser = MetricsParser()
            fp = open(f[0])
            data = parser.parse_fastq_screen_metrics(fp)
            fp.close()
            self.sample[s]["metrics"]["fastq_scr"] = data

    def parse_bc_metrics(self, re_str="*.bc_metrics"):
        log.info("parse_bc_metrics")
        for l in self["lane"].keys():
            glob_str = os.path.join(self.flowcell_dir, "%s_%s_*barcode" % (l, self.get_run_name()), re_str)
            f = glob.glob(glob_str)
            try:
                parser = MetricsParser()
                fp = open(f[0])
                data = parser.parse_bc_metrics(fp)
                fp.close()
                for key in data.keys():
                    s = "%s_%s" % (l, key)
                    self.sample[s]["bc_count"] = data[key]
                self["lane"][l]["bc_metrics"] = data
            except:
                log.warn("No bc_metrics info for lane %s: glob %s" % (s, glob_str))
                    
                
    def read_fastqc_metrics(self):
        log.info("read_fastqc_metrics")
        for s in self.sample:
            if s.endswith("unmatched"):
                continue
            self.sample[s]["metrics"]["fastqc"] = {'stats':None}
            glob_str = os.path.join(self.flowcell_dir, "fastqc", "%s_%s*%s-*" % (self.sample[s]["lane"], self.get_run_name(), self.sample[s]["barcode_id"]))
            d = glob.glob(glob_str)
            if len(d) == 0:
                log.warn("No fastqc metrics info for sample %s: glob %s" % (s, glob_str) )
                continue
            fastqc_dir=d[0]
            fqparser = ExtendedFastQCParser(fastqc_dir)
            stats = fqparser.get_fastqc_summary()
            self.sample[s]["metrics"]["fastqc"] = {'stats':stats}

    def parse_illumina_metrics(self, fullRTA):
        log.info("parse_illumina_metrics")
        fn = []
        for root, dirs, files in os.walk(self.archive_dir):
            for file in files:
                if file.endswith(".xml"):
                    fn.append(os.path.join(root, file))
        parser = IlluminaXMLParser()
        metrics = parser.parse(fn, fullRTA)
        self["metrics"]["illumina"] = metrics
