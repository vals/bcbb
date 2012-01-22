#!/usr/bin/env python
"""Setup analysis based on project identifiers

Usage:
  project_analysis_setup.py <YAML config file> <flow cell dir> <project dir>
                            [<YAML run information> --data_prefix=<data prefix>
                             --flowcell_alias=<flowcell alias>
                             --project_desc=<project_desc>
                             --lanes=<lanes> --barcode_ids=<barcode_ids>
                             --move_data --only_install_run_info
                             --only_install_fastq --dry_run --verbose]


Given a directory with demultiplexed flow cell data, and a project id
or a comma-separated list of lane numbers, project specific files will
be copied to a project directory. The optional <YAML run information>
file specifies details about the flowcell lanes, instead of retrieving
it from Galaxy. See automated_initial_analysis.py for details.

For a multiproject run_info file, only a subset of the lanes can be
used. The run_info file is therefore pruned, and the pruned file is
output to the project directory. The pruning is based on the options
<project_desc>, or <lanes>. Keyword ALL delivers all lanes.


Options:
  -d, --data_prefix=<data_prefix>               Install flowcells in <project_dir>/<data_prefix>
  -a, --flowcell_alias=<flowcell alias>         By default, samples are moved to a directory named
                                                <flowcell_id>. This option changes output directory to
                                                <flowcell_alias>.
  -y, --project_desc=<project_desc>             Project description in description field of run_info file, or ALL.
  -l, --lanes=<lanes>                           Comma-separated list of integers corresponding to lanes
  -b, --barcode_ids=<barcode_ids>               Comma-separated list of integers corresponding to barcode ids.
                                                Can only be used with one lane.
  -c, --customer_delivery                       Deliver data to customer. Delivers all demuxed fastq data to one directory.
  -i, --only_install_run_info                   Only install pruned run_info file.
  -f, --only_install_fastq                      Only install fastq files.
  -m, --move_data                               Move data instead of copying
  -n, --dry_run                                 Don't do anything samples, just list what will happen
  -v, --verbose                                 Print some more information
"""

import os
import sys
import re
from optparse import OptionParser

import yaml
import glob
import shutil
import logbook
from itertools import izip

from bcbio.log import create_log_handler
from bcbio.pipeline.run_info import get_run_info
from bcbio.pipeline.lane import get_flowcell_id
from bcbio.pipeline.fastq import get_single_fastq_files, get_barcoded_fastq_files, convert_barcode_id_to_name, get_fastq_files
from bcbio.pipeline.config_loader import load_config
from bcbio.pipeline.flowcell import Flowcell, Lane
from bcbio import utils

LOG_NAME = os.path.splitext(os.path.basename(__file__))[0]
log = logbook.Logger(LOG_NAME)

class DeliveryFlowcell(Flowcell):
    """A class for managing information about a flowcell.

    Adds functionality for keeping track of where data analysis 
    has taken place, fc alias for data delivery"""

    def __init__(self, fc_name, fc_date, data, fc_dir=None, fc_alias=None, fc_analysis_dir=None):
        #Flowcell.__init__(self,fc_name,fc_date,data, fc_dir)
        self.set_fc_dir(fc_dir)
        self.set_fc_date(fc_date)
        self.set_fc_name(fc_name)
        self.set_lanes(data)
        self.set_fc_alias(fc_alias)
        self.set_fc_analysis_dir(fc_analysis_dir)

    def get_fc_analysis_dir(self):
        return self.fc_analysis_dir
    def set_fc_analysis_dir(self,fc_analysis_dir):
        self.fc_analysis_dir = fc_analysis_dir

    def get_fc_alias(self):
        return self.fc_alias
    def set_fc_alias(self,fc_alias):
        self.fc_alias = fc_alias

    def get_fc_id(self):
        return "_".join([self.get_fc_date(), self.get_fc_name()])

    def prune_to_project(self,project, exclude_unmatched=True):
        """Return a new DeliveryFlowcell object just containing the lanes and samples belonging to a specific project"""
        lanes = []
        fc = None
        for lane in self.get_lanes():
            l = lane.prune_to_project(project, exclude_unmatched)
            if (l):
                lanes.append(l.to_structure())
        if (len(lanes)):
            fc = DeliveryFlowcell(self.get_fc_name(),self.get_fc_date(),lanes,self.get_fc_dir())
        return fc

    def __str__(self):
        s =  """DeliveryFlowcell: %s
    fc_dir:          %s
    fc_date:         %s
    fc_name:         %s
    fc_alias:        %s
    fc_analysis_dir: %s

    Number of lanes: %s
""" % (self.get_fc_id(), self.get_fc_dir(), self.get_fc_date(), self.get_fc_name(), self.get_fc_alias(), self.get_fc_analysis_dir(), len(self.lanes))
        return s

def main(config_file, fc_dir, project_dir, run_info_yaml=None, fc_alias=None, project_desc=None, lanes=None, barcodes=None):
    if project_desc is None and lanes is None:
        log.error("No project description or lanes provided: cannot deliver files without this information")
        sys.exit()

    config = load_config(config_file)
    ## Set log file in project output directory
    config.update(log_dir=os.path.join(project_dir, "log"))
    log_handler = create_log_handler(config)

    fc_dir = os.path.abspath(fc_dir)
    fc_name, fc_date, run_info = get_run_info(fc_dir, config, run_info_yaml)
    ## Alternativ with pipeline.flowcell.Flowcell
    fp = open(run_info_yaml)
    run_info_structure = yaml.load(fp)
    fc = DeliveryFlowcell(fc_name, fc_date, run_info_structure, fc_dir=fc_dir)
    with log_handler.applicationbound():
        #run_info = prune_run_info(run_info['details'], project_desc, lanes, barcodes)
        fc_src = fc.prune_to_project(project_desc, exclude_unmatched=True)
    if fc_src is None or len(fc.get_lanes()) == 0:
        if not project_desc is None:
            log.error("No lanes found with matching description %s: please check your flowcell run information" % project_desc)
            sys.exit()
        if not lanes  is None:
            log.error("No lanes found with numbers %s: please check your flowcell run information" % " ".join(lanes))
            sys.exit()
    # Set up a "target" flowcell that contains the delivery information (directory etc)
    fc_tgt = DeliveryFlowcell(fc_src.get_fc_name(), fc_src.get_fc_date(), fc_src.to_structure()['details'])
    fc_alias = "%s_%s" % (fc_date, fc_name) if not fc_alias else fc_alias
    fc_delivery_dir = os.path.abspath(os.path.join(project_dir, options.data_prefix, fc_alias))
    fc_analysis_dir = os.path.abspath(os.path.join(project_dir, options.data_prefix, "%s_%s" %(fc_date, fc_name)))
    if options.customer_delivery:
        fc_analysis_dir = os.path.abspath(os.path.join(project_dir, options.data_prefix))
    fc_tgt.set_fc_dir(fc_delivery_dir)
    fc_tgt.set_fc_alias(fc_alias)
    fc_tgt.set_fc_analysis_dir(fc_analysis_dir)

    with log_handler.applicationbound():
        _make_delivery_directory(fc_tgt)
        _save_run_info(fc_tgt, run_exit=options.only_run_info)
        run_main(fc_src, fc_tgt)

def run_main(fc_src, fc_tgt):
    for lane in fc_src.get_lanes():
        process_lane(lane, fc_src, fc_tgt)

def process_lane(lane, fc_src, fc_tgt):
    """Models bcbio process lane"""
    multiplex = lane.get_samples()
    log.info("Processing project: %s; lane %s; reference genome %s" %
             (lane.get_description(), lane.get_name(), lane.get_genome_build()))
    if multiplex:
        log.debug("Project %s is multiplexed as: %s" % (lane.get_description(), multiplex))
    fq = _get_barcoded_fastq_files(lane, multiplex, fc_src.get_fc_date(), fc_src.get_fc_name(), fc_src.get_fc_dir())
    ## Move data along with fastq files
    if not options.customer_delivery:
        fc_bc_dir = os.path.join(fc_tgt.get_fc_dir(), "%s_%s_%s_barcode" % (lane.get_name(), fc_tgt.get_fc_date(), fc_tgt.get_fc_name()))
    else:
        fc_bc_dir = fc_tgt.get_fc_dir()

    _make_dir(fc_bc_dir, "fastq.txt barcode directory")
    if not options.only_fastq:
        data, fastqc = _get_analysis_results(fc_src, lane)
        _deliver_data(data, fastqc, fc_tgt.get_fc_analysis_dir())


    #     [_deliver_fastq_file(fq_src, os.path.basename(fq_src), fc_bc_dir) for fq_src in fqpair]
    for fqpair in fq:
        for fastq_src in fqpair:
            fastq_tgt = fastq_src
            if options.customer_delivery:
                fastq_tgt = _convert_barcode_id_to_name(multiplex, fc_tgt.get_fc_name(), fastq_src)
            _deliver_fastq_file(fastq_src, os.path.basename(fastq_tgt), fc_bc_dir)


def _get_barcoded_fastq_files(lane, multiplex, fc_date, fc_name, fc_dir=None):
    fq = list()
    bc_dir = "%s_%s_%s_barcode" % (lane.get_name(), fc_date, fc_name)
    bc_dir = os.path.join(fc_dir, bc_dir)
    if multiplex is None:
        fq.append(get_fastq_files(bc_dir, {'lane':lane.get_name()}, fc_name))
    else:
        for bc in multiplex:
            if not os.path.exists(bc_dir):
                raise IOError("No barcode directory found: " + str(bc_dir))
            item = {'lane':lane.get_name()}
            fq.append(get_fastq_files(bc_dir, None, item, fc_name, bc_name=bc.get_barcode_id()))
    return fq

def _convert_barcode_id_to_name(multiplex, fc_name, fq):
    bcid2name = dict([(mp.get_barcode_id(), mp.get_barcode_name()) for mp in multiplex])
    bcid = re.search("_(\d+)_(\d+)_fastq.txt", fq)
    from_str = "%s_%s_fastq.txt" % (bcid.group(1), bcid.group(2))
    to_str   = "%s_%s.fastq" % (bcid2name[bcid.group(1)], bcid.group(2))
    return fq.replace(from_str, to_str)
 
def _deliver_fastq_file(fq_src, fq_tgt, outdir, fc_link_dir=None):
    _handle_data(fq_src, os.path.join(outdir, fq_tgt), f=shutil.move if options.move else shutil.copyfile)

def _make_delivery_directory(fc_tgt):
    """Make the output directory"""
    _make_dir(fc_tgt.get_fc_dir(), "flowcell delivery")
    _make_dir(fc_tgt.get_fc_analysis_dir(), "data delivery")
    if (os.path.basename(fc_tgt.get_fc_analysis_dir()) != fc_tgt.get_fc_alias() and not options.customer_delivery):
        _handle_data(fc_tgt.get_fc_analysis_dir(), os.path.join(os.path.dirname(fc_tgt.get_fc_analysis_dir()), fc_tgt.get_fc_alias()), os.symlink)

def _make_dir(dir, label):
    if not os.path.exists(dir):
        os.makedirs(dir)
        log.info("Creating %s directory %s" % (label, dir))
    else:
        log.warn("%s already exists: not creating new directory" % (dir))

def _handle_data(src, tgt, f=shutil.copyfile):
    if src is None:
        return
    if os.path.exists(tgt):
        log.warn("%s already exists: not doing anything!" %(tgt))
        return
    if options.dry_run:
        print "DRY_RUN: %s file %s to %s" % (f.__name__, src, tgt)
    else:
        log.info("%s file %s to %s" % (f.__name__, src, tgt))
        f(src, tgt)

def _deliver_data(data, fastqc, outdir):
    """Loop over data and fastqc and deliver files"""
    for src in data:
        tgt = os.path.join(outdir, os.path.basename(src))
        _handle_data(src, tgt, f=shutil.move if options.move else shutil.copyfile)

    for src in fastqc:
        tgt = os.path.join(outdir, "fastqc", os.path.basename(src))
        _handle_data(src, tgt, f=shutil.move if options.move else shutil.copytree)

def _get_analysis_results(fc, lane):
    """Get analysis results
    
    For now just glob the analysis directory for fastqc output and files with the give flowcell name
    """
    flowcell = "_".join([lane.get_name(), fc.get_fc_id()])
    glob_str = os.path.join(fc.get_fc_dir(), flowcell + "*.*")
    data = glob.glob(glob_str)
    glob_str = os.path.join(fc.get_fc_dir(), "fastqc", flowcell + "*")
    fastqc = glob.glob(glob_str)
    return data, fastqc

def _save_run_info(fc_tgt, run_exit=False):
    outfile = os.path.join(fc_tgt.get_fc_dir(), "project_run_info.yaml")
    if not options.dry_run:
        with open(outfile, "w") as out_handle:
            yaml.dump(fc_tgt.to_structure(), stream=out_handle)
    else:
        print "DRY_RUN:"
        yaml.dump(fc_tgt.to_structure(), stream=sys.stdout)
    if run_exit:
        sys.exit()

if __name__ == "__main__":
    usage = """
    project_analysis_setup.py <YAML config file> <flow cell dir> <project dir>
                            [<YAML run information> --data_prefix=<data prefix>
                             --flowcell_alias=<flowcell alias>
                             --project_desc=<project_desc>
                             --lanes=<lanes> --barcode_ids=<barcode_ids>
                             --move_data --only_install_run_info --only_install_fastq
                             --customer_delivery
                             --dry_run --verbose]

    For more extensive help type project_analysis_setup.py
"""

    parser = OptionParser(usage=usage)
    parser.add_option("-d", "--data_prefix", dest="data_prefix",
                      default="")
    parser.add_option("-a", "--flowcell_alias", dest="fc_alias")
    parser.add_option("-y", "--project_desc", dest="project_desc")
    parser.add_option("-l", "--lanes", dest="lanes")
    parser.add_option("-b", "--barcode_ids", dest="barcodes")

    parser.add_option("-i", "--only_install_fastq", dest="only_fastq", action="store_true",
                      default=False)
    parser.add_option("-f", "--only_install_run_info", dest="only_run_info", action="store_true",
                      default=False)
    parser.add_option("-m", "--move_data", dest="move", action="store_true",
                      default=False)
    parser.add_option("-c", "--customer_delivery", dest="customer_delivery", action="store_true",
                      default=False)

    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      default=False)
    parser.add_option("-n", "--dry_run", dest="dry_run", action="store_true",
                      default=False)
    (options, args) = parser.parse_args()
    if len(args) < 3:
        print __doc__
        sys.exit()
    kwargs = dict(
        fc_alias = options.fc_alias,
        project_desc = options.project_desc,
        lanes = options.lanes,
        barcodes = options.barcodes
        )
    main(*args, **kwargs)
