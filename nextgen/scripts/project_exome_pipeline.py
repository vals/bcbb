#!/usr/bin/env python
"""Provide snp and indel calling on project data.

This script currently performs the same analyses as
automated_initial_analysis.py except multiplexing. May start with
demultiplexed data, or continue where an automated_initial_analysis.py
left off.

Usage:
    project_exome_pipeline.py <YAML config file> <flow cell dir>
                              [<pruned YAML run information>]

The optional <pruned YAML run information> file specifies details
about the flowcell lanes relevant for this project. See
automated_initial_analysis.py and project_analysis_setup.py for more
information.
"""


import os
import sys
import copy
from optparse import OptionParser

import yaml

import subprocess
import glob
import collections

from bcbio.log import create_log_handler
from bcbio.pipeline import log
from bcbio.pipeline.lane import make_lane_items, get_flowcell_id
from bcbio.pipeline.run_info import get_run_info
from bcbio.solexa.flowcell import get_flowcell_info
from bcbio.galaxy.api import GalaxyApiAccess
from bcbio import utils
from bcbio.pipeline.demultiplex import add_multiplex_across_lanes
from bcbio.broad import BroadRunner
from bcbio.ngsalign import bwa
from bcbio.pipeline import lane
from bcbio.pipeline import sample
from bcbio.pipeline.merge import organize_samples

def main(config_file, fc_dir, run_info_yaml=None):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    log_handler = create_log_handler(config, log.name)
    with log_handler.applicationbound():
        run_main(config, config_file, fc_dir, run_info_yaml)

def run_main(config, config_file, fc_dir, run_info_yaml):
    # Working directory has to be identical to where (demultiplexed) fastq files are located
    fc_dir = os.path.normpath(fc_dir)
    work_dir = os.getcwd()
    align_dir = os.path.join(work_dir, "alignments")

    #(_, fastq_dir_label) = os.path.split(work_dir)
    #fastq_dir = os.path.join(project_dir, fastq_dir_label)

    #fc_name, fc_date = get_flowcell_info(fc_dir)
    fc_name, fc_date, run_info = get_run_info(fc_dir, config, run_info_yaml)
    #run_info = _get_run_info(fc_name, fc_date, config, run_info_yaml)
    #fastq_dir, galaxy_dir, config_dir = _get_full_paths(fastq_dir, config, config_file)
    galaxy_dir, config_dir = _get_full_paths(config, config_file)

    config_file = os.path.join(config_dir, os.path.basename(config_file))
#    fastq = fastq_dir,
    dirs = dict(galaxy= galaxy_dir, 
                align = align_dir, work = work_dir,
                config = config_dir, flowcell = fc_dir,
                fc_dir = fc_dir
                )
    # Since demultiplexing is already done, just extract run_items
    run_items = run_info['details']
    lane_items = []
    for info in run_items:
        print info
        config = _update_config_w_custom(config, info)
        lane_items.extend(make_lane_items(info, fc_date, fc_name, dirs, config))

    _run_parallel("process_alignment", lane_items, dirs, config)
    
    # Process samples
    sample_files, sample_fastq, sample_info = \
                  organize_samples(dirs, fc_name, fc_date, run_items)
    samples = ((n, sample_fastq[n], sample_info[n], bam_files, dirs, config, config_file)
               for n, bam_files in sample_files)
    _run_parallel("process_sample", samples, dirs, config)

# def _get_run_info(fc_name, fc_date, config, run_info_yaml):
#     """Retrieve run information from a passed YAML file or the Galaxy API.
#     """
#     if run_info_yaml and os.path.exists(run_info_yaml):
#         log.info("Found YAML samplesheet, using %s instead of Galaxy API" % run_info_yaml)
#         with open(run_info_yaml) as in_handle:
#             run_details = yaml.load(in_handle)
#         return dict(details=run_details, run_id="")
#     else:
#         log.info("Fetching run details from Galaxy instance")
#         #galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
#         #return galaxy_api.run_details(fc_name, fc_date)


def _run_parallel(fn_name, items, dirs, config):
    """Process a supplied function: single, multi-processor or distributed.
    """
    parallel = config["algorithm"]["num_cores"]
    if str(parallel).lower() == "messaging":
        runner = messaging.runner(dirs, config)
        return runner(fn_name, items)
    else:
        out = []
        fn = globals()[fn_name]
        with utils.cpmap(int(parallel)) as cpmap:
            for data in cpmap(fn, items):
                if data:
                    out.extend(data)
        return out

@utils.map_wrap
def process_alignment(*args):
    return lane.process_alignment(*args)

@utils.map_wrap
def process_sample(*args):
    return sample.process_sample(*args)

def _get_lane_info_from_fastq(fastq):
    fh = os.path.basename(fastq)
    lane = fh.split("_")[0]
    pu = "_".join(fh.split("_")[0:3])
    return lane, pu
    
def _get_full_paths(config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    #fastq_dir = utils.add_full_path(fastq_dir)
    config_dir = utils.add_full_path(os.path.dirname(config_file))
    galaxy_config_file = utils.add_full_path(config["galaxy_config"], config_dir)
    return os.path.dirname(galaxy_config_file), config_dir

def _update_config_w_custom(config, lane_info):
    """Update the configuration for this lane if a custom analysis is specified.
    """
    config = copy.deepcopy(config)
    analysis_type = lane_info.get("analysis", "")
    custom = config["custom_algorithms"].get(analysis_type, None)
    if custom:
        for key, val in custom.iteritems():
            config["algorithm"][key] = val
    # apply any algorithm details specified with the lane
    for key, val in lane_info.get("algorithm", {}).iteritems():
        config["algorithm"][key] = val
    return config


if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict()
    main(*args, **kwargs)

