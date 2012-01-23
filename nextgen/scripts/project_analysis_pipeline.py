#!/usr/bin/env python
"""Perform bcbb pipeline on (demultiplexed) project data.

This script currently performs the same analyses as
automated_initial_analysis.py except multiplexing. May start with
demultiplexed data, or continue where an automated_initial_analysis.py
left off.

Usage:
    project_analysis_pipeline.py <YAML config file> <flow cell dir>
                                 [<pruned YAML run information>]

The optional <pruned YAML run information> file specifies details
about the flowcell lanes relevant for this project. See
automated_initial_analysis.py and project_analysis_setup.py for more
information.
"""


import os
import sys
from optparse import OptionParser

import yaml

import subprocess
import glob
import collections

from bcbio.log import create_log_handler, setup_logging
from bcbio.distributed.messaging import parallel_runner
from bcbio.pipeline.run_info import get_run_info
from bcbio.pipeline.qcsummary import write_metrics, write_project_summary
from bcbio.variation.realign import parallel_realign_sample
from bcbio.variation.genotype import parallel_variantcall
from bcbio import utils
from bcbio.pipeline.lane import _update_config_w_custom
from bcbio.pipeline.merge import organize_samples
from bcbio.pipeline.config_loader import load_config

def main(config_file, fc_dir, run_info_yaml=None):
    config = load_config(config_file)
    work_dir = os.getcwd()
    if config.get("log_dir", None) is None:
        config["log_dir"] = os.path.join(work_dir, "log")
    setup_logging(config)
    run_main(config, config_file, fc_dir, work_dir, run_info_yaml)

def run_main(config, config_file, fc_dir, work_dir, run_info_yaml):
    # Working directory has to be identical to where (demultiplexed) fastq files are located
    align_dir = os.path.join(work_dir, "alignments")
    run_module = "bcbio.distributed"
    fc_name, fc_date, run_info = get_run_info(fc_dir, config, run_info_yaml)
    fc_dir = os.path.normpath(fc_dir)
    galaxy_dir, config_dir = _get_full_paths(config, config_file)

    config_file = os.path.join(config_dir, os.path.basename(config_file))

    dirs = dict(galaxy= galaxy_dir, 
                align = align_dir, work = work_dir,
                config = config_dir, flowcell = fc_dir,
                fc_dir = fc_dir,     fastq = fc_dir
                )
    run_parallel = parallel_runner(run_module, dirs, config, config_file)

    # Since demultiplexing is already done, just extract run_items
    run_items = run_info['details']
    lanes = ((info, fc_name, fc_date, dirs, config) for info in run_items)
    lane_items = []
    #lane_items = (fastq1, fastq2, barcode_info, barcode_prefix, description: name, dirs, post_process)
    for args in lanes:
        lane_items.extend(_make_lane_items(*args))
    
    align_items = run_parallel("process_alignment", lane_items)
    # Process samples
    samples = organize_samples(align_items, dirs, config_file)
    samples = run_parallel("merge_sample", samples)
    samples = run_parallel("recalibrate_sample", samples)
    samples = parallel_realign_sample(samples, run_parallel)
    samples = parallel_variantcall(samples, run_parallel)
    samples = run_parallel("process_sample", samples)
    samples = run_parallel("generate_bigwig", samples, {"programs": ["ucsc_bigwig"]})
    write_project_summary(samples)
    write_metrics(run_info, fc_name, fc_date, dirs)


def _make_lane_items(lane_items, fc_name, fc_date, dirs, config):
    """Return lane items equivalent to output from run_parallel("process_lane", lanes)"""
    lane_name = "%s_%s_%s" % (lane_items[0]['lane'], fc_date, fc_name)
    bc_files = _get_barcode_files(lane_items, lane_name, dirs, config)
    out = []
    for item in lane_items:
        config = _update_config_w_custom(config, item)
        if bc_files.has_key(item["barcode_id"]):
            fastq1, fastq2 = bc_files[item["barcode_id"]]
            cur_lane_name = lane_name
            cur_lane_desc = item["description"]
            if item.get("name", "") and config["algorithm"].get("include_short_name", True):
                cur_lane_desc = "%s : %s" % (item["name"], cur_lane_desc)
            if item["barcode_id"] is not None:
                cur_lane_name += "_%s" % (item["barcode_id"])
            out.append((fastq1, fastq2, item, cur_lane_name, cur_lane_desc,
                        dirs, config))
    return out

def _get_barcode_files(multiplex, base_name, dirs, config):
    bc_dir = os.path.join(dirs["fc_dir"], "%s_barcode" % base_name)
    out_files = []
    for info in multiplex:
        fq_fname = lambda x: os.path.abspath(os.path.join(bc_dir, "%s_%s_%s_fastq.txt" %
                             (base_name, info["barcode_id"], x)))
        bc_file1 = fq_fname("1")
        bc_file2 = fq_fname("2") if os.path.exists(fq_fname("2")) else None
        out_files.append((info["barcode_id"], bc_file1, bc_file2))
    out = {}
    for b, f1, f2 in out_files:
        if os.path.exists(f1):
            out[b] = (f1, f2)
    return out

def _get_full_paths(config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    config_dir = utils.add_full_path(os.path.dirname(config_file))
    galaxy_config_file = utils.add_full_path(config["galaxy_config"], config_dir)
    return os.path.dirname(galaxy_config_file), config_dir

if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict()
    main(*args, **kwargs)

