#!/usr/bin/env python
"""Perform an automated analysis on a sequencing run using Galaxy information.

Given a directory of solexa output, this retrieves details about the sequencing
run from the Galaxy description, and uses this to perform an initial alignment
and analysis.

Usage:
    automated_initial_analysis.py <YAML config file> <flow cell dir>
                                  [<YAML run information>]

The optional <YAML run information> file specifies details about the
flowcell lanes, instead of retrieving it from Galaxy. An example
configuration file is located in 'config/run_info.yaml'

Workflow:
    - Retrieve details on a run.
    - Align fastq files to reference genome.
    - Perform secondary analyses like SNP calling.
    - Generate summary report.
"""
import os
import sys
from optparse import OptionParser

import datetime
import yaml

from bcbio.solexa.flowcell import get_fastq_dir
from bcbio import utils
from bcbio.log import logger, setup_logging, version
from bcbio.distributed.messaging import parallel_runner
from bcbio.pipeline.run_info import get_run_info
from bcbio.pipeline.demultiplex import add_multiplex_across_lanes
from bcbio.pipeline.merge import organize_samples
from bcbio.pipeline.qcsummary import write_metrics, write_project_summary
from bcbio.variation.realign import parallel_realign_sample
from bcbio.variation.genotype import parallel_variantcall
from bcbio.pipeline.config_loader import load_config
from bcbio.google.sequencing_report import create_report_on_gdocs


def main(config_file, fc_dir, run_info_yaml=None):
    config = load_config(config_file)
    work_dir = os.getcwd()
    if config.get("log_dir", None) is None:
        config["log_dir"] = os.path.join(work_dir, "log")
    setup_logging(config)
    run_main(config, config_file, fc_dir, work_dir, run_info_yaml)


def run_main(config, config_file, fc_dir, work_dir, run_info_yaml):
    
    _record_sw_versions(config, os.path.join(work_dir,"bcbb_software_versions.txt"))
    
    align_dir = os.path.join(work_dir, "alignments")
    run_module = "bcbio.distributed"
    fc_name, fc_date, run_info = get_run_info(fc_dir, config, run_info_yaml)
    fastq_dir, galaxy_dir, config_dir = _get_full_paths(get_fastq_dir(fc_dir),
                                                        config, config_file)
    config_file = os.path.join(config_dir, os.path.basename(config_file))
    dirs = {"fastq": fastq_dir, "galaxy": galaxy_dir, "align": align_dir,
            "work": work_dir, "flowcell": fc_dir, "config": config_dir}
    run_parallel = parallel_runner(run_module, dirs, config, config_file)

    run_items = add_multiplex_across_lanes(run_info["details"], dirs["fastq"], fc_name)
    lanes = ((info, fc_name, fc_date, dirs, config) for info in run_items)
    lane_items = run_parallel("process_lane", lanes)
    
    # upload the sequencing report to Google Docs
    gdocs_indicator = os.path.join(work_dir,"gdocs_report_complete.txt")
    if not os.path.exists(gdocs_indicator) and create_report_on_gdocs(fc_date, fc_name, run_info_yaml, dirs, config):
        utils.touch_file(gdocs_indicator)

    # Remove spiked in controls, contaminants etc.
    lane_items = run_parallel("remove_contaminants",lane_items)
    
    align_items = run_parallel("process_alignment", lane_items)
    # process samples, potentially multiplexed across multiple lanes
    samples = organize_samples(align_items, dirs, config_file)
    samples = run_parallel("merge_sample", samples)
    run_parallel("screen_sample_contaminants", samples)
    samples = run_parallel("recalibrate_sample", samples)
    samples = parallel_realign_sample(samples, run_parallel)
    samples = parallel_variantcall(samples, run_parallel)
    samples = run_parallel("detect_sv", samples)
    samples = run_parallel("process_sample", samples)
    samples = run_parallel("generate_bigwig", samples, {"programs": ["ucsc_bigwig"]})
    write_project_summary(samples)
    write_metrics(run_info, fc_name, fc_date, dirs)


# ## Utility functions

def _record_sw_versions(config, sw_version_file):
    """Get the versions of software used in the pipeline and output to
       log and text file in working directory
    """
    sw_versions = version.get_versions(config)
    sw_versions['bcbb'] = version._get_git_commit()

    logger.info("bcbb pipeline is running with software versions: %s" % sw_versions)
    
    with open(sw_version_file,'w') as fh:
        fh.write("%s\n" % datetime.datetime.now().isoformat())
        for sw, ver in sw_versions.items():
            fh.write("%s\t%s\n" % (sw,ver))

def _get_full_paths(fastq_dir, config, config_file):
    """Retrieve full paths for directories in the case of relative locations.
    """
    fastq_dir = utils.add_full_path(fastq_dir)
    config_dir = utils.add_full_path(os.path.dirname(config_file))
    galaxy_config_file = utils.add_full_path(config["galaxy_config"], config_dir)
    return fastq_dir, os.path.dirname(galaxy_config_file), config_dir


def _get_run_info(fc_name, fc_date, config, run_info_yaml):
    """Retrieve run information from a passed YAML file or the Galaxy API.
    """
    if run_info_yaml and os.path.exists(run_info_yaml):
        logger.info("Found YAML samplesheet, using %s instead of Galaxy API" % run_info_yaml)
        with open(run_info_yaml) as in_handle:
            run_details = yaml.load(in_handle)
        return dict(details=run_details, run_id="")
    else:
        logger.info("Fetching run details from Galaxy instance")
        galaxy_api = GalaxyApiAccess(config['galaxy_url'], config['galaxy_api_key'])
        return galaxy_api.run_details(fc_name, fc_date)

if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict()
    main(*args, **kwargs)
