"""Task definitions for the Celery message queue (http://celeryproject.org/).
"""
import time

from celery.task import task

from bcbio.pipeline import sample, lane, toplevel, storage, shared, variation
from bcbio.variation import realign, genotype
from bcbio.google import sequencing_report

# TODO: Make things work without this temporary fix.
# See issue #191 on GitHub: https://github.com/SciLifeLab/bcbb/issues/191
import sys
sys.path.insert(0, "")

# Global configuration for tasks in the main celeryconfig module
import celeryconfig


@task(queue="google_docs")
def create_report_on_gdocs(*args):
    [fc_date, fc_name, run_info_yaml, dirs, config] = args
    return sequencing_report.create_report_on_gdocs(fc_date,fc_name,run_info_yaml,dirs,config)

@task(ignore_results=True, queue="toplevel")
def analyze_and_upload(*args):
    """Run full analysis and upload results to Galaxy instance.

    Workers need to run on the machine with Galaxy installed for upload,
    but the actual processing can be distributed to multiple nodes.
    """
    config_file = celeryconfig.BCBIO_CONFIG_FILE
    remote_info = args[0]
    toplevel.analyze_and_upload(remote_info, config_file)


@task(ignore_results=True, queue="toplevel")
def fetch_data(*args):
    """Transfer sequencing data from a remote machine. Could be e.g. a sequencer
    or a pre-processing machine.
    """
    config_file = celeryconfig.BCBIO_CONFIG_FILE
    remote_info = args[0]
    toplevel.fetch_data(remote_info, config_file)


@task(ignore_results=True, queue="toplevel")
def backup_data(*args):
    """Backup sequencing data from a remote machine. Could be e.g. a sequencer
    or a pre-processing machine.
    """
    config_file = celeryconfig.BCBIO_CONFIG_FILE
    remote_info = args[0]
    toplevel.backup_data(remote_info, config_file)


@task(ignore_results=True, queue="storage")
def long_term_storage(*args):
    config_file = celeryconfig.BCBIO_CONFIG_FILE
    remote_info = args[0]
    storage.long_term_storage(remote_info, config_file)


@task
def process_lane(*args):
    return lane.process_lane(*args)


@task
def remove_contaminants(*args):
    return lane.remove_contaminants(*args)


@task
def process_alignment(*args):
    return lane.process_alignment(*args)


@task
def mark_duplicates_sample(*args):
    return sample.mark_duplicates_sample(*args)


@task
def merge_sample(*args):
    return sample.merge_sample(*args)


@task
def recalibrate_sample(*args):
    return sample.recalibrate_sample(*args)


@task
def realign_sample(*args):
    return realign.realign_sample(*args)


@task
def process_sample(*args):
    return sample.process_sample(*args)


@task
def generate_bigwig(*args):
    return sample.generate_bigwig(*args)


@task
def combine_bam(*args):
    return shared.combine_bam(*args)


@task
def screen_sample_contaminants(*args):
    return sample.screen_sample_contaminants(*args)


@task
def variantcall_sample(*args):
    return genotype.variantcall_sample(*args)


@task
def combine_variant_files(*args):
    return genotype.combine_variant_files(*args)


@task
def detect_sv(*args):
    return variation.detect_sv(*args)


@task
def test(x):
    print x
    time.sleep(5)
    return x
