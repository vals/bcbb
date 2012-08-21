"""High level entry point for processing a sample.

Samples may include multiple lanes, or barcoded subsections of lanes,
processed together.
"""
import os
import subprocess


from bcbio import broad
from bcbio.broad.picardrun import picard_mark_duplicates
from bcbio.utils import file_exists, save_diskspace
from bcbio.distributed.transaction import file_transaction
from bcbio.pipeline.lane import _update_config_w_custom
from bcbio.log import logger
from bcbio.pipeline.merge import (combine_fastq_files, merge_bam_files)
from bcbio.pipeline.qcsummary import generate_align_summary, screen_for_contamination
from bcbio.pipeline.variation import (recalibrate_quality, finalize_genotyper,
                                      variation_effects)
from bcbio.rnaseq.cufflinks import assemble_transcripts
from bcbio.pipeline.shared import ref_genome_info


def merge_sample(data):
    """Merge fastq and BAM files for multiple samples.
    """
    logger.info("Combining fastq and BAM files %s" % str(data["name"]))
    config = _update_config_w_custom(data["config"], data["info"])
    genome_build, sam_ref = ref_genome_info(data["info"], config, data["dirs"])
    fastq1, fastq2 = combine_fastq_files(data["fastq_files"], data["dirs"]["work"],
                                         config)
    sort_bam = merge_bam_files(data["bam_files"], data["dirs"]["work"], config)

    return [[{"name": data["name"],
              "genome_build": genome_build, "sam_ref": sam_ref,
              "work_bam": sort_bam, "fastq1": fastq1, "fastq2": fastq2,
              "dirs": data["dirs"], "config": config,
              "config_file": data["config_file"]}]]


def recalibrate_sample(data):
    """Recalibrate quality values from aligned sample BAM file.
    """
    logger.info("Recalibrating {} with GATK".format(str(data["name"])))
    if data["config"]["algorithm"]["recalibrate"]:
        recal_bam = recalibrate_quality(data["work_bam"], data["fastq1"],
                                        data["fastq2"], data["sam_ref"],
                                        data["dirs"], data["config"])
        save_diskspace(data["work_bam"], \
                       "Recalibrated to {}".format(recal_bam), data["config"])
        data["work_bam"] = recal_bam

    return [[data]]


def mark_duplicates_sample(data):
    """Mark duplicate molecules in sample BAM file.
    """
    if not data["config"]["algorithm"].get("mark_duplicates", False):
        return [[data]]

    logger.info("Marking duplicates in {} with Picard".format(str(data["name"])))
    picard = broad.runner_from_config(data["config"])
    dup_bam, _ = picard_mark_duplicates(picard, data["work_bam"])
    reason = "Marked duplicates of {0} in {1}, so {0} is no longer needed" \
             "".format(data["work_bam"], dup_bam)
    save_diskspace(data["work_bam"], reason, data["config"])
    data["work_bam"] = dup_bam

    return [[data]]


def screen_sample_contaminants(data):
    """Screen the sample fastq files for contaminants
    """
    if data["config"]["algorithm"]["screen_contaminants"]:
        logger.info("Screening for contaminants on sample %s with genome %s" \
        % (str(data["name"]), str(data["genome_build"])))
        screen_for_contamination(data["fastq1"],
                                 data["fastq2"],
                                 data["config"])


# General processing

def process_sample(data):
    """Finalize processing for a sample, potentially multiplexed.
    """
    if data["config"]["algorithm"]["snpcall"]:
        logger.info("Finalizing variant calls: %s" % str(data["name"]))
        data["vrn_file"] = finalize_genotyper(data["vrn_file"], data["work_bam"],
                                              data["sam_ref"], data["config"])
        logger.info("Calculating variation effects for %s" % str(data["name"]))
        ann_vrn_file, effects_file = variation_effects(data["vrn_file"],
                                                       data["sam_ref"],
                                                       data["genome_build"],
                                                       data["config"])
        if ann_vrn_file:
            data["vrn_file"] = ann_vrn_file
            data["effects_file"] = effects_file

    if data["config"]["algorithm"].get("transcript_assemble", False):
        data["tx_file"] = assemble_transcripts(data["work_bam"],
                                               data["sam_ref"],
                                               data["config"])
    if data["sam_ref"] is not None:
        logger.info("Generating summary files: %s" % str(data["name"]))
        generate_align_summary(data["work_bam"], data["fastq2"] is not None,
                               data["sam_ref"], data["name"],
                               data["config"], data["dirs"])

    return [[data]]


def generate_bigwig(data):
    """Provide a BigWig coverage file of the sorted alignments.
    """
    logger.info("Preparing BigWig file %s" % str(data["name"]))
    bam_file = data["work_bam"]
    wig_file = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if not file_exists(wig_file):
        with file_transaction(wig_file) as tx_file:
            cl = [data["config"]["analysis"]["towig_script"], bam_file,
                  data["config_file"], "--outfile=%s" % tx_file]
            subprocess.check_call(cl)

    data["bigwig_file"] = wig_file
    return [[data]]
