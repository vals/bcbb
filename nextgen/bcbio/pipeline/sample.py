"""High level entry point for processing a sample.

Samples may include multiple lanes, or barcoded subsections of lanes,
processed together.
"""
import os
import subprocess


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
    logger.info("Recalibrating %s with GATK" % str(data["name"]))
    if data["config"]["algorithm"]["recalibrate"]:
        recal_bam = recalibrate_quality(data["work_bam"], data["fastq1"],
                                        data["fastq2"], data["sam_ref"],
                                        data["dirs"], data["config"])
        save_diskspace(data["work_bam"], "Recalibrated to %s" % recal_bam,
                       data["config"])
        data["work_bam"] = recal_bam
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


# ## General processing

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


#def _filter_out_genomes(data):
#    """ Filters out genomes found in run_info.yaml
#    """
#    #genome_build, sam_ref = ref_genome_info(data["info"], config, data["dirs"])
#    sam_ref = data["sam_ref"]
#    
#    log.info("Removing genome from sample %s" % str(data["name"]))
#    try:
#        # XXX should be able to reach run_info.yaml info instead of data here
#        if data["filter_out_genomes"]:
#            for genome in data["filter_out_genomes"].split(","):
#                (out_file, ext) = os.path.splitext(os.path.basename(fastq1))
#                out_file = out_file+"-stripped-"+genome+ext
#                cl = ["bowtie", "--solexa1.3-quals", "--un", out_file, sam_ref, "-1", data["fastq1"], "-2", data["fastq2"], "/dev/null"]
#                log.info("Running %s" % cl)
#                subprocess.check_call(cl)
#    except KeyError:
#        log.error("Not removing genomes, directive filter_out_genomes undefined in run_info.yaml")
#        pass
