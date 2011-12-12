"""High level entry point for processing a sample.

Samples may include multiple lanes, or barcoded subsections of lanes,
processed together.
"""
import os
import subprocess

from bcbio.utils import file_transaction
from bcbio.pipeline.lane import _update_config_w_custom
from bcbio.pipeline import log
from bcbio.pipeline.alignment import get_genome_ref
from bcbio.pipeline.merge import (combine_fastq_files, merge_bam_files)
from bcbio.pipeline.qcsummary import generate_align_summary, screen_for_contamination
from bcbio.pipeline.variation import (recalibrate_quality, run_genotyper,
                                      variation_effects)
from bcbio.rnaseq.cufflinks import assemble_transcripts


# ## General processing

def process_sample(sample_name, fastq_files, info, bam_files, dirs,
                   config, config_file):
    """Finalize processing for a sample, potentially multiplexed.
    """
    config = _update_config_w_custom(config, info)

    genome_build = info.get("genome_build", None)
    (_, sam_ref) = get_genome_ref(genome_build, config["algorithm"]["aligner"],
                                  dirs["galaxy"])
    fastq1, fastq2 = combine_fastq_files(fastq_files, dirs["work"])
    log.info("Combining and preparing wig file %s" % str(sample_name))
    sort_bam = merge_bam_files(bam_files, dirs["work"], config)
    (gatk_bam, vrn_file, effects_file) = ("", "", "")
    if config["algorithm"]["screen_contaminants"]:
        log.info("Screening for contaminants on sample %s with genome %s" % \
        (str(sample_name), str(genome_build)))
        screen_for_contamination(fastq1, fastq2, config, genome_build)

    # _filter_out_genomes(data)

    if config["algorithm"]["recalibrate"]:
        log.info("Recalibrating %s with GATK" % str(sample_name))
        gatk_bam = recalibrate_quality(sort_bam, fastq1, fastq2, sam_ref,
                                       dirs, config)
        if config["algorithm"]["snpcall"]:
            log.info("SNP genotyping %s with GATK" % str(sample_name))
            vrn_file = run_genotyper(gatk_bam, sam_ref, config)
            log.info("Calculating variation effects for %s" % str(sample_name))
            effects_file = variation_effects(vrn_file, genome_build, config)

    if config["algorithm"].get("transcript_assemble", False):
        tx_file = assemble_transcripts(sort_bam, sam_ref, config)

    if sam_ref is not None:
        log.info("Generating summary files: %s" % str(sample_name))
        generate_align_summary(sort_bam, fastq2 is not None, sam_ref,
                               sample_name, config, dirs)
    bam_to_wig(sort_bam, config, config_file)

    return [sample_name, fastq_files, info, sort_bam, gatk_bam, vrn_file,
            effects_file]


def bam_to_wig(bam_file, config, config_file):
    """Provide a BigWig coverage file of the sorted alignments.
    """
    wig_file = "%s.bigwig" % os.path.splitext(bam_file)[0]
    if not (os.path.exists(wig_file) and os.path.getsize(wig_file) > 0):
        cl = [config["analysis"]["towig_script"], bam_file, config_file]
        with file_transaction(wig_file):
            cl = [os.path.expandvars(command) for command in cl]
            subprocess.check_call(cl)
    return wig_file


def _filter_out_genomes(data):
    """ Filters out genomes found in run_info.yaml
    """
    print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    print data
    #genome_build, sam_ref = ref_genome_info(data["info"], config, data["dirs"])
    sam_ref = data["sam_ref"]

    log.info("Removing genome from sample %s" % str(data["name"]))
    try:
        # not data ! should reach run_info.yaml somehow from here
        if data["filter_out_genomes"]:
            for genome in data["filter_out_genomes"].split(","):
                (out_file, ext) = os.path.splitext(os.path.basename(fastq1))
                out_file = out_file+"-stripped-"+genome+ext
                cl = ["bowtie", "--solexa1.3-quals", "--un", out_file, sam_ref, "-1", data["fastq1"], "-2", data["fastq2"], "/dev/null"]
                log.info("Running %s" % cl)
                subprocess.check_call(cl)
    except KeyError:
        log.error("Not removing genomes, directive filter_out_genomes undefined in run_info.yaml")
        pass
