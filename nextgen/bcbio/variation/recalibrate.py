"""Perform quality score recalibration with the GATK toolkit.

Corrects read quality scores post-alignment to provide improved estimates of
error rates based on alignments to the reference genome.

http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
"""
import os
import shutil

from bcbio import broad
from bcbio.utils import curdir_tmpdir, file_exists, save_diskspace
from bcbio.distributed.transaction import file_transaction
from bcbio.variation.realign import has_aligned_reads


def gatk_recalibrate(dup_align_bam, ref_file, config, snp_file=None):
    """Perform a GATK recalibration of the sorted aligned BAM,
    producing recalibrated BAM.
    """
    broad_runner = broad.runner_from_config(config)
    platform = config["algorithm"]["platform"]
    broad_runner.run_fn("picard_index_ref", ref_file)
    recal_file = _gatk_count_covariates(broad_runner, dup_align_bam, \
                                        ref_file, platform, snp_file)
    recal_bam = _gatk_table_recalibrate(broad_runner, dup_align_bam, ref_file, \
                                        recal_file, platform)
    broad_runner.run_fn("picard_index", recal_bam)

    return recal_bam


def _gatk_table_recalibrate(broad_runner, dup_align_bam, ref_file, recal_file, platform):
    """Step 2 of GATK recalibration -- use covariates to re-write output file.
    """
    out_file = "%s-gatkrecal.bam" % os.path.splitext(dup_align_bam)[0]
    if file_exists(out_file):
        return out_file

    if _recal_available(recal_file):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_file) as tx_out_file:
                params = ["-T", "TableRecalibration",
                          "-recalFile", recal_file,
                          "-R", ref_file,
                          "-I", dup_align_bam,
                          "--out", tx_out_file,
                          "-baq",  "RECALCULATE",
                          "-l", "INFO",
                          "-U",
                          "-OQ",
                          "--default_platform", platform,
                          ]
                broad_runner.run_gatk(params, tmp_dir)
    else:
        shutil.copy(dup_align_bam, out_file)

    return out_file


def _recal_available(recal_file):
    """Determine if it's possible to do a recalibration; do we have data?
    """
    if os.path.exists(recal_file):
        with open(recal_file) as in_handle:
            while 1:
                line = in_handle.next()
                if not line.startswith("#"):
                    break

            test_line = in_handle.next()
            if test_line and not test_line.startswith("EOF"):
                return True

    return False


def _gatk_count_covariates(broad_runner, dup_align_bam, ref_file, platform,
        snp_file):
    """Step 1 of GATK recalibration process -- counting covariates.
    """
    out_file = "%s.recal" % os.path.splitext(dup_align_bam)[0]
    if not file_exists(out_file):
        if has_aligned_reads(dup_align_bam):
            with curdir_tmpdir() as tmp_dir:
                with file_transaction(out_file) as tx_out_file:
                    params = ["-T", "CountCovariates",
                              "-cov", "ReadGroupCovariate",
                              "-cov", "QualityScoreCovariate",
                              "-cov", "CycleCovariate",
                              "-cov", "DinucCovariate",
                              "-recalFile", tx_out_file,
                              "-I", dup_align_bam,
                              "-R", ref_file,
                              "-l", "INFO",
                              "-U",
                              "-OQ",
                              "--default_platform", platform,
                              ]
                    if snp_file:
                        params += ["--knownSites", snp_file]
                    broad_runner.run_gatk(params, tmp_dir)
        else:
            with open(out_file, "w") as out_handle:
                out_handle.write("# No aligned reads")
    return out_file
