"""Multiprocessing ready entry points for sample analysis.
"""
from bcbio import utils
from bcbio.pipeline import sample, lane, shared, variation
from bcbio.variation import realign, genotype


@utils.map_wrap
def process_lane(*args):
    return lane.process_lane(*args)


@utils.map_wrap
def remove_contaminants(*args):
    return lane.remove_contaminants(*args)


@utils.map_wrap
def process_alignment(*args):
    return lane.process_alignment(*args)


@utils.map_wrap
def mark_duplicates_sample(*args):
    return sample.mark_duplicates_sample(*args)


@utils.map_wrap
def merge_sample(*args):
    return sample.merge_sample(*args)


@utils.map_wrap
def recalibrate_sample(*args):
    return sample.recalibrate_sample(*args)


@utils.map_wrap
def realign_sample(*args):
    return realign.realign_sample(*args)


@utils.map_wrap
def screen_sample_contaminants(*args):
    return sample.screen_sample_contaminants(*args)


@utils.map_wrap
def process_sample(*args):
    return sample.process_sample(*args)


@utils.map_wrap
def generate_bigwig(*args):
    return sample.generate_bigwig(*args)


@utils.map_wrap
def combine_bam(*args):
    return shared.combine_bam(*args)


@utils.map_wrap
def variantcall_sample(*args):
    return genotype.variantcall_sample(*args)


@utils.map_wrap
def combine_variant_files(*args):
    return genotype.combine_variant_files(*args)

@utils.map_wrap
def detect_sv(*args):
    return variation.detect_sv(*args)
