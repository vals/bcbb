"""Compares the Illumina HiSeq bar codes in a fastq file with a given
samplesheet file or with bar codes given in Illumina documentation.
# This script uses a few assumptions about how the barcoding works.
# Since the aim is to look at Illumina HiSeq barcodes, the assumption is that
# the barcode is located at the 3' end of read 1 for each paired-end read.

Usage:
    count_barcodes.py <fastq file> [<run info yaml file>]
        -l, --length
        -b, --back
        -m, --mismatch
"""
import sys
from Bio import pairwise2
from optparse import OptionParser
import yaml

from bcbio.solexa import INDEX_LOOKUP


def main(fastq, run_info_file, length, offset, mismatch):
    bcodes = {}  # collect counts for all observed barcodes

    mismatch = 1
    length = 6

    last_was_header = False
    for line in open(fastq):
        if not last_was_header:
            if line[0] == "@":
                last_was_header = True
            continue

        bcode = line[-(offset + 1 + length):-(offset + 1)].strip()
        if bcode in bcodes:
            bcodes[bcode] += 1
        else:
            bcodes[bcode] = 1
        if last_was_header:
            last_was_header = False

    if run_info_file != None:
        with open(run_info_file) as in_handle:
            run_info = yaml.load(in_handle)
            given_bcodes = [bc["sequence"] for bc in run_info[0]["multiplex"]]

        matched_bc_grouping = approximate_matching(bcodes, given_bcodes, mismatch)

    else:
        matched_bc_grouping = approximate_matching(bcodes, list(set(INDEX_LOOKUP.values())), mismatch)

    print matched_bc_grouping
    with open("dumped.yaml", "w") as out_handle:
        yaml.dump(matched_bc_grouping, out_handle, width=70)


def approximate_matching(bcodes, given_bcodes, mismatch):
    """ Returns a dectionary with found barcodes along with info
    """
    matched_bc_grouping = {}
    found_bcodes = set()

    assert mismatch >= 0, "Amount of mismatch cannot be negative."
    for bc, count in bcodes.items():
        for bc_given in given_bcodes:
            aligns = pairwise2.align.globalms(bc, bc_given,
                        5.0, -4.0, -9.0, -0.5, one_alignment_only=True)
            bc_aligned, bc_g_aligned = aligns[0][:2]
            matches = sum(1 for i, base in enumerate(bc_aligned) \
                                                if base == bc_g_aligned[i])
            gaps = bc_aligned.count("-")
            cur_mismatch = len(bc) - matches + gaps

            if cur_mismatch <= mismatch:
                if bc_given not in matched_bc_grouping:
                    matched_bc_grouping[bc_given] = {"variants": [], "count": 0}
                matched_bc_grouping[bc_given]["variants"].append(bc)
                matched_bc_grouping[bc_given]["count"] += count
                found_bcodes.add(bc)

    for bc, matches in matched_bc_grouping.items():
        for illumina_index, illumina_bc in INDEX_LOOKUP.items():
            if illumina_bc == bc:
                if "indexes" not in matches:
                    matches["indexes"] = []
                matches["indexes"].append(illumina_index)

    matched_bc_grouping["unmatched"] = list(set(bcodes) - found_bcodes)

    return matched_bc_grouping

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-l", "--length", dest="length", default=6)
    parser.add_option("-b", "--back", dest="offset", default=0)
    parser.add_option("-m", "--mismatch", dest="mismatch", default=1)
    options, args = parser.parse_args()
    if len(args) == 1:
        fastq, = args
        run_info = None
    elif len(args) == 2:
        fastq, run_info = args
    else:
        print __doc__
        sys.exit()
    main(fastq, run_info, int(options.length), int(options.offset),
                                                int(options.mismatch))
