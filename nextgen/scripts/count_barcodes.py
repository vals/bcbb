"""Compares the Illumina HiSeq bar codes in a fastq file with a given
samplesheet file or with bar codes given in Illumina documentation.
# This script uses a few assumptions about how the barcoding works.
# Since the aim is to look at Illumina HiSeq barcodes, the assumption is that
# the barcode is located at the 3' end of read 1 for each paired-end read.

Usage:
    count_barcodes.py <fastq file> [<run info yaml file>]
        -o, --out_file <name of yaml file to be written (if not given, the out
                        file will be named as the fastq file, though with
                        .txt replaced by _barcodes.yaml)>
        -l, --length <number of characters in the bar code (default is 6)>
        -b, --back <number of steps back from the end of each line where the
                    bar code ends, 1 for Illumina (defualt is 0).>
        -m, --mismatch <number of erroneous characters to consider when
                        matching bar codes (default is 1).>
        -v, --verbose (sets the script to print out what is written to
                        the file)

    If called only with the fastq file, the bar codes will be matched to and
    grouped with bar codes from the Illumina documentation.
    If also the <run info yaml file> is given, the bar codes in the fastq file
    will be matched to the bar codes in the run info file. Index names will
    still be extracted from the Illumina documentation table.

Example:
    count_barcodes.py 1_110106_FC70BUKAAXX_1_fastq.txt -v -b 1 -m 0
    will create a file named "1_110106_FC70BUKAAXX_1_fastq_barcodes.yaml"
    where the bar codes matched to bar codes in the Illumina documentation
    will be listed. Every bar code node in the yaml has a field with count,
    the number (of variations) of the bar code found; indexes, the names of
    the Illumina bar codes. And a list of the variations of the bar code
    considered to match the bar code.
    The yaml will also have a node "unmatched", listing the bar codes which
    could not be matched with the given mismatch setting, along with the count
    of each of those bar codes.

    Essentially, it looks somewhat like this:
        TGACCA:
          count: 902
          indexes: [rpi4, INDEX4, IDX4, IN4, '4', R4, index4, r4, idx4, RPI4,
          in4]
          variants: [AGACCA, GGACCA, TCACCA, TGACCG, TAACCA, TGACCA, TGATCA,
          TGACCC, TTACCA, TGACCT, TGCCCA]
        TTAGGC:
          count: 1066
          indexes: [rpi3, IDX3, IN3, '3', R3, INDEX3, index3, r3, idx3, RPI3,
          in3]
          variants: [TTACGC, TTAAGC, TAAGGC, TTAGGC, TTAGGA, TTCGGC, TTAGAC,
          CTAGGC, TTAGTC, TGAGGC, TTATGC, ATAGGC, TCAGGC]
        unmatched: {AAAGTC: 1, AACAGT: 1, AATCAC: 1, ACAAGT: 1, ACACCA: 1,
          ACAGGA: 1, ACAGGC: 2, ACATGC: 1, ACCGAG: 1, ACCGCG: 1, ACGATC: 2,
          ACTCGG: 1, ACTCTC: 1}
"""
import os
import sys
from Bio import pairwise2
from optparse import OptionParser
import yaml
import collections

from bcbio.solexa import INDEX_LOOKUP
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def main(fastq, run_info_file, lane, out_file,
    length, offset, mismatch, verbose, cutoff):
    if run_info_file:
        compare_run_info_and_index_lookup(run_info_file)

    # Collect counts for all observed barcodes
    bcodes = collections.defaultdict(int)
    in_handle = open(fastq)
    for _, sequence, _ in FastqGeneralIterator(in_handle):
        bcode = sequence[-(offset + 1 + length):-(offset + 1)].strip()
        bcodes[bcode] += 1

    # Seperate out the most common barcodes
    total = float(sum(bcodes.itervalues()))
    bc_matched = []
    for bc, num in bcodes.iteritems():
        if float(num) / total >= cutoff:
            bc_matched.append(bc)

    # Check with mismatch against most common

    matched_bc_grouping = approximate_matching(fastq, dict(bcodes), \
                                                bc_matched, mismatch, offset, length)

    if not out_file:
        out_file = fastq.split(".txt")[0] + "_barcodes.yaml"

    with open(out_file, "w+") as out_handle:
        yaml.dump(matched_bc_grouping, out_handle, width=70)


def match_against_run_info(bcodes, run_info_file, mismatch, lane):
    given_bcodes = []
    with open(run_info_file) as in_handle:
        run_info = yaml.load(in_handle)
        given_bcodes += [bc["sequence"] for bc in run_info[lane - 1]["multiplex"]]

    return approximate_matching(bcodes, given_bcodes, mismatch)


def approximate_matching(fastq, bcodes, given_bcodes, mismatch, offset, length):
    """Returns a dectionary with matched barcodes along with info.
    """
    # TODO: Splitting - need some way to copy the entire name / sequence / quality
    # three lines from where the barcode was collected in to a seperate file
    # for each barcode.

    # Strategy:
    # Count -> Sort -> Split -> Match all against most common with mistmatch

    # We will have the barcode count after the split. So we don't need to
    # count when matching. So in the matching step we could iterate over the file,
    # storing the three relevant lines in memory before being written to the
    # correct file.
    # Use FastqGeneralIterator to get the name / sequence / quality triple!
    matched_bc_grouping = {}
    found_bcodes = set()
    number = dict(matched=0., unmatched=0.)
    out_format = "out/out_--b--_--r--_fastq.txt"
    out_writer = output_to_fastq(out_format)

    assert mismatch >= 0, "Amount of mismatch cannot be negative."
    handle = open(fastq)
    for title, sequence, quality in FastqGeneralIterator(handle):
    #for bc, count in bcodes.items():
        bc = sequence[-(offset + 1 + length):-(offset + 1)].strip()
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
                    matched_bc_grouping[bc_given] = {"variants": [], \
                                                        "count": 0}
                count = bcodes[bc]
                if bc not in matched_bc_grouping[bc_given]["variants"]:
                    matched_bc_grouping[bc_given]["variants"].append(bc)

                matched_bc_grouping[bc_given]["count"] += count
                found_bcodes.add(bc)
                number["matched"] += count

                out_writer(bc, title, sequence, quality, None, None, None)

    for bc, matches in matched_bc_grouping.items():
        for illumina_index, illumina_bc in INDEX_LOOKUP.items():
            if illumina_bc == bc:
                if "indexes" not in matches:
                    matches["indexes"] = []
                matches["indexes"].append(illumina_index)

    matched_bc_grouping["unmatched"] = \
    dict((code, bcodes[code]) for code in set(bcodes) - found_bcodes)
    number["unmatched"] = float(sum(matched_bc_grouping["unmatched"].values()))
    percentage = 100. * number["matched"] / sum(number.values())
    print("Percentage matched: %.3f%%" % percentage)
    return matched_bc_grouping


def _write_to_handles(name, seq, qual, fname, out_handles):
    try:
        out_handle = out_handles[fname]
    except KeyError:
        out_handle = open(fname, "w")
        out_handles[fname] = out_handle
    out_handle.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))


def output_to_fastq(output_base):
    """Write a set of paired end reads as fastq, managing output handles.
    """
    work_dir = os.path.dirname(output_base)
    if not os.path.exists(work_dir) and work_dir:
        try:
            os.makedirs(work_dir)
        except OSError:
            assert os.path.isdir(work_dir)
    out_handles = dict()

    def write_reads(barcode, name1, seq1, qual1, name2, seq2, qual2):
        read1name = output_base.replace("--r--", "1").replace("--b--", barcode)
        _write_to_handles(name1, seq1, qual1, read1name, out_handles)
        if name2:
            read2name = output_base.replace("--r--", "2").replace("--b--", barcode)
            _write_to_handles(name2, seq2, qual2, read2name, out_handles)
    return write_reads


def compare_run_info_and_index_lookup(run_info):
    # TODO: Check run_info against INDEX_LOOKUP
    pass

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("--lane", dest="lane", default=0)
    parser.add_option("-o", "--out_file", dest="out_file", default=None)
    parser.add_option("-l", "--length", dest="length", default=6)
    parser.add_option("-b", "--back", dest="offset", default=0)
    parser.add_option("-m", "--mismatch", dest="mismatch", default=1)
    parser.add_option("-v", "--verbose", dest="verbose", default=False, \
                                                        action="store_true")
    parser.add_option("-c", "--cutoff", dest="cutoff", default=0.02)
    options, args = parser.parse_args()
    if len(args) == 1:
        fastq, = args
        run_info = None
    elif len(args) == 2:
        fastq, run_info = args
    else:
        print __doc__
        sys.exit()

    main(fastq, run_info, int(options.lane), options.out_file, \
            int(options.length), int(options.offset), int(options.mismatch), \
            options.verbose, float(options.cutoff))
