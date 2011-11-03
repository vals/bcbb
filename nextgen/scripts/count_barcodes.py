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
import sys
from Bio import pairwise2
from optparse import OptionParser
import yaml
from operator import itemgetter

from bcbio.solexa import INDEX_LOOKUP


def main(fastq, run_info_file, lane, out_file, length, offset, mismatch, verbose):
    # TODO: Check run_info against INDEX_LOOKUP

    # Collect counts for all observed barcodes
    bcodes = {}
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

    # Seperate out the most common barcodes
    total = sum(bcodes.itervalues())
    bc_out = sorted(bcodes.iteritems(), key=itemgetter(1), reverse=True)

    # Splitting away those with less than 2% seem like a reasonable thing to do
    # (according to one single data point!)

    # Check 1 mismatch against most common

    # Check 2 mismatch against most common

    # ...

    # if run_info_file != None:
    #     matched_bc_grouping = \
    #     match_against_run_info(bcodes, run_info_file, mismatch, lane)
    # else:
    #     matched_bc_grouping = approximate_matching(bcodes, \
    #                                 list(set(INDEX_LOOKUP.values())), mismatch)

    # if verbose:
    #     print yaml.dump(matched_bc_grouping, width=70)

    if not out_file:
        out_file = fastq.split(".txt")[0] + "_barcodes.yaml"

    with open(out_file, "w+") as out_handle:
        out_handle.writelines("%s %f \n" % (bc, float(num) / float(total)) for bc, num in bc_out)
        #yaml.dump(bc_out, out_handle, width=70)

    print


def match_against_run_info(bcodes, run_info_file, mismatch, lane):
    given_bcodes = []
    with open(run_info_file) as in_handle:
        run_info = yaml.load(in_handle)
        given_bcodes += [bc["sequence"] for bc in run_info[lane - 1]["multiplex"]]

    return approximate_matching(bcodes, given_bcodes, mismatch)


def approximate_matching(bcodes, given_bcodes, mismatch):
    """Returns a dectionary with matched barcodes along with info.
    """
    matched_bc_grouping = {}
    found_bcodes = set()
    number = dict(matched=0., unmatched=0.)

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
                    matched_bc_grouping[bc_given] = {"variants": [], \
                                                        "count": 0}
                matched_bc_grouping[bc_given]["variants"].append(bc)
                matched_bc_grouping[bc_given]["count"] += count
                found_bcodes.add(bc)
                number["matched"] += count

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

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("--lane", dest="lane", default=0)
    parser.add_option("-o", "--out_file", dest="out_file", default=None)
    parser.add_option("-l", "--length", dest="length", default=6)
    parser.add_option("-b", "--back", dest="offset", default=0)
    parser.add_option("-m", "--mismatch", dest="mismatch", default=1)
    parser.add_option("-v", "--verbose", dest="verbose", default=False, \
                                                        action="store_true")
    options, args = parser.parse_args()
    if len(args) == 1:
        fastq, = args
        run_info = None
    elif len(args) == 2:
        fastq, run_info = args
    else:
        print __doc__
        sys.exit()

    main(fastq, run_info, int(options.lane), options.out_file, int(options.length), \
            int(options.offset), int(options.mismatch), options.verbose)
