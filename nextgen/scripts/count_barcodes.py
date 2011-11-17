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
    length, offset, mismatch, verbose, cutoff, dry_run, mode):
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
    print("Total after read:\t" + str(total))
    bc_matched = []
    for bc, num in bcodes.iteritems():
        if float(num) / total >= cutoff:
            bc_matched.append(bc)

    if mode == "demultiplex":
        # Check with mismatch against most common, split fastq files.
        bc_grouping = match_and_split(fastq, dict(bcodes), \
                                bc_matched, mismatch, offset, length, dry_run)
    elif mode == "count":
        bc_grouping = match_and_count(dict(bcodes), bc_matched, mismatch)

    if not out_file:
        out_file = fastq.split(".txt")[0] + "_barcodes.yaml"

    with open(out_file, "w+") as out_handle:
        yaml.dump(bc_grouping.__dict__, out_handle, width=70)


class BarcodeGrouping(object):
    """Stores and represents the groupings of matched barcodes.
    """
    def __init__(self):
        self.matched = dict()
        self.unmatched = dict()

    def add_illumina_indexes(self):
        """Attaches matched illumina barcodes to barcode groups
        """
        for bc, matches in self.matched.items():
            for illumina_index, illumina_bc in INDEX_LOOKUP.items():
                if illumina_bc == bc:
                    if "indexes" not in matches:
                        matches["indexes"] = []
                    matches["indexes"].append(illumina_index)

    def add_unmatched_barcodes(self, bcodes, found_bcodes):
        """Appends a subdictionary with the barcodes not considered to match
        the primary barcodes.
        """
        self.unmatched = dict((code, dict(count=bcodes[code], \
        variants=[code])) for code in set(bcodes) - found_bcodes)

# def match_against_run_info(bcodes, run_info_file, mismatch, lane):
#     given_bcodes = []
#     with open(run_info_file) as in_handle:
#         run_info = yaml.load(in_handle)
#         given_bcodes += [bc["sequence"] for \
#         bc in run_info[lane - 1]["multiplex"]]

#     return approximate_matching(bcodes, given_bcodes, mismatch)


def match_and_count(bcodes, given_bcodes, mismatch):
    """Returns a dictionary with matched barcodes along with info.
    """
    bc_grouping = BarcodeGrouping()
    number = dict(matched=0., unmatched=0.)
    found_bcodes = set()

    assert mismatch >= 0, "Amount of mismatch cannot be negative."
    if mismatch == 0:
        for bc, count in bcodes.items():
            if bc in given_bcodes:
                if bc not in bc_grouping.matched:
                    bc_grouping.matched[bc] = {"variants": [bc], "count": 0}

                bc_grouping.matched[bc]["count"] += count
                found_bcodes.add(bc)
                number["matched"] += count

    else:
        for bc, count in bcodes.items():
            for bc_given in given_bcodes:
                current_mismatch = bc_mismatch(bc, bc_given)
                if current_mismatch <= mismatch:
                    if bc_given not in bc_grouping.matched:
                        bc_grouping.matched[bc_given] = {"variants": [], \
                                                            "count": 0}
                    if bc not in bc_grouping.matched[bc_given]["variants"]:
                        bc_grouping.matched[bc_given]["variants"].append(bc)

                    bc_grouping.matched[bc_given]["count"] += count
                    found_bcodes.add(bc)
                    number["matched"] += count
                    break

    bc_grouping.add_illumina_indexes()
    bc_grouping.add_unmatched_barcodes(bcodes, found_bcodes)

    total = sum(bcodes.values())

    number["unmatched"] = \
    float(sum(value["count"] for value in bc_grouping.unmatched.values()))
    percentage = 100. * number["matched"] / sum(number.values())
    print("Hard numbers:\t\t" + str(number))
    print("Sum:\t\t\t" + str(sum(number.values())))
    print("Total:\t\t\t" + str(total))
    print("Percentage matched:\t%.3f%%" % percentage)
    return bc_grouping


def match_and_split(fastq, bcodes, given_bcodes, \
    mismatch, offset, length, dry_run):
    """Matches barcodes in the fastq file 'fastq' and prints out the lines
    corresponding to the matched barcode groupes in to seperate files.
    When done it returns a dictionary with matched barcodes along with info.
    """
    bc_grouping = BarcodeGrouping()
    found_bcodes = set()
    number = dict(matched=0., unmatched=0.)
    out_format = fastq.split(".")[0] + "_out/out_--b--_--r--_fastq.txt"
    out_writer = output_to_fastq(out_format)

    assert mismatch >= 0, "Amount of mismatch cannot be negative."
    handle = open(fastq)
    if mismatch == 0:
        for title, sequence, quality in FastqGeneralIterator(handle):
            bc = sequence[-(offset + 1 + length):-(offset + 1)].strip()
            if bc in given_bcodes:
                if bc not in bc_grouping.matched:
                    bc_grouping.matched[bc] = {"variants": [], "count": 0}
                    if bc not in bc_grouping.matched[bc]["variants"]:
                        bc_grouping.matched[bc]["variants"].append(bc)

                    bc_grouping.matched[bc]["count"] += 1
                    found_bcodes.add(bc)
                    number["matched"] += 1

                    if not dry_run:
                        out_writer(bc, title, sequence, quality, \
                        None, None, None)
    else:
        for title, sequence, quality in FastqGeneralIterator(handle):
            bc = sequence[-(offset + 1 + length):-(offset + 1)].strip()
            for bc_given in given_bcodes:
                cur_mismatch = bc_mismatch(bc, bc_given)
                if cur_mismatch <= mismatch:
                    if bc_given not in bc_grouping.matched:
                        bc_grouping.matched[bc_given] = {"variants": [], \
                                                            "count": 0}
                    if bc not in bc_grouping.matched[bc_given]["variants"]:
                        bc_grouping.matched[bc_given]["variants"].append(bc)

                    bc_grouping.matched[bc_given]["count"] += 1
                    found_bcodes.add(bc)
                    number["matched"] += 1

                    if not dry_run:
                        out_writer(bc, title, sequence, quality, \
                        None, None, None)

    bc_grouping.add_illumina_indexes()
    bc_grouping.add_unmatched_barcodes(bcodes, found_bcodes)

    number["unmatched"] = float(sum(bc_grouping.unmatched.values()))
    percentage = 100. * number["matched"] / sum(number.values())
    print("Hard numbers:" + number)
    print("Percentage matched: %.3f%%" % percentage)
    return bc_grouping


def bc_mismatch(bc, bc_given):
    """Calculate and return the number of mismatches between two barcodes.
    """
    aligns = pairwise2.align.globalms(bc, bc_given,
                5.0, -4.0, -9.0, -0.5, one_alignment_only=True)
    bc_aligned, bc_g_aligned = aligns[0][:2]
    matches = sum(1 for i, base in enumerate(bc_aligned) \
                                        if base == bc_g_aligned[i])
    gaps = bc_aligned.count("-")
    cur_mismatch = len(bc) - matches + gaps
    return cur_mismatch


def _write_to_handles(name, seq, qual, fname, out_handles):
    """Defines format and location to write the fastq triple.
    """
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
            read2name = \
            output_base.replace("--r--", "2").replace("--b--", barcode)
            _write_to_handles(name2, seq2, qual2, read2name, out_handles)

    return write_reads


def compare_run_info_and_index_lookup(run_info_file):
    """A simple check to see that the barcodes in the given run_info matches
    barcodes specified by Illumina documentation.
    """
    known = INDEX_LOOKUP.values()
    unknown = set()
    run_info = yaml.load(run_info_file)
    for lanes in run_info:
        for b_ids in lanes["multiplex"]:
            bc = b_ids["sequence"]
            if bc not in known:
                unknown.add(bc)
    
    return unknown

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
    parser.add_option("-n", "--dryrun", dest="dry_run", default=False, \
                                                        action="store_true")
    parser.add_option("--mode", dest="mode", default="demultiplex")
    options, args = parser.parse_args()
    if len(args) == 1:
        fastq, = args
        run_info = None
    elif len(args) == 2:
        fastq, run_info = args
    else:
        print __doc__
        sys.exit()

    # import cProfile
    # cProfile.run("main(fastq, run_info, int(options.lane), options.out_file,\
    # int(options.length), int(options.offset), int(options.mismatch), \
    # options.verbose, float(options.cutoff), options.dry_run, options.mode)",\
    # sort='cumulative')
    main(fastq, run_info, int(options.lane), options.out_file, \
    int(options.length), int(options.offset), int(options.mismatch), \
    options.verbose, float(options.cutoff), options.dry_run, options.mode)
