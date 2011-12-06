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

    matched:
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
    unmatched:
        AAAAAA:
            count: 3
            variants: [AAAAAA]
        AAAAAC:
            count: 1
            variants: [AAAAAC]
        AAAAAG:
            count: 2
            variants: [AAAAAG]
        AAAAAT:
            count: 1
            variants: [AAAAAT]
"""
from __future__ import with_statement
import collections
from itertools import izip
import os
from optparse import OptionParser
import sys
import shutil
import yaml

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from bcbio.solexa import INDEX_LOOKUP


def main(fastq, run_info_file, lane, out_file,
    length, offset, mismatch, verbose, cutoff, dry_run, mode):
    
    bc_matched = []
    if run_info_file:
        bc_matched = _get_run_info_barcodes(run_info_file,lane)
        compare_run_info_and_index_lookup(bc_matched)
        
    out_format = fastq.split(".")[0] + "_out/out_--b--_--r--_fastq.txt"
    if mode == "demultiplex":
        out_writer = output_to_fastq(out_format)

    # Collect counts for all observed barcodes
    bcodes = collections.defaultdict(int)
    in_handle = open(fastq)
    if offset == 0:
        minus_offset = None
    else:
        minus_offset = -offset
    for title, sequence, quality in FastqGeneralIterator(in_handle):
        bcode = sequence[-(offset + length):minus_offset].strip()
        if mode == "demultiplex":
            out_writer(bcode, title, sequence, quality, None, None, None)
        bcodes[bcode] += 1
        # Count number not A in the last one

    # Automatically determine which barcodes to use unless a run_info file was specified
    if not run_info_file:
        bc_matched = _get_common_barcodes(bcodes,cutoff)
        
    # Get the barcode statistics
    # bcm_nums, bcm_parts = _get_barcode_statistics(bc_matched,bcodes)
    
    # TODO: Match matched bcs against each other.

    if mode == "demultiplex":
        # Check with mismatch against most common, split fastq files.
        bc_grouping = \
        match_and_merge(dict(bcodes), bc_matched, mismatch, out_format)

    elif mode == "count":
        bc_grouping = match_and_count(dict(bcodes), bc_matched, mismatch)

    if not out_file:
        out_file = fastq.split(".txt")[0] + "_barcodes.yaml"

    with open(out_file, "w+") as out_handle:
        yaml.dump(bc_grouping.__dict__, out_handle, width=70)

def _get_barcode_statistics(exp_barcodes,obs_barcodes):
    """Get the prevalence of a particular barcode, as a number and a frequency"""
    total = float(sum(obs_barcodes.itervalues()))
    print("Total after read:\t%.0f" % (total,))
    bcm_nums = []
    bcm_parts = []
    for bc in exp_barcodes:
        num = obs_barcodes.get(bc,0)
        part = float(num)/total
        bcm_nums.append(num)
        bcm_parts.append(part)
    return (bcm_nums,bcm_parts)
 
def _get_common_barcodes(bcodes,cutoff):
    """Get the barcodes that occur at a frequency above cutoff"""
    # Seperate out the most common barcodes
    total = float(sum(bcodes.itervalues()))
    bc_matched = []
    for bc, num in bcodes.iteritems():
        part = float(num) / total
        if part >= cutoff:
            bc_matched.append(bc)
            
    return bc_matched

def _get_run_info_barcodes(run_info_file,lane=0):
    """Extract the barcodes to demultiplex against from run_info file"""

    barcodes = {}
    with open(run_info_file) as fh:
        run_info = yaml.load(fh)
        for lane_info in run_info:
            if (lane != 0 and int(lane_info.get("lane",0)) != lane):
                continue
            for bc in lane_info.get("multiplex",{}):
                barcodes[bc.get("sequence","")] = 1
    return barcodes.keys()


class BarcodeGrouping(object):
    """Stores and represents the groupings of matched barcodes.
    """
    def __init__(self):
        self.matched = dict()
        self.unmatched = dict()

    def add_illumina_indexes(self):
        """Attaches matched illumina barcodes to barcode groups
        """
        for bc, matches in self.matched.iteritems():
            for illumina_index, illumina_bc in INDEX_LOOKUP.iteritems():
                if illumina_bc == bc:
                    if "indexes" not in matches:
                        matches["indexes"] = []
                    matches["indexes"].append(illumina_index)

    def add_given_matched_barcodes(self, given_bcodes):
        """Appends the barcodes we wish to match against, with count 0.
        """
        for bcode in given_bcodes:
            self.matched[bcode] = dict(count=0, variants=[bcode])

    def add_unmatched_barcodes(self, bcodes, found_bcodes):
        """Appends a subdictionary with the barcodes not considered to match
        the primary barcodes.
        """
        self.unmatched = dict((code, dict(count=bcodes[code], \
        variants=[code])) for code in set(bcodes) - found_bcodes)

    def load_from_yaml(self, yaml_file):
        """Loads the matched and unmatched elements from the given yaml file.
        """
        with open(yaml_file, "r") as in_handle:
            bc_grouping_dict = yaml.load(in_handle)
            self.matched = bc_grouping_dict["matched"]
            self.unmatched = bc_grouping_dict["unmatched"]

    def handle_Ns(self):
        """Goes through the matched barcodes, looking which ones contains 'N'
        and tries to replace that with a matched barcode which does not
        contain 'N'. If that is not possible, the barcode will be moved to the
        unmatched catergory.
        """
        for barcode, info in self.matched.items():
            if 'N' in barcode:
                for matched_barcode in info["variants"]:
                    if 'N' not in matched_barcode:
                        self.matched[matched_barcode] = info
                        break
                else:
                    self.unmatched[barcode] = info

                del self.matched[barcode]

# def match_against_run_info(bcodes, run_info_file, mismatch, lane):
#     given_bcodes = []
#     with open(run_info_file) as in_handle:
#         run_info = yaml.load(in_handle)
#         given_bcodes += [bc["sequence"] for \
#         bc in run_info[lane - 1]["multiplex"]]

#     return approximate_matching(bcodes, given_bcodes, mismatch)


def _match_barcodes(bcode,given_bcodes,mismatch,masked=False):
    """Logic for matching a barcode against the given barcodes"""
    
    # First check for perfect matches
    matched = ""
    if bcode in given_bcodes:
        matched = bc    
    # If a perfect match could not be found, do a finer matching but only if we allow mismatches or the given barcodes contain masked positions
    elif mismatch > 0 or masked:
        for gbc in given_bcodes:
            current_mismatch = bc_mismatched(bcode,gbc,mismatch)
            if current_mismatch <= mismatch:
                matched = gbc
                break
    return matched

def match_and_count(bcodes, given_bcodes, mismatch):
    """Returns a dictionary with matched barcodes along with info.
    """
    bc_grouping = BarcodeGrouping()
    number = dict(matched=0., unmatched=0.)
    found_bcodes = set()

    assert mismatch >= 0, "Amount of mismatch cannot be negative."
    # Set a flag indicating whether the given barcodes contain masked nucleotides
    masked = False
    for gbc in given_bcodes:
        if 'N' in gbc:
            masked = True
            break
        
    for bc, count in bcodes.iteritems():
        match = _match_barcodes(bc,given_bcodes,mismatch,masked)
        if len(match) > 0:
            if match not in bc_grouping.matched:
                bc_grouping.matched[match] = {"variants": [], "count": 0}                                      
            if bc not in bc_grouping.matched[match]["variants"]:
                bc_grouping.matched[match]["variants"].append(bc)

            bc_grouping.matched[match]["count"] += count
            found_bcodes.add(bc)
            number["matched"] += count
    
    bc_grouping.add_unmatched_barcodes(bcodes, found_bcodes)
    bc_grouping.handle_Ns()
    bc_grouping.add_illumina_indexes()

    total = sum(bcodes.itervalues())

    number["unmatched"] = \
    float(sum(value["count"] for value in bc_grouping.unmatched.itervalues()))
    percentage = 100. * number["matched"] / sum(number.values())
    print("Hard numbers:\t\t" + str(number))
    print("Sum:\t\t\t" + str(sum(number.values())))
    print("Total:\t\t\t" + str(total))
    print("Percentage matched:\t%.3f%%" % percentage)
    return bc_grouping


def match_and_merge(bcodes, given_bcodes, mismatch, format):
    """Matches barcodes and merges the corresponding splitted fastq files in to
    the large fastq files.
    Returns a dictionary with matched barcodes along with info.
    """
    bc_grouping = BarcodeGrouping()
    bc_grouping.add_given_matched_barcodes(given_bcodes)
    found_bcodes = set()
    merger = FileMerger()

    old_unmatched = \
    format.replace("--r--", "1").replace("--b--", "unmatched")
    if os.path.exists(old_unmatched):
        os.remove(old_unmatched)

    assert mismatch >= 0, "Amount of mismatch cannot be negative."
    
    # Set a flag indicating whether the given barcodes contain masked nucleotides
    masked = False
    for gbc in given_bcodes:
        if 'N' in gbc:
            masked = True
            break
    
    for bc, count in bcodes.iteritems():
        match = _match_barcodes(bc,given_bcodes,mismatch,masked)
        if len(match) > 0:
            if match not in bc_grouping.matched:
                bc_grouping.matched[match] = {"variants": [], "count": 0}                                      
            if bc not in bc_grouping.matched[match]["variants"]:
                bc_grouping.matched[match]["variants"].append(bc)

            bc_grouping.matched[match]["count"] += count
            found_bcodes.add(bc)
            number["matched"] += count
    
            if match != bc:
                merge_matched_files(match,bc,format,merger)
        else:
            merge_matched_files("unmatched",bc,format,merger)
    
    merger.close()

    bc_grouping.add_unmatched_barcodes(bcodes, found_bcodes)
    bc_grouping.handle_Ns()
    bc_grouping.add_illumina_indexes()

    total = sum(bcodes.itervalues())

    number = dict(matched=0., unmatched=0.)
    number["matched"] = \
    float(sum(value["count"] for value in bc_grouping.matched.itervalues()))
    number["unmatched"] = \
    float(sum(value["count"] for value in bc_grouping.unmatched.itervalues()))
    percentage = 100. * number["matched"] / sum(number.values())
    print("Hard numbers:\t\t" + str(number))
    print("Sum:\t\t\t%.0f" % (sum(number.values()),))
    print("Total:\t\t\t%d" % (total,))
    print("Percentage matched:\t%.3f%%" % percentage)
    return bc_grouping


def merge_matched_files(primary_bc, matched_bc, format, merger):
    """Merges the fastq file corresponding to matched_bc in to the fastq file
    corresponding to primary_bc.
    """
    matched_file = format.replace("--r--", "1").replace("--b--", matched_bc)
    primary_file = format.replace("--r--", "1").replace("--b--", primary_bc)

    merger(matched_file, primary_file)

    os.remove(matched_file)


def bc_mismatched(bc, bc_given, mismatch):
    """Calculate the number of mismatches between two barcodes until
    the given mismatch allowance is reached.
    Return number of mismatches if there are less than the allowance, otherwise
    return (mismatch + 1).
    """
    current_mismatch = 0
    for c1, c2 in izip(bc, bc_given):
        if c1 != c2 and c2 != 'N':
            current_mismatch += 1
            if current_mismatch > mismatch:
                break

    return current_mismatch


class FileMerger(dict):
    """Returns an object which handles merging files.
    """
    def __call__(self, source_file, target_file):
        self.append_to_handles(source_file, target_file)

    def append_to_handles(self, source_file, target_file):
        """Appends the contents of source_file to the target_file, storing
        file objects in the target_handles dictionary to speed up IO.
        """
        try:
            out_handle = self[target_file]
        except KeyError:
            out_handle = open(target_file, "a")
            self[target_file] = out_handle

        in_handle = open(source_file, "r")
        shutil.copyfileobj(in_handle, out_handle)
        in_handle.close()

    def close(self):
        """Close all the opened files.
        """
        for handle in self.itervalues():
            handle.close()


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


def compare_run_info_and_index_lookup(run_info_barcodes):
    """A simple check to see that the barcodes in the given run_info matches
    barcodes specified by Illumina documentation.
    """
    known = INDEX_LOOKUP.values()
    unknown = set()
    for bc in run_info_barcodes:
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
