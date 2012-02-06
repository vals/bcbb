"""Performs demultiplexing or just gather barcode statistics of given input.

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
import re
import unittest
import glob

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from bcbio.solexa import INDEX_LOOKUP


def main(fastq1, fastq2, barcode_file, lane, out_file,
    length, offset, mismatch, verbose, cutoff, dry_run, mode, out_format):

    bc_matched = []
    if barcode_file:
        if barcode_file.split(".")[-1] == "yaml":
            bc_matched = _get_run_info_barcodes(barcode_file, lane)
            print(bc_matched)
            compare_run_info_and_index_lookup(bc_matched)
        elif barcode_file.split(".")[-1] == "cfg":
            barcodes = {}
            with open(barcode_file) as in_handle:
                for line in (l for l in in_handle if not l.startswith("#")):
                    name, seq = line.rstrip("\r\n").split()
                    barcodes[seq] = name
            bc_matched = barcodes.keys()

    # Collect counts for all observed barcodes
    bcodes = collections.defaultdict(int)
    in_handle = open(fastq1)

    if mode == "demultiplex":
        out_writer = FileWriter(out_format)

        if fastq2 is not None:
            in_handle2 = open(fastq2)

    old_fastq_header_re = re.compile(r"(?P<instrument>[\d\w-]*)"
    ":(?P<fc_lane>\d*):(?P<tile>\d*):(?P<x>\d*)"
    ":(?P<y>\d*)#(?P<multiplex_id>[\d\w]*)/(?P<pair>\d)")
    reformat_title = False
    if fastq2 is None:
        fastq1_iterator = FastqGeneralIterator(in_handle)
        title, sequence, quality = fastq1_iterator.next()
        _, bcode = trim_sequence(sequence, offset, length)

        if mode == "demultiplex":
            if old_fastq_header_re.match(title):
                reformat_title = True
            out_writer(bcode, title, sequence, quality, None, None, None)

        for title, sequence, quality in FastqGeneralIterator(in_handle):
            _, bcode = trim_sequence(sequence, offset, length)

            if mode == "demultiplex":
                out_writer(bcode, title, sequence, quality, None, None, None)

            bcodes[bcode] += 1

    else:
        fastqs_iterator = izip(FastqGeneralIterator(in_handle), FastqGeneralIterator(in_handle2))
        (title1, sequence1, quality1), (title2, sequence2, quality2) = \
        fastqs_iterator.next()
        _, bcode = trim_sequence(sequence1, offset, length)

        if mode == "demultiplex":
            if old_fastq_header_re.match(title1):
                reformat_title = True
            out_writer(bcode, title1, sequence1, quality1, title2, sequence2, quality2)

        for (title1, sequence1, quality1), (title2, sequence2, quality2) in \
        fastqs_iterator:
            _, bcode = trim_sequence(sequence1, offset, length)

            if mode == "demultiplex":
                out_writer(bcode, title1, sequence1, quality1, title2, sequence2, quality2)

            bcodes[bcode] += 1

    in_handle.close()

    if mode == "demultiplex":
        out_writer.close()
        if fastq2 is not None:
            in_handle2.close()

    # Automatically determine which barcodes to use unless a run_info file was specified
    if not barcode_file:
        bc_matched = _get_common_barcodes(bcodes, cutoff)

    if mode == "demultiplex":
        # Check with mismatch against most common, split fastq files.
        bc_grouping = match_and_merge(dict(bcodes), bc_matched, mismatch, \
        out_format, paired=bool(fastq2))

        if reformat_title:
            out_glob = out_format.replace("--b--", "*").replace("--r--", "1")
            out_files = [f for f in glob.glob(out_glob) if "unmatched" not in f]
            if fastq2 is not None:
                out_glob2 = out_format.replace("--b--", "*").replace("--r--", "2")
                out_files2 = [f for f in glob.glob(out_glob2) if "unmatched" not in f]

            if fastq2 is None:
                for out_file in out_files:
                    old_handle = open(out_file)
                    new_handle = open(out_file + ".tmp", "w")
                    for title, seq, qual in FastqGeneralIterator(old_handle):
                        trimmed_seq, barcode = trim_sequence(seq, offset, length)
                        trimmed_qual, _ = trim_sequence(qual, offset, length)
                        new_title = convert_old_fastq_header(title, barcode, \
                        old_fastq_header_re)
                        new_handle.write("@%s\n%s\n+\n%s\n" % \
                        (new_title, trimmed_seq, trimmed_qual))
                    old_handle.close()
                    new_handle.close()
                    shutil.move(out_file + ".tmp", out_file)

            else:
                for out_file1, out_file2 in zip(out_files, out_files2):
                    old_handle1 = open(out_file1)
                    old_handle2 = open(out_file2)
                    new_handle1 = open(out_file1 + ".tmp", "w")
                    new_handle2 = open(out_file2 + ".tmp", "w")
                    for (title, seq1, qual1), (_, seq2, qual2) in \
                    izip(FastqGeneralIterator(old_handle1), FastqGeneralIterator(old_handle2)):
                        trimmed_seq1, barcode = trim_sequence(seq1, offset, length)
                        trimmed_qual1, _ = trim_sequence(qual1, offset, length)
                        new_title1 = convert_old_fastq_header(title, barcode, \
                        old_fastq_header_re, pair="1")
                        new_title2 = convert_old_fastq_header(title, barcode, \
                        old_fastq_header_re, pair="2")
                        new_handle1.write("@%s\n%s\n+\n%s\n" % \
                        (new_title1, trimmed_seq1, trimmed_qual1))
                        new_handle2.write("@%s\n%s\n+\n%s\n" % \
                        (new_title2, seq2, qual2))
                    old_handle1.close()
                    old_handle2.close()
                    new_handle1.close()
                    new_handle2.close()
                    shutil.move(out_file1 + ".tmp", out_file1)
                    shutil.move(out_file2 + ".tmp", out_file2)

        rename_masked_filenames(out_format, bool(fastq2))

    elif mode == "count":
        bc_grouping = match_and_count(dict(bcodes), bc_matched, mismatch)

    if not out_file:
        out_file = fastq1.split(".txt")[0] + "_barcodes.yaml"

    with open(out_file, "w+") as out_handle:
        yaml.dump(bc_grouping.__dict__, out_handle, width=70)


def trim_sequence(sequence, offset, length):
    """Returns the trimmed sequence as well as the barcode.
    """
    if offset == 0:
        minus_offset = None
    else:
        minus_offset = -offset

    barcode = sequence[-(offset + length):minus_offset].strip()
    trimmed_sequence = sequence[:-(offset + length)]

    return trimmed_sequence, barcode


def put_back_unmatched_barcodes_in_sequences(out_format):
    """One might want to run the unmatched sequences through some other
    demultiplexing script, then these must be in the correct format.
    """
    unmatched_filename = out_format.replace("--b--", "unmatched").replace("--r--", "1")
    new_unmatched_filename = out_format.replace("--b--", "formattedunmatched").replace("--r--", "1")
    unmatched_fastq = open(unmatched_filename)
    new_unmatched_fastq = open(new_unmatched_filename, "w")
    for title, sequence, quality in FastqGeneralIterator(unmatched_fastq):
        new_unmatched_fastq.write("@%s\n%s\n+\n%s\n" % (title[:-6], sequence + title[-6:] + "A", quality))

    unmatched_fastq.close()
    new_unmatched_fastq.close()


def convert_old_fastq_header(title, barcode, old_fastq_header_re, pair="1"):
    """Convert a header/title from a fastq file from the old format to the
    new format (after CASAVA 1.8).
    """
    m = old_fastq_header_re.match(title.rstrip())
    if m:
        matchdict = m.groupdict()
        matchdict[pair] = pair
        new_title = \
        "%(instrument)s:::%(fc_lane)s:%(tile)s:%(x)s:%(y)s %(pair)s:::" \
        % matchdict
        new_title += barcode
    else:
        new_title = title

    return new_title


def _get_barcode_statistics(exp_barcodes, obs_barcodes):
    """Get the prevalence of a particular barcode, as a number and a frequency
    """
    total = float(sum(obs_barcodes.itervalues()))
    print("Total after read:\t%.0f" % (total,))
    bcm_nums = []
    bcm_parts = []
    for bc in exp_barcodes:
        num = obs_barcodes.get(bc, 0)
        part = float(num) / total
        bcm_nums.append(num)
        bcm_parts.append(part)
    return (bcm_nums, bcm_parts)


def _get_common_barcodes(bcodes, cutoff):
    """Get the barcodes that occur at a frequency above cutoff
    """
    # Seperate out the most common barcodes
    total = float(sum(bcodes.itervalues()))
    bc_matched = []
    for bc, num in bcodes.iteritems():
        part = float(num) / total
        if part >= cutoff:
            bc_matched.append(bc)

    return bc_matched


def _get_run_info_barcodes(run_info_file, lane=0, length=6):
    """Extract the barcodes to demultiplex against from run_info file"""

    barcodes = set([])
    with open(run_info_file) as fh:
        run_info = yaml.load(fh)
        for lane_info in run_info:
            if (lane != 0 and int(lane_info.get("lane", 0)) != lane):
                continue
            for bc in lane_info.get("multiplex", {}):
                barcodes.add(bc["sequence"][:length])

    return list(barcodes)


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


def _match_barcodes(bcode, given_bcodes, mismatch, masked=False):
    """Logic for matching a barcode against the given barcodes"""

    # First check for perfect matches
    matched = ""
    if bcode in given_bcodes:
        matched = bcode
    # If a perfect match could not be found, do a finer matching but only if
    # we allow mismatches or the given barcodes contain masked positions.
    elif mismatch > 0 or masked:
        for gbc in given_bcodes:
            current_mismatch = bc_mismatched(bcode, gbc, mismatch)
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
    # Set a flag indicating whether the given
    # barcodes contain masked nucleotides
    masked = False
    for gbc in given_bcodes:
        if 'N' in gbc:
            masked = True
            break

    for bc, count in bcodes.iteritems():
        match = _match_barcodes(bc, given_bcodes, mismatch, masked)
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


def match_and_merge(bcodes, given_bcodes, mismatch, format, paired=False):
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

    # Set a flag indicating whether the given
    # barcodes contain masked nucleotides
    masked = False
    for gbc in given_bcodes:
        if 'N' in gbc:
            masked = True
            break

    for bc, count in bcodes.iteritems():
        match = _match_barcodes(bc, given_bcodes, mismatch, masked)
        if len(match) > 0:
            if match not in bc_grouping.matched:
                bc_grouping.matched[match] = {"variants": [], "count": 0}
            if bc not in bc_grouping.matched[match]["variants"]:
                bc_grouping.matched[match]["variants"].append(bc)

            bc_grouping.matched[match]["count"] += count
            found_bcodes.add(bc)

            if match != bc:
                merge_matched_files(match, bc, format, merger, paired)
        else:
            merge_matched_files("unmatched", bc, format, merger, paired)

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


def merge_matched_files(primary_bc, matched_bc, format, merger, paired=False):
    """Merges the fastq file corresponding to matched_bc in to the fastq file
    corresponding to primary_bc.
    """
    matched_file = format.replace("--r--", "1").replace("--b--", matched_bc)
    primary_file = format.replace("--r--", "1").replace("--b--", primary_bc)

    merger(matched_file, primary_file)

    os.remove(matched_file)

    if paired:
        matched_file = format.replace("--r--", "2").replace("--b--", matched_bc)
        primary_file = format.replace("--r--", "2").replace("--b--", primary_bc)

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
            opened = False
            while not opened:
                try:
                    out_handle = open(target_file, "a")
                    self[target_file] = out_handle
                    opened = True
                except IOError as e:
                    if e.errno == 24:
                        print("Closing 100 files")
                        for (name, handle) in self.items()[:100]:
                            handle.close()
                            del self[name]
                        pass
                    else:
                        raise e

        with open(source_file, "r") as in_handle:
            shutil.copyfileobj(in_handle, out_handle)

    def close(self):
        """Close all the opened files.
        """
        for handle in self.itervalues():
            handle.close()


class FileWriter(dict):
    """A class for efficientely writing to several file
    """
    def __init__(self, output_base):
        self.output_base = output_base
        self.work_dir = os.path.dirname(output_base)
        if not os.path.exists(self.work_dir) and self.work_dir:
            try:
                os.makedirs(self.work_dir)
            except OSError:
                assert os.path.isdir(self.work_dir)

    def __call__(self, barcode, name1, seq1, qual1, name2, seq2, qual2):
        self.write_reads(barcode, name1, seq1, qual1, name2, seq2, qual2)

    def write_reads(self, barcode, name1, seq1, qual1, name2, seq2, qual2):
        read1name = self.output_base.replace("--r--", "1")
        read1name = read1name.replace("--b--", barcode)
        self.write_to_handles(name1, seq1, qual1, read1name)
        if name2:
            read2name = self.output_base.replace("--r--", "2")
            read2name = read2name.replace("--b--", barcode)
            self.write_to_handles(name2, seq2, qual2, read2name)

    def write_to_handles(self, name, seq, qual, fname):
        try:
            out_handle = self[fname]
        except KeyError:
            opened = False
            while not opened:
                try:
                    out_handle = open(fname, "a")
                    self[fname] = out_handle
                    opened = True
                except IOError as e:
                    if e.errno == 24:
                        print("Closing some of the %s files" % (len(self),))
                        for (name, handle) in self.items()[:100]:
                            handle.close()
                            del self[name]
                        pass
                    else:
                        raise e

        out_handle.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))

    def close(self):
        """Close all the opened files.
        """
        for handle in self.itervalues():
            handle.close()


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


def rename_masked_filenames(out_format, paired=False):
    """Output files with N in the filename should be renamed to match an
    unmasked barcode.
    """
    out_glob = out_format.replace("--b--", "*")
    out_glob = out_glob.replace("--r--", "1")
    created_files = glob.glob(out_glob)
    out_glob_re = out_format.replace("--b--", r"(?P<multiplex_id>[\d\w]*)")
    out_glob_re = out_glob_re.replace("--r--", r"(?P<pair>\d)")
    glob_re = re.compile(out_glob_re)
    for filename in created_files:
        match_dict = glob_re.match(filename).groupdict()
        barcode = match_dict["multiplex_id"]
        if "N" not in barcode:
            continue

        for title, _, _ in FastqGeneralIterator(open(filename)):
            other_barcode = title[-6:]
            if "N" not in other_barcode:
                break

        os.rename(filename, out_format.replace("--b--", other_barcode).replace("--r--", "1"))
        print("Changed %s -> %s" % (filename, other_barcode))
        if paired:
            filename2 = out_format.replace("--b--", barcode).replace("--r--", "2")
            new_filename2 = out_format.replace("--b--", other_barcode).replace("--r--", "2")
            os.rename(filename2, new_filename2)

if __name__ == "__main__":
    usage = \
    """Usage: %prog <fastq1> [<fastq2>] [<run info yaml file>] [options]

    If called only with the fastq file, the bar codes will be matched to and
    grouped with bar codes from the Illumina documentation.
    If also the <run info yaml file> is given, the bar codes in the fastq file
    will be matched to the bar codes in the run info file. Index names will
    still be extracted from the Illumina documentation table.
    """
    parser = OptionParser(usage=usage)
    parser.add_option("--lane", dest="lane", default=0, \
    help="Specifies a lane (from the run_info.yaml) to use demultiplex \
    information from. Setting this to 0 means 'look at all lanes', this is \
    also the default.")
    parser.add_option("-o", "--out_file", dest="out_file", default=None, \
     help="The name of the file where the statistics about the barcodes will \
     be stored. If not given, the statistics will be saved in a file with \
     the same name as the fastq file but with '.txt' replaced by \
     '_barcodes.yaml'")
    parser.add_option("-l", "--length", dest="length", default=6, \
    help="The length of the barcode sequences, default is 6.")
    parser.add_option("-b", "--back", dest="offset", default=0, \
    help="The number of characters there are after the barcode in the \
    sequence.\n For example, sometimes Illumina barcodes have an 'A' after \
    the barcode, if this is not taken in consideration for the barcodes used \
    to match against there will be a mismatch, and in that case one would \
    supply '-b 1' as an option.")
    parser.add_option("-m", "--mismatch", dest="mismatch", default=1, \
    help="The number of mismatches to allow. Default is 1.")
    parser.add_option("-v", "--verbose", dest="verbose", default=False, \
    action="store_true", help="Print the statistics (which are stored in the \
    out_file) to stdout after the script is complete.")
    parser.add_option("-c", "--cutoff", dest="cutoff", default=0.02, \
    help="The frequency of a barcode which is needed for the barcode to be \
    considered matched, and used to match other barcodes against. This is only \
    used when not supplying a run_info or barcode-file.\n Default is 0.02 (2%).")
    parser.add_option("-n", "--dryrun", dest="dry_run", default=False, \
    action="store_true", help="Not implemented.")
    parser.add_option("--mode", dest="mode", default="demultiplex", \
    help="The mode for the script to run in. Default is 'demultiplex'.\n \
    There are two available modes: 'demultiplex', where files for each barcode \
    are saved. As well as the mode 'count', where no actual demultiplexing \
    takes place, but statistics and matching information is gathered and saved \
    in the out_file.")
    parser.add_option("--DEBUG", dest="debug", default=False, \
    action="store_true", help="Start pdb when an error occurs. \
    False if not supplied.")
    parser.add_option("--out_format", dest="out_format", \
    default="demultiplex_out/out_--b--_--r--_fastq.txt", \
    help="The format of the demultiplexed files which will be saved during \
    the run of the script. This should contain '--b--', where the barcode \
    sequence will be substituted in, and '--r--', where the pair number will \
    be substituted (will be '1' when run on fastq which isn't paired).")
    parser.add_option("-p", "--paired", dest="paired", default=False, \
    action="store_true", help="Indicate when supplying two arguments \
    whether they refer to a pair of fastq files or a fastq file along with \
    file with barcode information. If not given, two files are considered to \
    be a single fastq file and a barcode information file.")
    options, args = parser.parse_args()
    if len(args) == 1:
        fastq1, = args
        fastq2 = None
        run_info = None
    elif len(args) == 2 and options.paired:
        fastq1, fastq2 = args
        run_info = None
    elif len(args) == 2 and not options.paired:
        fastq1, run_info = args
        fastq2 = None
    elif len(args) == 3:
        fastq1, fastq2, run_info = args
    else:
        parser.print_help()
        sys.exit()

    # import cProfile
    # cProfile.run("main(fastq, run_info, int(options.lane), options.out_file,\
    # int(options.length), int(options.offset), int(options.mismatch), \
    # options.verbose, float(options.cutoff), options.dry_run, options.mode)",\
    # sort='cumulative')
    try:
        main(fastq1, fastq2, run_info, int(options.lane), options.out_file, \
        int(options.length), int(options.offset), int(options.mismatch), \
        options.verbose, float(options.cutoff), options.dry_run, \
        options.mode, \
        options.out_format)
    except Exception as e:
        if options.debug:
            import pdb
            import traceback
            type, value, tb = sys.exc_info()
            traceback.print_exc()
            pdb.post_mortem(tb)
        else:
            raise e


# --- TESTS

class BarcodeTest(unittest.TestCase):
    """Test the methods in the demultiplexing script.
    """
    def setUp(self):
        self.example_sequence = \
        "AGAGGAAAAGAGGAAGAGAGGAATATAAAATTTATTTTTGCTACTATTTA" + \
        "TTTTACTATTTGACCATTTTTACTGCTATTTTCACTCACTCAAATATGGTTGCCAATA"

    def test_trim_sequence(self):
        """Test trimming the barcode from the sequence.
        """
        desired_sequence = \
        "AGAGGAAAAGAGGAAGAGAGGAATATAAAATTTATTTTTGCTACTATTTA" + \
        "TTTTACTATTTGACCATTTTTACTGCTATTTTCACTCACTCAAATATGGTT"
        desired_barcode = "GCCAAT"
        trimmed_sequence, barcode = trim_sequence(self.example_sequence, 1, 6)
        assert desired_sequence == trimmed_sequence, "Sequence not trimmed correctly"
        assert desired_barcode == barcode, "Barcode not extracted correctly"

    def test_convert_old_fastq_header(self):
        """Test converting a fastq header of the old format in
        to the new format.
        """
        old_fastq_header_re = re.compile(r"(?P<instrument>[\d\w-]*)"
        ":(?P<fc_lane>\d*):(?P<tile>\d*):(?P<x>\d*)"
        ":(?P<y>\d*)#(?P<multiplex_id>[\d\w]*)/(?P<pair>\d)")
        old_header = "@HWI-ST1018:8:1101:13101:38091#0/1"
        desired_header = "HWI-ST1018:::8:1101:13101:38091 1:::AAAAAA"
        new_header = convert_old_fastq_header(old_header[1:], "AAAAAA", old_fastq_header_re)
        assert new_header == desired_header, "Header conversion failed: \n%s\n%s" % \
        (new_header, desired_header)
