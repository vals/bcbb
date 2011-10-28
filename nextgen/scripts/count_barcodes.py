
# This script uses a few assumptions about how the barcoding works.
# Since the aim is to look at Illumina HiSeq barcodes, the assumption is that
# the barcode is located at the 3' end of read 1 for each paired-end read.

import sys
from operator import itemgetter
from Bio import pairwise2

import yaml

from bcbio.solexa import INDEX_LOOKUP

if len(sys.argv) < 3:
    print "python", sys.argv[0], "<fastq file containing barcode sequence>" + \
    " <position where bar code starts> <offset from end> <run_info yaml>"
    sys.exit(0)

offset = int(sys.argv[3])
bcodes = {}  # collect counts for all observed barcodes

mismatch = 1

cntr = 0
last_was_header = False
pos = int(sys.argv[2])

with open(sys.argv[4]) as in_handle:
    run_info = yaml.load(in_handle)

for line in open(sys.argv[1]):
    if not last_was_header:
        if line[0] == "@":
            last_was_header = True
        continue

    bcode = line[pos:-(offset + 1)].strip()
    if bcode in bcodes:
        bcodes[bcode] += 1
    else:
        bcodes[bcode] = 1
    if last_was_header:
        last_was_header = False

for e in sorted(bcodes.items(), key=itemgetter(1)):
    print e[0] + "\t" + str(e[1])

given_bcodes = [bc["sequence"] for bc in run_info[0]["multiplex"]]

matched_bc = {}
for bc, c in bcodes.items():
    if bc in given_bcodes:
        if bc in matched_bc:
            matched_bc[bc] += c
        else:
            matched_bc[bc] = c

print matched_bc


def approximate_matching(bcodes, given_bcodes, mismatch):
    """ Returns a dectionary with found barcodes along with info
    """
    matched_bc_grouping = {}
    found_bcodes = set()
    if mismatch >= 0:
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
                    print cur_mismatch,
                    matched_bc_grouping[bc_given]["variants"].append(bc)
                    matched_bc_grouping[bc_given]["count"] += count
                    found_bcodes.add(bc)

    print "%s unmatched barcodes" % len(set(bcodes) - found_bcodes)

    return matched_bc_grouping

matched_bc_grouping = approximate_matching(bcodes, given_bcodes, mismatch)

print
print matched_bc_grouping

matched_against_index = approximate_matching(bcodes, list(set(INDEX_LOOKUP.values())), mismatch)
print matched_against_index
