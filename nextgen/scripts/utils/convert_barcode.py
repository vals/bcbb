#!/usr/bin/env python
"""
Rename barcodes in file names.

Usage:
  convert_barcode.py <YAML run information> <glob>
                     --dry_run --verbose --recursive
                     --name_to_barcode
"""

import os
import sys
import re

from optparse import OptionParser
from bcbio.pipeline.flowcell import Flowcell, Lane

import yaml
import fnmatch
import shutil

def main(run_info_yaml, glob_str):
    fp = open(run_info_yaml)
    run_info_structure = yaml.load(fp)
    fp.close()
    fc = Flowcell(None, None, run_info_structure)
    if options.name_to_barcode:
        bcmap = _name_to_barcode_id(fc)
    else:
        bcmap = _barcode_id_to_name(fc)
    
    if options.recursive:
        for root, dirnames, filenames in os.walk("./"):
            for filename in fnmatch.filter(filenames, glob_str):
                rename_file(filename, bcmap)
    else:
        filenames = os.listdir("./")
        for filename in fnmatch.filter(filenames, glob_str):
            rename_file(filename, bcmap)

def rename_file(filename, bcmap):
    if options.verbose:
        print filename
    lane, date, fc, bc = _get_flowcell_info(filename)
    from_str = "%s_%s_%s_%s" % (lane, date, fc, bc)
    to_str  = "%s_%s_%s_%s" % (lane, date, fc, bcmap[lane][str(bc)])
    print from_str
    print to_str
    src = filename
    tgt = src.replace(from_str, to_str)
    if options.dry_run:
        print "DRY_RUN: renaming %s to %s" % (src, tgt)
    else:
        print "renaming %s to %s" % (src, tgt)
        shutil.move(src, tgt)
        
def _barcode_id_to_name(fc):
    bcid2name = dict()
    for lane in fc.get_lanes():
        multiplex = lane.get_samples()
        bcid2name[lane.get_name()] = dict([(str(mp.get_barcode_id()), mp.get_barcode_name()) for mp in multiplex])
    return bcid2name

def _name_to_barcode_id(fc):
    name2bcid = dict()
    for lane in run_info:
        multiplex = lane.get_samples()
        bcid2name[lane] = dict([(mp.get_barcode_name(), str(mp.get_barcode_id())) for mp in multiplex])
    return name2bcid

def _get_flowcell_info(filename):
    fn = os.path.basename(filename)
    regexp = r'^([0-9])_([0-9]{6})_([A-Za-z0-9]+)_([\_A-Za-z0-9]+)\.'
    m = re.search(regexp, fn, re.I)
    if not m or len(m.groups()) == 0:
        print "regexp %s failed!" % regexp
        sys.exit()
    lane, date, fc, bc = m.group(1), m.group(2), m.group(3), m.group(4)

    assert int(lane) > 0 and int(lane) < 9, "lane is not between 1 and 8"
    return (str(lane), date, fc, str(bc))

if __name__ == "__main__":
    usage = """
    convert_barcode.py <YAML run information> <glob>
                       --dry_run --verbose --recursive
                       --name_to_barcode

    For more extensive help, type convert_barcode.py
    """

    parser = OptionParser(usage=usage)
    parser.add_option("-n", "--dry_run", dest="dry_run", action="store_true",
                      default=False)
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      default=False)
    parser.add_option("-r", "--recursive", dest="recursive", action="store_false",
                      default=True)
    parser.add_option("-b", "--name_to_barcode", dest="name_to_barcode", action="store_true",
                      default=False)
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print __doc__
        sys.exit()
    kwargs = dict(
        )
    main(*args, **kwargs)
