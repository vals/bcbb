"""Pipeline utilities to retrieve FASTQ formatted files for processing.
"""
import os
import glob
import subprocess
import contextlib
import collections

import pysam

from bcbio import broad
from bcbio.utils import file_exists, safe_makedir
from bcbio.distributed.transaction import file_transaction


def get_fastq_files(directory, work_dir, item, fc_name, bc_name=None, glob_ext="_fastq.txt",
                    config=None, unpack=True):
    """Retrieve fastq files for the given lane, ready to process.
    """
    if "files" in item and bc_name is None:
        names = item["files"]
        if isinstance(names, basestring):
            names = [names]
        files = [x if os.path.isabs(x) else os.path.join(directory, x) for x in names]

    else:
        assert fc_name is not None
        lane = item["lane"]
        if bc_name:
            glob_str = "%s_*%s_%s_*%s" % (lane, fc_name, bc_name, glob_ext)
        else:
            glob_str = "%s_*%s*%s" % (lane, fc_name, glob_ext)
        files = glob.glob(os.path.join(directory, glob_str))
        
        # Include gzipped files
        glob_str = "%s.gz" % glob_str
        files.extend(glob.glob(os.path.join(directory, glob_str)))
        
        files.sort()
        if len(files) > 2 or len(files) == 0:
            raise ValueError("Did not find correct files for %s %s %s %s" %
                    (directory, lane, fc_name, files))
    ready_files = []
    for fname in files:
        if fname.endswith(".gz") and unpack:
            # TODO: Parallelize using pgzip
            ready_name = os.path.splitext(fname)[0]
            ready_files.append(ready_name)
            if not os.path.exists(ready_name):
                cl = ["gunzip", fname]
                subprocess.check_call(cl)
        elif fname.endswith(".bam"):
            ready_files = convert_bam_to_fastq(fname, work_dir, config)
        else:
            assert os.path.exists(fname), fname
            ready_files.append(fname)
    ready_files = [x for x in ready_files if x is not None]
    return ready_files[0], (ready_files[1] if len(ready_files) > 1 else None)


def convert_bam_to_fastq(in_file, work_dir, config):
    """Convert BAM input file into FASTQ files.
    """
    out_dir = safe_makedir(os.path.join(work_dir, "fastq_convert"))
    out_files = [os.path.join(out_dir, "{0}_{1}.fastq".format(
                 os.path.splitext(os.path.basename(in_file))[0], x))
                 for x in ["1", "2"]]
    if _is_paired(in_file):
        out1, out2 = out_files
    else:
        out1 = out_files[0]
        out2 = None
    if not file_exists(out1):
        broad_runner = broad.runner_from_config(config)
        broad_runner.run_fn("picard_bam_to_fastq", in_file, out1, out2)
    if os.path.getsize(out2) == 0:
        out2 = None
    return [out1, out2]

def _is_paired(bam_file):
    # XXX need development version of pysam for this to work on
    # fastq files without headers (ie. FastqToSam)
    # Instead return true by default and then check after output
    return True
    with contextlib.closing(pysam.Samfile(bam_file, "rb")) as work_bam:
        for read in bam_file:
            return read.is_paired


def get_single_fastq_files(lane, fc_dir, fc_name):
    return get_fastq_files(fc_dir, lane, fc_name)

# TODO: these two could probably also be handled much more efficiently
def get_barcoded_project_files(multiplex, lane, fc_dir, fc_name):
    fq = list()
    for bc in multiplex:
        fq.append(get_fastq_files(fc_dir, lane, fc_name, bc['name'], ".fastq"))
    return fq

def get_barcoded_fastq_files(multiplex, item, fc_dir, fc_name, fc_date):
    fq = list()
    bc_dir = "%s_%s_%s_barcode" % (item["lane"], fc_date, fc_name)
    bc_dir = os.path.join(fc_dir, bc_dir)
    if multiplex is None:
        fq.append(get_fastq_files(bc_dir, item, fc_name))
    else:
        for bc in multiplex:
            if not os.path.exists(bc_dir):
                raise IOError("No barcode directory found: " + str(bc_dir))
            fq.append(get_fastq_files(bc_dir, item, fc_name, bc_name=bc['barcode_id']))
    return fq


# TODO: these two could probably be handled much more efficiently
def convert_barcode_id_to_name(multiplex, fc_name, fq):
    """Convert barcode id to sample description, changing extension from _fastq.txt to .fastq in the process"""
    fqout = list([None, None])
    if multiplex is None:
        fqout[0] = fq[0]
        if not fq[1] == None:
            fqout[1] = fq[1]
    else:
        bcid2name = dict([(mp['barcode_id'], mp['name']) for mp in multiplex])
        for bcid in bcid2name.keys():
            mstr = "%s_%s_" % (fc_name, bcid) 
            if fq[0].find(mstr) != -1:
                from_str = "%s_%s_" %(fc_name, bcid)
                to_str   = "%s_%s_" %(fc_name, bcid2name[bcid])
                fqout[0] = fq[0].replace(from_str, to_str)
                if not fq[1] == None:
                    fqout[1] = fq[1].replace(from_str, to_str)
    fqout[0] = fqout[0].replace("_fastq.txt", ".fastq")
    if not fqout[1] == None:
        fqout[1] = fqout[1].replace("_fastq.txt", ".fastq")
    return os.path.basename(fqout[0]), (os.path.basename(fqout[1]) if len(fqout) > 1 else None)

def convert_name_to_barcode_id(multiplex, fc_name, fq):
    """Convert sample description to barcode id, changing extension from .fastq to _fastq.txt in the process"""
    fqout = list([None, None])
    name2bcid = dict([(mp['name'], mp['barcode_id']) for mp in multiplex])
    for name in name2bcid.keys():
        mstr = "%s_%s_" % (fc_name, name) 
        if fq[0].find(mstr) != -1:
            from_str = "%s_%s_" %(fc_name, name)
            to_str   = "%s_%s_" %(fc_name, name2bcid[name])
            fqout[0] = fq[0].replace(from_str, to_str)
            if not fq[1] == None:
                fqout[1] = fq[1].replace(from_str, to_str)
    fqout[0] = fqout[0].replace(".fastq", "_fastq.txt")
    if not fqout[1] == None:
        fqout[1] = fqout[1].replace(".fastq", "_fastq.txt")
    return os.path.basename(fqout[0]), (os.path.basename(fqout[1]) if len(fqout) > 1 else None)

def get_multiplex_items(multiplex, lane, fc_dir, fc_name, fc_date):
    mitems = list()
    bc_dir = "%s_%s_%s_barcode" % (lane, fc_date, fc_name)
    bc_dir = os.path.join(fc_dir, bc_dir)
    lane_name = "%s_%s_%s" % (lane, fc_date, fc_name)
    for bc in multiplex:
        mname = bc['barcode_id']
        mlane_name = "%s_%s" % (lane_name, mname) if mname else lane_name
        msample = bc['name']
        if msample is None:
            msample = "%s---%s" % (sample_name, mname)
        if not os.path.exists(bc_dir):
            raise IOError("No barcode directory found: " + str(bc_dir))
        fastq1, fastq2 = get_fastq_files(bc_dir, {"lane" : lane}, fc_name, bc_name=bc['barcode_id'])
        mitems.append((fastq1, fastq2 , mlane_name, msample))
    return mitems
