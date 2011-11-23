"""Pipeline support for barcode analysis and de-mulitplexing.
"""
import os
import copy
import subprocess

from Bio import SeqIO

from bcbio import utils
from bcbio.pipeline.fastq import get_fastq_files

def split_by_barcode(fastq1, fastq2, multiplex, base_name, dirs, config):
    """Split a fastq file into multiplex pieces using barcode details.
    """
    if not multiplex:
        return [("", "", fastq1, fastq2)]
    bc_dir = os.path.join(dirs["work"], "%s_barcode" % base_name)
    nomatch_file = "%s_unmatched_1_fastq.txt" % base_name
    metrics_file = "%s_bc.metrics" % base_name
    out_files = []
    for info in multiplex:
        fq_fname = lambda x: os.path.join(bc_dir, "%s_%s_%s_fastq.txt" %
                            (base_name, info["barcode_id"], x))
        bc_file1 = fq_fname("1")
        bc_file2 = fq_fname("2") if fastq2 else None
        out_files.append((info["barcode_id"], info["name"], bc_file1, bc_file2))
    with utils.chdir(bc_dir):
        if not os.path.exists(nomatch_file) and not os.path.exists(metrics_file):
            tag_file = _make_tag_file(multiplex,config)
            cl = [config["program"]["barcode"], tag_file,
                  "%s_--b--_--r--_fastq.txt" % base_name,
                  fastq1]
            if fastq2:
                cl.append(fastq2)
            cl.append("--mismatch=%s" % config["algorithm"]["bc_mismatch"])
            cl.append("--metrics=%s" % metrics_file)
            if int(config["algorithm"]["bc_read"]) == 2:
                cl.append("--second")
            if int(config["algorithm"]["bc_position"]) == 5:
                cl.append("--five")
            if config["algorithm"].get("bc_allow_indels", True) is False:
                cl.append("--noindel")
            if "bc_offset" in config["algorithm"]:
                cl.append("--bc_offset=%s" % config["algorithm"]["bc_offset"])
            with utils.file_transaction(out_files + [nomatch_file, metrics_file]):
                cl = [os.path.expandvars(command) for command in cl]
                subprocess.check_call(cl)

    # Prunes files that are not present in the filesystem
    out_files = [(b, n, f1, f2) for (b, n, f1, f2) in out_files if os.path.exists(f1)]
    return out_files

def _make_tag_file(barcodes,config):
    tag_file = "%s-barcodes.cfg" % barcodes[0].get("barcode_type", "barcode")
    barcodes = _adjust_illumina_tags(barcodes,config)
    with open(tag_file, "w") as out_handle:
        for bc in barcodes:
            out_handle.write("%s %s\n" % (bc["barcode_id"], bc["sequence"]))
    return tag_file

def _adjust_illumina_tags(barcodes,config):
    """Handle additional trailing A in Illumina barocdes.

    Illumina barcodes are listed as 6bp sequences but have an additional
    A base when coming off on the sequencer. This checks for this case and
    adjusts the sequences appropriately if needed. In case the 
    bc_illumina_no_trailing configuration option is set to true, this method
    will instead make sure that the barcodes do not include the trailing A
    and that demultiplexing will not attempt to match it.
    """
    
    skip_a = config["algorithm"].get("bc_illumina_no_trailing",False)
    # Set the bc_offset parameter if we're skipping the trailing A
    if skip_a:
        config["algorithm"]["bc_offset"] = 1
    illumina_size = 7
    all_illumina = True
    need_a = False
    for bc in barcodes:
        # Will only process in case all barcodes are illumina
        if bc.get("barcode_type", "illumina").lower().find("illumina") == -1:
            return barcodes
    new_bc = copy.deepcopy(barcodes)
    for bc in new_bc:
        if (not skip_a and (
            not bc["sequence"].upper().endswith("A") or
            len(bc["sequence"]) < illumina_size)):
            bc["sequence"] = "%sA" % bc["sequence"]
        elif (skip_a and bc["sequence"].upper().endswith("A") and
              len(bc["sequence"]) == illumina_size):
            bc["sequence"] = bc["sequence"][:illumina_size-1]
    barcodes = new_bc
    return barcodes

def add_multiplex_across_lanes(run_items, fastq_dir, fc_name):
    """Add multiplex information to control and non-multiplexed lanes.

    Illumina runs include barcode reads for non-multiplex lanes, and the
    control, when run on a multiplexed flow cell. This checks for this
    situation and adds details to trim off the extra bases.
    """
    fastq_dir = utils.add_full_path(fastq_dir)
    # determine if we have multiplexes and collect expected size
    fastq_sizes = []
    tag_sizes = []
    has_barcodes = False
    for item in run_items:
        if item.get("multiplex", None):
            has_barcodes = True
            tag_sizes.extend([len(b["sequence"]) for b in item["multiplex"]])
            fastq_sizes.append(_get_fastq_size(item, fastq_dir, fc_name))
    if not has_barcodes:  # nothing to worry about
        return run_items
    fastq_sizes = list(set(fastq_sizes))

    # discard 0 sizes to handle the case where lane(s) are empty or failed
    try:
        fastq_sizes.remove(0)
    except ValueError: pass

    tag_sizes = list(set(tag_sizes))
    final_items = []
    for item in run_items:
        if item.get("multiplex", None) is None:
            assert len(fastq_sizes) == 1, \
                   "Multi and non-multiplex reads with multiple sizes"
            expected_size = fastq_sizes[0]
            assert len(tag_sizes) == 1, \
                   "Expect identical tag size for a flowcell"
            tag_size = tag_sizes[0]
            this_size = _get_fastq_size(item, fastq_dir, fc_name)
            if this_size == expected_size:
                item["multiplex"] = [{"name" : item.get("name", item["description"]),
                                      "barcode_id": "trim",
                                      "sequence" : "N" * tag_size}]
            else:
                assert this_size == expected_size - tag_size, \
                       "Unexpected non-multiplex sequence"
        final_items.append(item)
    return final_items

def _get_fastq_size(item, fastq_dir, fc_name):
    """Retrieve the size of reads from the first flowcell sequence.
    """
    (fastq1, _) = get_fastq_files(fastq_dir, item, fc_name)
    with open(fastq1) as in_handle:
        try:
            rec = SeqIO.parse(in_handle, "fastq").next()
            size = len(rec.seq)
        except StopIteration:
            size = 0
    return size

