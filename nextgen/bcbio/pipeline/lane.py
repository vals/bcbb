"""Top level driver functionality for processing a sequencing lane.
"""
import os
import copy
import glob
import socket
from difflib import SequenceMatcher

from bcbio.log import logger
from bcbio.pipeline.fastq import get_fastq_files, get_multiplex_items
from bcbio.pipeline.demultiplex import split_by_barcode
from bcbio.pipeline.alignment import align_to_sort_bam, remove_contaminants as rc
from bcbio.solexa.flowcell import get_flowcell_info
from bcbio.bam.trim import brun_trim_fastq

def _get_base_name(file1, file2):

    m = SequenceMatcher(None, file1, file2)
    (a, b, l) = m.find_longest_match(0,len(file1),0,len(file2))
    if a != 0 or b != 0 or l == 0:
        return None
    return file1[0:l].rstrip(" _")

def process_lane(lane_items, fc_name, fc_date, dirs, config):
    """Prepare lanes, potentially splitting based on barcodes.
    """
        
    lane_name = "%s_%s_%s" % (lane_items[0]['lane'], fc_date, fc_name)
    full_fastq1, full_fastq2 = get_fastq_files(dirs["fastq"], dirs["work"],
                                               lane_items[0], fc_name, config=config)
    
    # Filter phiX
    if config["algorithm"].get("filter_phix",False):
        logger.info("Filtering phiX from %s" % lane_name)
        info = {"genomes_filter_out": "phix", "description": lane_items[0].get("description",lane_name)}
        processed = remove_contaminants(full_fastq1, full_fastq2, info, lane_name, info["description"], dirs, config)
        (full_fastq1, full_fastq2, _, lane_name) = processed[0:4]
        
    logger.info("Demultiplexing %s" % base_name)
    bc_files = split_by_barcode(full_fastq1, full_fastq2, lane_items,
                                lane_name, dirs, config)
    out = []
    for item in lane_items:
        config = _update_config_w_custom(config, item)
        # Can specify all barcodes but might not have actual sequences
        # Would be nice to have a good way to check this is okay here.
        if item["barcode_id"] in bc_files:
            fastq1, fastq2 = bc_files[item["barcode_id"]]
            cur_lane_name = lane_name
            cur_lane_desc = item["description"]
            if item.get("name", "") and config["algorithm"].get("include_short_name", True):
                cur_lane_desc = "%s : %s" % (item["name"], cur_lane_desc)
            if item["barcode_id"] is not None:
                cur_lane_name += "_%s" % (item["barcode_id"])
            if config["algorithm"].get("trim_reads", False):
                trim_info = brun_trim_fastq([x for x in [fastq1, fastq2] if x is not None],
                                            dirs, config)
                fastq1 = trim_info[0]
                if fastq2 is not None:
                    fastq2 = trim_info[1]
            out.append((fastq1, fastq2, item, cur_lane_name, cur_lane_desc,
                        dirs, config))
    return out


def remove_contaminants(fastq1, fastq2, info, lane_name, lane_desc,
                      dirs, config):
    """Remove reads mapping to the specified contaminating reference
    """
    
    base_name = None
    genome_build = info.get("genomes_filter_out",None)
    # Skip filtering of phix in case we have already done that for the lane
    if genome_build is not None and not (genome_build == "phix" and config["algorithm"].get("filter_phix",False)) and os.path.exists(fastq1):
        program = config["algorithm"].get("remove_contaminants","bowtie")
        logger.info("Removing %s contaminants on sample %s in lane %s, using %s" % (genome_build,info["description"],lane_name,program))
        fastq1, fastq2 = rc(fastq1,fastq2,genome_build,program,lane_name,dirs,config)
        base_name = _get_base_name(os.path.basename(fastq1), os.path.basename(fastq2))
    return [fastq1, fastq2, info, (base_name or lane_name), lane_desc, dirs, config]

def process_alignment(fastq1, fastq2, info, lane_name, lane_desc,
                      dirs, config):
    """Do an alignment of fastq files, preparing a sorted BAM output file.
    """
    aligner = config["algorithm"].get("aligner", None)
    out_bam = ""
    if os.path.exists(fastq1) and aligner:
        logger.info("Aligning lane %s with %s aligner" % (lane_name, aligner))
        out_bam = align_to_sort_bam(fastq1, fastq2, info["genome_build"], \
                                aligner, lane_name, lane_desc, dirs, config)
    return [{"fastq": [fastq1, fastq2], "out_bam": out_bam, "info": info,
             "config": config}]


def _update_config_w_custom(config, lane_info):
    """Update the configuration for this lane if a custom analysis is specified.
    """
    config = copy.deepcopy(config)
    analysis_type = lane_info.get("analysis", "")
    custom = config["custom_algorithms"].get(analysis_type, None)
    if custom:
        for key, val in custom.iteritems():
            config["algorithm"][key] = val
    # apply any algorithm details specified with the lane
    for key, val in lane_info.get("algorithm", {}).iteritems():
        config["algorithm"][key] = val
    return config


def make_lane_items(info, fc_date, fc_name, dirs, config):
    sample_name = info.get("description", "")
    if (config["algorithm"].get("include_short_name", True) and
            info.get("name", "")):
        sample_name = "%s---%s" % (info.get("name", ""), sample_name)
    genome_build = info.get("genome_build", None)
    multiplex = info.get("multiplex", "")
    logger.info("Processing sample: %s; lane %s; reference genome %s; " \
             "researcher %s; analysis method %s" %
             (sample_name, info["lane"], genome_build,
              info.get("researcher", ""), info.get("analysis", "")))
    lane_items = []
    if multiplex:
        logger.debug("Sample %s is multiplexed as: %s" % (sample_name, multiplex))
        mitems = get_multiplex_items(multiplex, info['lane'], dirs['fc_dir'], fc_name, fc_date)
        for fastq1, fastq2, mlane_name, msample in mitems:
            lane_items.append((fastq1, fastq2, genome_build, mlane_name, msample, dirs, config))
    else:
        # TODO: Not multiplex: what to do?
        pass
    return lane_items


def get_flowcell_id(run_info, fc_dir, check_bc=True, glob_ext="_fastq.txt"):
    for lane in run_info:
        for bc in lane:
            if check_bc:
                glob_str = "%s_*_barcode/*%s" % (bc['lane'], glob_ext)
            else:
                glob_str = "%s_*%s" % (lane, glob_ext)
                next
    files = glob.glob(os.path.join(fc_dir, glob_str))
    try:
        (name, date) = get_flowcell_info(os.path.basename(files[0]))
    except:
        raise StandardError("No flowcell information found in " + str(fc_dir))
    return name, date
