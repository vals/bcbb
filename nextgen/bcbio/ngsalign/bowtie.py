"""Next gen sequence alignments with Bowtie (http://bowtie-bio.sourceforge.net).
"""
import os
import subprocess

from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction

galaxy_location_file = "bowtie_indices.loc"

def _bowtie_args_from_config(config):
    """Configurable high level options for bowtie.
    """
    qual_format = config["algorithm"].get("quality_format", None)
    if qual_format is None or qual_format.lower() == "illumina":
        qual_flags = ["--phred64-quals"]
    else:
        qual_flags = []
    multi_mappers = config["algorithm"].get("multiple_mappers", True)
    multi_flags = ["-M", 1] if multi_mappers else ["-m", 1]
    cores = config.get("resources", {}).get("bowtie", {}).get("cores", None)
    core_flags = ["-p", str(cores)] if cores else []
    return core_flags + qual_flags + multi_flags

def align(fastq_file, pair_file, ref_file, out_base, align_dir, config,
          extra_args=None, rg_name=None):
    """Do standard or paired end alignment with bowtie.
    """
    out_file = os.path.join(align_dir, "%s.sam" % out_base)
    if not file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            cl = [config["program"]["bowtie"]]
            cl += _bowtie_args_from_config(config)
            cl += extra_args if extra_args is not None else []
            cl += ["-q",
                   "-v", config["algorithm"]["max_errors"],
                   "-k", 1,
                   "-X", 2000, # default is too selective for most data
                   "--best",
                   "--strata",
                   "--sam",
                   ref_file]
            if pair_file:
                cl += ["-1", fastq_file, "-2", pair_file]
            else:
                cl += [fastq_file]
            cl += [tx_out_file]
            cl = [str(i) for i in cl]
            subprocess.check_call(cl)
    return out_file


def remove_contaminants(fastq_file, pair_file, ref_file, out_base, fastq_dir, config,
                        extra_args=None, rg_name=None):
    """Remove reads aligning to the contaminating reference genome 
    """
    
    tmp_file = os.path.join(fastq_dir, "%s_clean" % out_base)
    out_file = []
    if not len(glob.glob("%s*" % tmp_file)) > 0:
        with file_transaction(tmp_file) as tx_out_file:
            cl = [config["program"]["bowtie"]]
            cl += _bowtie_args_from_config(config)
            cl += extra_args if extra_args is not None else []
            cl += ["-un", tx_out_file,
                   ref_file]
            if pair_file:
                cl += ["-1", fastq_file, "-2", pair_file]
            else:
                cl += [fastq_file]
            cl += ["/dev/null"]
            cl = [str(i) for i in cl]
            subprocess.check_call(cl)
            
            # Rename the temporary files
            if pair_file:
                for i in ("1","2"):
                    tfile = "%s_%s" % (tmp_file,i)
                    out_file.append("%s_fastq.txt" % tfile)
                    os.rename(tfile,out_file[-1])
            else:
                out_file = ["%s_fastq.txt" % tmp_file,None]
                os.rename(tmp_file,out_file[0])
            
    return out_file
