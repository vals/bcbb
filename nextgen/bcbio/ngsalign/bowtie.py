"""Next gen sequence alignments with Bowtie (http://bowtie-bio.sourceforge.net).
"""
import os
import subprocess
import glob

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
    
    out_root = os.path.join(fastq_dir,out_base)
    out_files = ["%s_1.ext" % out_root,
                 "%s_2.ext" % out_root,
                 "%s.filter_metrics" % out_root]
    suffix = "_fastq.txt"
    
    if not len(glob.glob("%s_[12]%s" % (out_root,suffix))) > 0:
        with file_transaction(out_files) as (tx_out_file1, tx_out_file2, tx_metrics_file):
            out = tx_out_file1
            if pair_file:
                out = out.replace("_1.ext",".ext")
            
            cl = [config["program"]["bowtie"]]
            cl += _bowtie_args_from_config(config)
            cl += extra_args if extra_args is not None else []
            # Allow for read pairs mapping at opposite ends of e.g. the phiX genome
            # Trim 7bp from 3'end corresponding to the barcode
            cl += ["--best", "-X", "6000","-3","7"]
            cl += ["--un", out,
                   ref_file]
            if pair_file:
                cl += ["-1", fastq_file, "-2", pair_file]
            else:
                cl += [fastq_file]
            cl += ["/dev/null"]
            cl = [str(i) for i in cl]
            
            # Get the output, echo it as well as write it to the metrics file
            output = subprocess.check_output(cl,stderr=subprocess.STDOUT)
            print output
            
            with open(tx_metrics_file,"w") as fh:
                fh.write("%s\n" % str(output))
            
        dest_files = []
        for i, out_file in enumerate(out_files):
            if not out_file.endswith(".ext"): continue
            if not os.path.exists(out_file):
                if i == 1 and not pair_file:
                    dest_files.append(pair_file)
                    continue 
                open(out_file,"w").close()
            dest_file = out_file.replace(".ext",suffix)
            os.rename(out_file,dest_file)
            dest_files.append(dest_file)
    else:
        dest_files = glob.glob("%s_[12]%s" % (out_root,suffix))    
    
    dest_files.append(out_base)
    return dest_files
