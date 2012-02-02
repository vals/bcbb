"""A wrapper for providing a command line argument interface to
count_barcode_demultiplex.py which is compatible with the demultiplex module
of BCBB.
"""
import subprocess
from optparse import OptionParser


def main(fastq1, fastq2, barcode_info_file, mismatch, outfile, back, length):
    cl = ["count_barcode_demultiplex.py"]
    cl.append(fastq1)
    if fastq2:
        cl.append(fastq2)
    cl.append(barcode_info_file)
    cl.append("--mismatch=%s" % mismatch)
    cl.append("--out_file=%s" % outfile)
    cl.appenf("--length=%s" % length)
    cl.append("--back=%s" % back)
    if fastq2:
        cl.append("--paired")

    subprocess.check_call(cl)

if __name__ == "__main__":
    usage = "Wrapper for count_barcode_demultiplex.py"
    parser = OptionParser(usage=usage)

    parser.add_option("--mismatch", dest="mismatch", default=1)
    parser.add_option("--metrics", dest="metrics_file", default=None)
    parser.add_option("--second", dest="first_read", action="store_false", \
    default=True)
    parser.add_option("--five", dest="three_end", action="store_false", \
    default=True)
    parser.add_option("--noindel", dest="indels", action="store_false", \
    default=True)
    parser.add_option("--bc_offset", dest="bc_offset", default=0)

    options, args = parser.parse_args()

    if len(args) == 3:
        barcode_info_file, out_format, fastq1 = args
        fastq2 = None
    elif len(args) == 4:
        barcode_info_file, out_format, fastq1, fastq2 = args

    length = 6

    if options.three_end:
        back = int(options.bc_offset)
    else:
        back = -length - int(options.bc_offset)

    mismatch = int(options.mismatch)
    outfile = options.metrics_file

    main(fastq1, fastq2, barcode_info_file, mismatch, outfile, back, length)
