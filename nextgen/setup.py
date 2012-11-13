#!/usr/bin/env python
"""Setup file and install script for NextGen sequencing analysis scripts.
"""
from setuptools import setup, find_packages
import os

setup(name="bcbio-nextgen",
      version="0.3a",
      author="Brad Chapman",
      author_email="chapmanb@50mail.com",
      description="Automated, distributed next-gen sequencing pipeline; includes Galaxy interaction",
      license="MIT",
      url="https://github.com/chapmanb/bcbb/tree/master/nextgen",
      namespace_packages=["bcbio"],
      packages=find_packages(),
      scripts=[
           'scripts/analyze_quality_recal.py',
           'scripts/automated_initial_analysis.py',
           'scripts/barcode_sort_trim.py',
           'scripts/bcbb_log_subscriber.py',
           'scripts/distributed_nextgen_pipeline.py',
           'scripts/bam_to_wiggle.py',
           'scripts/barcode_sort_trim.py',
           'scripts/monthly_billing_report.py',
           'scripts/monthly_billing_report.py',
           'scripts/illumina_finished_msg.py',
           'scripts/nextgen_analysis_server.py',
           'scripts/solexa_qseq_to_fastq.py',
           'scripts/store_finished_sqn.py',
           'scripts/upload_to_galaxy.py',
           'scripts/variant_effects.py',
           'scripts/count_barcodes.py'  # temporary name
          ],
      package_data={
          'config': ['*.yaml'],
          },
      install_requires=[
          "biopython >= 1.58",
          "Mako >= 0.3.6",
          "PyYAML >= 3.09",
          "amqplib >= 0.6.1",
          "Logbook >= 0.3",
          "pysam >= 0.4.1",
          "fabric >= 1.2",
          "pycrypto >= 2.1, < 2.4",
          "paramiko >= 1.7.7.1",
          "setuptools",
          "celery == 2.5",
          "nose >= 1.0.0",
          "gdata >= 2.0.14",
          "drmaa >= 0.5",
          "couchdb >= 0.8",
          "BeautifulSoup4"
      ])

os.system("git rev-parse --short --verify HEAD > ~/.bcbb_pipeline_version")
