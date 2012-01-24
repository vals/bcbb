"""A module for writing lane QC data (typically from RTA) to google docs
"""

import os
import logbook 
from bcbio.log import create_log_handler
from bcbio.pipeline.qcsummary import RTAQCMetrics

def create_qc_report_on_gdocs(base_dir, run_info, config):
    """Get the RTA QC metrics for a run and upload to google docs"""
    
    encoded_credentials = get_credentials(config)
    if not encoded_credentials:
        log.warn("Could not find Google Docs account credentials. No QC report was written")
        return
    
    # Get the required parameters from the post_process.yaml configuration file
    gdocs = config.get("gdocs_upload",None)
    
    # Get the GDocs demultiplex result file title
    gdocs_spreadsheet = gdocs.get("gdocs_qc_file",None)
    if not gdocs_spreadsheet:
        log.warn("Could not find Google Docs QC file title in config. No QC data were written to Google Docs")
        return

    # Parse the QC metrics
    qc = RTAQCMetrics(base_dir)
    
    # Get some metadata from the configuration
    run_cfg = qc.configuration()
    fc_date = run_cfg.date()
    fc_name = run_cfg.flowcell()
    
    
    # Add email notification
    email = gdocs.get("gdocs_email_notification",None)
    log_handler = create_log_handler({'email': email},"")
    with log_handler.applicationbound():
        
        # Inject the fc_date and fc_name in the email subject
        with logbook.Processor(lambda record: record.extra.__setitem__('run', "%s_%s" % (fc_date,fc_name))):
        
            # Get the barcode statistics. Get a deep copy of the run_info since we will modify it
            bc_metrics = get_bc_stats(fc_date,fc_name,work_dir,copy.deepcopy(run_info))
            
            # Upload the data
            write_run_report_to_gdocs(fc_date,fc_name,bc_metrics,gdocs_spreadsheet,encoded_credentials)
            
            # Get the projects parent folder
            projects_folder = gdocs.get("gdocs_projects_folder",None)
            
            # Write the bc project summary report
            if projects_folder:
                write_project_report_to_gdocs(fc_date,fc_name,bc_metrics,encoded_credentials,projects_folder)

def write_run_report_to_gdocs(fc_date,fc_name,qc,gdocs_spreadsheet,encoded_credentials):
    
    