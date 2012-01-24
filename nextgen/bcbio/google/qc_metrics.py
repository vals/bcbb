"""A module for writing lane QC data (typically from RTA) to google docs
"""

import os
import logbook 
from bcbio.log import create_log_handler
from bcbio.pipeline import log
from bcbio.pipeline.qcsummary import RTAQCMetrics
from bcbio.pipeline.flowcell import Flowcell
from bcbio.google import (_from_unicode,_to_unicode,get_credentials)
from bcbio.google.bc_metrics import (_write_to_worksheet,get_spreadsheet)

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
    
    # Get a flowcell object 
    fc = Flowcell(fc_name,fc_date,run_info.get("details",[]))
    
    write_run_report_to_gdocs(fc,qc,gdocs_spreadsheet,encoded_credentials)
    
    # Add email notification
#    email = gdocs.get("gdocs_email_notification",None)
#    log_handler = create_log_handler({'email': email},"")
#    with log_handler.applicationbound():
#        
#        # Inject the fc_date and fc_name in the email subject
#        with logbook.Processor(lambda record: record.extra.__setitem__('run', "%s_%s" % (fc_date,fc_name))):
#        
#            # Get the barcode statistics. Get a deep copy of the run_info since we will modify it
#            bc_metrics = get_bc_stats(fc_date,fc_name,work_dir,copy.deepcopy(run_info))
#            
#            # Upload the data
#            write_run_report_to_gdocs(fc_date,fc_name,bc_metrics,gdocs_spreadsheet,encoded_credentials)
#            
#            # Get the projects parent folder
#            projects_folder = gdocs.get("gdocs_projects_folder",None)
#            
#            # Write the bc project summary report
#            if projects_folder:
#                write_project_report_to_gdocs(fc_date,fc_name,bc_metrics,encoded_credentials,projects_folder)

def write_run_report_to_gdocs(fc,qc,ssheet_title,encoded_credentials,wsheet_title=None):
    
    # Connect to google and get the spreadsheet
    client, ssheet = get_spreadsheet(ssheet_title,encoded_credentials)
    if not client or not ssheet:
        return False

    qc_metrics = qc.metrics()
    qc_stats = qc.getQCstats()
    
    # Create the header row
    header = [["Date"],["Flowcell"],["Lane"],["Description"]]
    # Add the metric labels
    metric_lbl = []
    for metric in qc_metrics:
        if metric[0].endswith('_sd'):
            continue
        if metric[0] not in qc_stats:
            continue
        metric_lbl.append(metric[0])
        header.append([metric[1]])
        
    # Iterate over the lanes of the flowcell and collect the data
    rows = []
    for lane in fc.get_lanes():
        # First the meta data
        row = [fc.get_fc_date(),fc.get_fc_name(),lane.get_name(),lane.get_description()]
        # Then the QC data
        for metric in metric_lbl:
            value = qc_stats[metric]
            sd_key = "%s_sd" % metric
            if sd_key in qc_stats:
                sd = qc_stats[sd_key]
            else:
                sd = None
            cell = ""
            for read in value.keys():
                val = "%s" % value[read][lane.get_name()]
                if sd is not None:
                    val += " +/- %s" % sd[read][lane.get_name()]
                if metric.find('cluster') >= 0:
                    cell = val
                    break
                cell += "%s: %s\n" % (read,val)
            row.append(cell)
        rows.append(row)
    
    success = True
    if wsheet_title is None:
        wsheet_title = "%s_%s_%s" % (fc.get_fc_date(),fc.get_fc_name(),"QC")
        success &= _write_to_worksheet(client,ssheet,wsheet_title,rows,header,False)

    return success

                
        
    
    