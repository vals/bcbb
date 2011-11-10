#!/usr/bin/env python
"""Functions for getting barcode statistics from demultiplexing"""

import os
import re
import copy
import glob
from bcbio.utils import UnicodeReader
import bcbio.google.connection
import bcbio.google.document
import bcbio.google.spreadsheet
from bcbio.google import (_from_unicode,_to_unicode)
from bcbio.pipeline import log


# The structure of the demultiplex result
BARCODE_STATS_HEADER = [
                 ['Project name','project_name'],
                 ['Date','date'],
                 ['Flowcell','flowcell_id'],
                 ['Lane','lane'],
                 ['Description','description'],
                 ['Sample name','sample_name'],
                 ['bcbb internal barcode index','barcode_id'],
                 ['Barcode name','name'],
                 ['Barcode sequence','sequence'],
                 ['Barcode type','barcode_type'],
                 ['Demultiplexed read (pair) count','barcode_read_count'],
                 ['Demultiplexed read (pair) count (millions)','barcode_read_count_millions'],
                 ['Comment','comment']
                ]
  
# The structure of the sequencing result
SEQUENCING_RESULT_HEADER = [
                 ['Sample name','sample_name'],
                 ['Run','run'],
                 ['Lane','lane'],
                 ['Read (pair) count','read_count'],
                 ['Read (pair) count (millions)','read_count_millions'],
                 ['Comment','comment'],
                 ['Pass','pass']
                ]
  
def create_bc_report_on_gdocs(fc_date, fc_name, work_dir, run_info, config):
    """Get the barcode read distribution for a run and upload to google docs"""
    
    # Get the required parameters from the post_process.yaml configuration file
    gdocs = config.get("gdocs_upload",None)
    if not gdocs:
        log.info("No GDocs upload section specified in config file, will not upload demultiplex data")
        return
    
    # Get the GDocs demultiplex result file title
    gdocs_spreadsheet = gdocs.get("gdocs_dmplx_file",None)
    if not gdocs_spreadsheet:
        log.warn("Could not find Google Docs demultiplex results file title in config. No demultiplex counts were written to Google Docs")
        return
    
    # Get the account credentials
    encoded_credentials = ""
    encoded_credentials_file = gdocs.get("gdocs_credentials",None)
    if not encoded_credentials_file:
        log.warn("Could not find Google Docs account credentials. No demultiplex report was written")
        return
    # Check if the credentials file exists
    if not os.path.exists(encoded_credentials_file):
        log.warn("The Google Docs credentials file could not be found. No demultiplex data was written")
        return
    with open(encoded_credentials_file) as fh:
        encoded_credentials = fh.read().strip()
    
    # Get the barcode statistics. Get a deep copy of the run_info since we will modify it
    #bc_metrics = get_bc_stats(fc_date,fc_name,work_dir,copy.deepcopy(run_info))
    fc = Flowcell(fc_name,fc_date,run_info,work_dir)
    
    # Upload the data
    write_run_report_to_gdocs(fc,fc_date,fc_name,bc_metrics,gdocs_spreadsheet,encoded_credentials)
    
    # Get the projects parent folder
    projects_folder = gdocs.get("gdocs_projects_folder",None)
    
    # Write the bc project summary report
    if projects_folder:
        write_project_report_to_gdocs(fc,fc_date,fc_name,bc_metrics,encoded_credentials,projects_folder)

def get_spreadsheet(ssheet_title,encoded_credentials):
    """Connect to Google docs and get a spreadsheet"""
    
    # Convert the spreadsheet title to unicode
    ssheet_title = _to_unicode(ssheet_title)
    
    # Create a client class which will make HTTP requests with Google Docs server.
    client = bcbio.google.spreadsheet.get_client()
    bcbio.google.connection.authenticate(client,encoded_credentials)
    
    # Locate the spreadsheet
    ssheet = bcbio.google.spreadsheet.get_spreadsheet(client,ssheet_title)
    
    # Check that we got a result back
    if not ssheet:
        log.warn("No document with specified title '%s' found in GoogleDocs repository" % ssheet_title)
        return (None,None)
    
    log.info("Found spreadsheet matching the supplied title: '%s'" % (ssheet.title.text))
    
    return (client,ssheet)
               
def write_project_report_to_gdocs(flowcell,encoded_credentials,gdocs_folder=""):
    """Upload the sample read distribution for a project to google docs"""
    
    # Create a client class which will make HTTP requests with Google Docs server.
    client = bcbio.google.spreadsheet.get_client(encoded_credentials)
    doc_client = bcbio.google.document.get_client(encoded_credentials)
    
    # Get a reference to the parent folder
    parent_folder = bcbio.google.document.get_folder(doc_client,gdocs_folder)
    
    # Get the projects on the flowcell
    projects = flowcell.get_project_names()
    
    # Loop over the projects and write the project summary for each
    for project_name in projects:
        
        pruned_fc = flowcell.prune_to_project(project_name)
        ssheet_title = project_name + "_sequencing_results"
        ssheet = bcbio.google.spreadsheet.get_spreadsheet(client,ssheet_title)
        if not ssheet:
            bcbio.google.document.add_spreadsheet(doc_client,ssheet_title)
            ssheet = bcbio.google.spreadsheet.get_spreadsheet(client,ssheet_title)
    
        _write_project_report_to_gdocs(client,ssheet,pruned_fc)
        _write_project_report_summary_to_gdocs(client,ssheet)
        
        # Just to make it look a bit nicer, remove the default 'Sheet1' worksheet
        wsheet = bcbio.google.spreadsheet.get_worksheet(client,ssheet,'Sheet 1')
        if wsheet:
            client.DeleteWorksheet(wsheet)
            
        folder_name = project_name
        folder = bcbio.google.document.get_folder(doc_client,folder_name)
        if not folder:
            log.info("creating folder '%s'" % folder_name)
            folder = bcbio.google.document.add_folder(doc_client,folder_name,parent_folder)
            
        ssheet = bcbio.google.document.move_to_folder(doc_client,ssheet,folder)
        log.info("'%s' spreadsheet written to folder '%s'" % (ssheet.title.text,folder_name))
        

def _write_project_report_to_gdocs(client, ssheet, flowcell):

    # Get the spreadsheet if it exists
    # Otherwise, create it
    wsheet_title = "%s_%s" % (flowcell.get_fc_date(),flowcell.get_fc_name())
    
    # Flatten the project_data structure into a list
    rows = []
    for sample in flowcell.get_samples():
        row = (sample.get_name(),wsheet_title,sample["lane"],sample.get_read_count(),sample.get_rounded_read_count(),sample.get("comment",""),"")
        rows.append(row)
    
    # Write the data to the worksheet
    return _write_to_worksheet(client,ssheet,wsheet_title,rows,SEQUENCING_RESULT_HEADER,False)
    
def _write_project_report_summary_to_gdocs(client, ssheet):
    """Summarize the data from the worksheets and write them to a "Summary" worksheet"""
    
    # Summary data
    summary_data = {}
    
    # Get the list of worksheets in the spreadsheet
    wsheet_feed = bcbio.google.spreadsheet.get_worksheets_feed(client,ssheet)
    # Loop over the worksheets and parse the data from the ones that contain flowcell data
    for wsheet in wsheet_feed.entry:
        if not re.match(r'^\d{2}[01]\d[0-3]\d_\S+XX$',wsheet.title.text):
            continue
        wsheet_data = bcbio.google.spreadsheet.get_cell_content(client,ssheet,wsheet,'2')
        # Add the results from the worksheet to the summarized data    
        for wsheet_row in wsheet_data:
            sample_name = wsheet_row[0]
            summary_row = summary_data.get(sample_name,None)
            if not summary_row:
                summary_row = [""]*len(wsheet_row)
                summary_row[0] = sample_name
                for i in range(1,len(wsheet_row)):
                    summary_row[i] = []
                summary_data[sample_name] = summary_row 
            for i in range(1,len(wsheet_row)):
                summary_row[i].append(wsheet_row[i])
    
    # Concatenate the data in the summary structure
    for sample in summary_data.values():
        # Concatenate the non-number columns using ';' as separator
        for i in (1,2,5,6):
            sample[i] = ";".join(sample[i])
        
        # Sum up the read counts
        counts = sample[3]
        sum = 0
        for count in counts:
            try:
                c = int(count)
            except ValueError:
                c = 0
            sum += c
        # Count the millions
        msum = round(sum/1000000.,2)
        sample[3] = unicode(sum)
        sample[4] = unicode(msum)
        
    # Write the summary data to the worksheet
    wsheet_title = "Summary"
    return _write_to_worksheet(client,ssheet,wsheet_title,summary_data.values(),SEQUENCING_RESULT_HEADER,False)
            

def write_run_report_to_gdocs(fc, fc_date, fc_name, bc_metrics, ssheet_title, encoded_credentials, wsheet_title=None, append=False, split_project=False):
    """Upload the barcode read distribution for a run to google docs"""
    
    # Connect to google and get the spreadsheet
    client, ssheet = get_spreadsheet(ssheet_title,encoded_credentials)
    if not client or not ssheet:
        return False
    
    # Get the projects in the run
    projects = fc.get_project_names()
    log.info("The run contains data from: '%s'" % "', '".join(projects))
    
    # If we will split the worksheet by project, use the project names as worksheet titles
    success = True
    if split_project:
        # Filter away the irrelevent project entries and write the remaining to the appropriate worksheet
        for project in projects:
            pruned_fc = fc.prune_to_project(project)
            success &= _write_to_worksheet(client,ssheet,project,pruned_fc.to_rows(),BARCODE_STATS_HEADER,append)
            
    # Else, set the default title of the worksheet to be a string of concatenated date and flowcell id
    else:
        if wsheet_title is None:
            wsheet_title = "%s_%s" % (fc_date,fc_name)
        success &= _write_to_worksheet(client,ssheet,wsheet_title,fc.to_rows(),BARCODE_STATS_HEADER,append)

    return success

def _write_to_worksheet(client,ssheet,wsheet_title,rows,header,append):
    """Generic method to write a set of rows to a worksheet on google docs"""
    
    # Convert the worksheet title to unicode
    wsheet_title = _to_unicode(wsheet_title)
    
    # Add a new worksheet, possibly appending or replacing a pre-existing worksheet according to the append-flag
    wsheet = bcbio.google.spreadsheet.add_worksheet(client,ssheet,wsheet_title,len(rows)+1,len(header),append)
    if wsheet is None:
        log.info("Could not add a worksheet '%s' to spreadsheet '%s'" % (wsheet_title,ssheet.title.text))
        return False
    
    # Write the data to the worksheet
    log.info("Adding data to the '%s' worksheet" % (wsheet_title))
    return bcbio.google.spreadsheet.write_rows(client,ssheet,wsheet,[col_header[0] for col_header in header],rows)
     
    
    