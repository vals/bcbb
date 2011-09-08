#!/usr/bin/env python
"""Upload the information on the number of reads per barcode after demultiplex generated by the analysis pipeline to a spreadsheet on Google docs.

Given a tab-separated file of data generated by the analysis pipeline (typically WORKDIR/demultiplexed_read_counts.txt),
the title of a (pre-existing) spreadsheet on a GoogleDocs account (typically "Demultiplexed read counts") and the account 
username and password, will create a new worksheet in the spreadsheet and enter the data into it. 

Usage:
    demultiplex_upload.py <demultiplex report file> <spreadsheet title> <Google user name> <Google password>

The columns in the report file are assumed to be:

    date
    fcid
    lane_id
    description
    barcode_id
    barcode_name
    barcode_sequence
    demultiplexed_reads
"""
import os
import sys
import csv
import base64
import yaml
from optparse import OptionParser
from bcbio.pipeline import log

import gdata.spreadsheet.service
import gdata.docs.service

# The structure of the demultiplex result file on the form {title: column index}
header = {
          'date': 0,
          'flowcell_id': 1,
          'lane_id': 2,
          'description': 3,
          'barcode_id': 4,
          'barcode_name': 5,
          'barcode_sequence': 6,
          'barcode_type': 7,
          'demultiplexed_reads': 8
          }

def decode_credentials(credentials):
    
    if not credentials:
        return None
    
    # Split the username and password
    return base64.b64decode(credentials).split(':',1);
    

def main(dmplx_report_file, gdocs_spreadsheet_title, credentials):

    # Get the local demultiplex result file
    if not os.path.exists(dmplx_report_file):
        log.info("Could not find local demultiplex results file %s. No demultiplex counts were written to Google Docs" % dmplx_report_file)
        return
    
    # Split the username and password
    credentials = decode_credentials(credentials)
    if not credentials:
        log.info("Could not find GDocs credentials in config. No demultiplex counts were written to Google Docs")
        return
    
    # Get the GDocs demultiplex result file title
    if not gdocs_spreadsheet_title:
        log.info("Could not find Google Docs demultiplex results file title in config. No demultiplex counts were written to Google Docs")
        return
    
    upload_demultiplex_data(dmplx_report_file,gdocs_spreadsheet_title,credentials)
    
    
def upload_demultiplex_data(dmplx_report_file, gdocs_spreadsheet_title, credentials):
       
    # Create a client class which will make HTTP requests with Google Docs server.
    client = gdata.spreadsheet.service.SpreadsheetsService()
    client.email = credentials[0]
    client.password = credentials[1]
    client.source = 'demultiplex_upload.py'
    client.ProgrammaticLogin()
    
    # Create a query that restricts the spreadsheet feed to documents having the supplied title
    q = gdata.spreadsheet.service.DocumentQuery(params={'title':gdocs_spreadsheet_title, 'title-exact':'true'})
    # Query the server for an Atom feed containing a list of your documents.
    feed = client.GetSpreadsheetsFeed(query=q)
    
    # Check that we got a result back
    if len(feed.entry) == 0:
        log.info("No document with specified title '%s' found in GoogleDocs repository" % gdocs_spreadsheet_title)
        return
    
    # If we get more than one feed item back, will have to implement some way of resolving these
    if len(feed.entry) > 1:
        log.info("More than one document match the specified title '%s', will have to implement some way of resolving these!" % gdocs_spreadsheet_title)
        return
    
    # Get the matching spreadsheet entry from the feed    
    spreadsheet = feed.entry[0]
    ss_key = spreadsheet.id.text.split('/')[-1]
    log.info("Found spreadsheet matching the supplied title: '%s'" % spreadsheet.title.text)
    
    # Parse the result file to get the necessary parameters to determine worksheet names etc.
    with open(dmplx_report_file,"rb") as f:
        csvr = csv.reader(f,dialect='excel-tab')
        
        # Read the entire file
        rows = []
        for row in csvr:
            rows.append(row)
        
        # Check that we got any rows
        if len(rows) == 0:
            log.info("Did not read any data from input file '%s'" % result_file)
            return
        
        # Get the date and flowcell id and base the name of the worksheet on these
        date = rows[0][header.get('date',0)]
        fcid = rows[0][header.get('flowcell_id',1)]
        
        # Check if there already is a worksheet present matching the data
        ws_title = "%s_%s" % (date,fcid)
        q = gdata.spreadsheet.service.DocumentQuery(params={'title':ws_title})
        feed = client.GetWorksheetsFeed(key=ss_key,query=q)
        # If there already exists one or more matching worksheet(s), add an enumerating suffix to the ws title
        if len(feed.entry) > 0:
            ws_title += "(%s)" % (len(feed.entry) + 1)
        
        # Create a new worksheet for this run and add it to the end of the spreadsheet
        log.info("Adding a new worksheet, '%s', to the spreadsheet" % ws_title)
        ws = client.AddWorksheet(ws_title,len(rows)+1,len(header),ss_key)
        if ws is None:
            log.info("Could not add a worksheet '%s' to spreadsheet with key '%s'" % (ws_title,ss_key))
            return
        
        ws_id = ws.id.text.split('/')[-1]

        log.info("Adding data to the '%s' worksheet" % ws_title)
        # First, print the header
        for val, j in header.items():
             client.UpdateCell(1,j+1,val,ss_key,ws_id)
        
        # Iterate over the rows and add the data to the worksheet
        i=2
        for row in rows:
            j=1
            for val in row:
                client.UpdateCell(i,j,val,ss_key,ws_id)
                j += 1
            i += 1
 


if __name__ == "__main__":
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) != 3:
        print "Incorrect arguments"
        print __doc__
        sys.exit()
    kwargs = dict()
    main(*args, **kwargs)