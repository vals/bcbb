"""Create reports on google docs
"""
# import copy
import logbook
import time
import yaml

from bcbio.google import (_from_unicode, _to_unicode, get_credentials)
import bcbio.google.bc_metrics
import bcbio.google.qc_metrics
from bcbio.pipeline.qcsummary import RTAQCMetrics
from bcbio.pipeline.flowcell import Flowcell
from bcbio.log import create_log_handler
from bcbio.log import logger2 as log


def create_report_on_gdocs(fc_date, fc_name, run_info_yaml, dirs, config):
    """
    Create reports on gdocs containing both demultiplexed read counts and QC data
    """

    success = True
    try:
        # Parse the run_info.yaml file
        with open(run_info_yaml, "r") as fh:
            run_info = yaml.load(fh)

        # Get the gdocs account credentials
        encoded_credentials = get_credentials(config)
        if not encoded_credentials:
            log.warn("Could not find Google Docs account credentials in configuration. \
                No sequencing report was written")
            return False

        # Get the required parameters from the post_process.yaml configuration file
        gdocs = config.get("gdocs_upload", None)

        # Add email notification
        email = gdocs.get("gdocs_email_notification", None)
        smtp_host = config.get("smtp_host", "")
        smtp_port = config.get("smtp_port", "")
        log_handler = create_log_handler({'email': email, 'smtp_host': smtp_host, 'smtp_port': smtp_port}, True)

        with log_handler.applicationbound():

            # Inject the fc_date and fc_name in the email subject
            with logbook.Processor(lambda record: record.extra.__setitem__('run', "%s_%s" % (fc_date, fc_name))):

                try:
                    log.info("Started creating sequencing report on Google docs for %s_%s on %s" % (fc_date, fc_name, time.strftime("%x @ %X")))

                    # Get a flowcell object
                    fc = Flowcell(fc_name, fc_date, run_info, dirs.get("work", None))

                    # Get the GDocs demultiplex result file title
                    gdocs_dmplx_spreadsheet = gdocs.get("gdocs_dmplx_file", None)
                    # Get the GDocs QC file title
                    gdocs_qc_spreadsheet = gdocs.get("gdocs_qc_file", None)

                    # FIXME: Make the bc stuff use the Flowcell module
                    if gdocs_dmplx_spreadsheet is not None:
                        # Upload the data
                        success &= bcbio.google.bc_metrics.write_run_report_to_gdocs(fc, fc_date, fc_name, gdocs_dmplx_spreadsheet, encoded_credentials)
                    else:
                        log.warn("Could not find Google Docs demultiplex results file title in configuration. \
                            No demultiplex counts were written to Google Docs for %s_%s" \
                            % (fc_date, fc_name))

                    # Parse the QC metrics
                    qc = RTAQCMetrics(dirs.get("flowcell", None))
                    if gdocs_qc_spreadsheet is not None:
                        success &= bcbio.google.qc_metrics.write_run_report_to_gdocs(fc, qc, gdocs_qc_spreadsheet, encoded_credentials)
                    else:
                        log.warn("Could not find Google Docs QC file title in configuration. \
                            No QC data were written to Google Docs for %s_%s" \
                            % (fc_date, fc_name))

                    # Get the projects parent folder
                    projects_folder = gdocs.get("gdocs_projects_folder", None)

                    # Write the bc project summary report
                    if projects_folder is not None:
                        success &= create_project_report_on_gdocs(fc, qc, \
                            encoded_credentials, projects_folder)

                except Exception as e:
                    success = False
                    raise

                if success:
                    log.info("Sequencing report successfully created on \
                        Google docs for %s_%s on %s" \
                        % (fc_date, fc_name, time.strftime("%x @ %X")))
                else:
                    log.warn("Encountered exception when writing sequencing \
                        report for %s_%s to Google docs on %s" \
                        % (fc_date, fc_name, time.strftime("%x @ %X")))

    except Exception as e:
        success = False
        log.warn("Encountered exception when writing sequencing report to Google Docs: %s" % e)

    return success


def create_project_report_on_gdocs(fc, qc, encoded_credentials, gdocs_folder):
    """Upload the sample read distribution for a project to google docs"""

    success = True

    # Create a client class which will make HTTP requests with Google Docs server.
    client = bcbio.google.spreadsheet.get_client(encoded_credentials)
    doc_client = bcbio.google.document.get_client(encoded_credentials)

    # Get a reference to the parent folder
    parent_folder = bcbio.google.document.get_folder(doc_client, gdocs_folder)

    if not parent_folder:
        parent_folder_title = "root directory"
    else:
        parent_folder_title = _from_unicode(parent_folder.title.text)

    # Loop over the projects
    for project_name in fc.get_project_names():

        # Get a flowcell object containing just the data for the project
        project_fc = fc.prune_to_project(project_name)

        folder_name = project_name
        folder = bcbio.google.document.get_folder(doc_client, folder_name)
        if not folder:
            folder = bcbio.google.document.add_folder(doc_client, folder_name, parent_folder)
            log.info("Folder '%s' created under '%s'" \
                % (_from_unicode(folder_name), parent_folder_title))

        ssheet_title = project_name + "_sequencing_results"
        ssheet = bcbio.google.spreadsheet.get_spreadsheet(client, ssheet_title)
        if not ssheet:
            bcbio.google.document.add_spreadsheet(doc_client, ssheet_title)
            ssheet = bcbio.google.spreadsheet.get_spreadsheet(client, ssheet_title)
            ssheet = bcbio.google.document.move_to_folder(doc_client, ssheet, folder)
            ssheet = bcbio.google.spreadsheet.get_spreadsheet(client, ssheet_title)
            log.info("Spreadsheet '%s' created in folder '%s'" \
                % (_from_unicode(ssheet.title.text), _from_unicode(folder_name)))

        success &= bcbio.google.bc_metrics._write_project_report_to_gdocs( \
            client, ssheet, project_fc)
        success &= bcbio.google.bc_metrics._write_project_report_summary_to_gdocs( \
            client, ssheet)
        success &= bcbio.google.qc_metrics.write_run_report_to_gdocs( \
            project_fc, qc, ssheet_title, encoded_credentials)
        log.info("Sequencing results report written to spreadsheet '%s'" \
            % _from_unicode(ssheet.title.text))

    return success
