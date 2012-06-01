"""
Create statusdb report
"""

import logbook
import time
import yaml
import couchdb

from bcbio.log import create_log_handler
from bcbio.log import logger2 as log
from bcbio.qc import FlowcellQCMetrics

def _save_obj(db, obj, url):
    dbobj = db.get(obj.get_db_id())
    if dbobj is None:
        log.info("Creating entity type %s with id %s in url %s" % (obj["entity_type"], obj.get_db_id(), url))
        db.save(obj)
    else:
        obj["_rev"] = dbobj.get("_rev")
        if obj != dbobj:
            log.info("Updating %s object with id %s in url %s" % (obj["entity_type"], obj.get_db_id(), url))
            db.save(obj)
        else:
            log.info("Object %s already present in url %s and not in need of updating" % (obj.get_db_id(), url))
    return True


## Do the actual reporting
## This part could be run through a handler
def report_to_statusdb(fc_name, fc_date, run_info_yaml, dirs, config):
    """
    Create statusdb report on a couchdb server.

    A FlowcellQCMetrics object holds information about a flowcell. QC
    results are stored at the flowcell level and sample level
    depending on analysis. Lane level QC data are stored in the
    FlowcellQCMetrics object.
    """
    success = True
    try:
        statusdb_config = config.get("statusdb", None)
        if statusdb_config is None:
            log.info("Could not find statusdb section in configuration. No statusdb reporting will be done")
            return False
        statusdb_couchdb =  statusdb_config.get("url", None)
        if statusdb_couchdb is None:
            log.warn("No url field found in statusdb configuration section.")
            return False

        # Add email notification
        email = statusdb_config.get("statusdb_email_notification", None)
        smtp_host = config.get("smtp_host", "")
        smtp_port = config.get("smtp_port", "")
        log_handler = create_log_handler({'email': email, 'smtp_host': smtp_host, 'smtp_port': smtp_port}, True)

        with log_handler.applicationbound():
            with logbook.Processor(lambda record: record.extra.__setitem__('run', "%s_%s" % (fc_date, fc_name))):
                log.info("Started creating QC Metrics report on statusdb for %s_%s on %s" % (fc_date, fc_name, time.strftime("%x @ %X")))

                # Create object and parse all available metrics; no checking
                # is currently done for missing files
                try:
                    qc_obj = FlowcellQCMetrics(dirs.get("work", None), dirs.get("flowcell", None), url=statusdb_couchdb)
                    qc_obj.parse_run_info_yaml()
                    qc_obj.read_picard_metrics()
                    qc_obj.read_fastqc_metrics()
                    qc_obj.parse_filter_metrics()
                    qc_obj.parse_fastq_screen()
                    qc_obj.parse_bc_metrics()
                    qc_obj.parse_illumina_metrics()
                except:
                    qc_obj = None
                # FIXME: error checking!
                if qc_obj is not None:
                    try:
                        # Save data at a sample level
                        log.info("Connecting to server at %s" % qc_obj.url)
                        try:
                            couch = couchdb.Server(url="http://%s" % qc_obj.url)
                        except:
                            log.warn("Connecting to server at %s failed" % qc_obj.url)
                        log.info("Connecting to server at %s succeeded" % qc_obj.url)
                        db=couch['qc']
                        # Save samples
                        for s in qc_obj.sample.keys():
                            obj = qc_obj.sample[s]
                            log.info("Saving sample %s" % obj.name())
                            _save_obj(db, obj, qc_obj.url)
                        # Save flowcell object
                        _save_obj(db, qc_obj, qc_obj.url)
                    except Exception as e:
                        success = False
                else:
                    log.warn("Couldn't populate FlowcellQCMetrics object. No QC data written to statusdb for %s_%s" % (fc_date, fc_name))
            if success:
                log.info("QC Metrics report successfully written to statusdb for %s_%s on %s" \
                         % (fc_date, fc_name, time.strftime("%x @ %X")))
            else:
                log.warn("Encountered exception when writing to statusdb for %s_%s on %s" \
                         % (fc_date, fc_name, time.strftime("%x @ %X")))

    except Exception as e:
        success = False
        log.warn("Encountered exception when writing QC metrics to statusdb: %s" % e)

    return success
