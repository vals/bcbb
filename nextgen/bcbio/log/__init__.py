"""Utility functionality for logging.
"""
import os
import sys
from datetime import datetime

import logging

from bcbio import utils

LOG_NAME = "nextgen_pipeline"

logger = logging.getLogger(LOG_NAME)


def setup_logging(config):
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        formatter = logging.Formatter('[%(asctime)s] %(message)s')
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        log_dir = config.get("log_dir", None)
        if log_dir:
            logfile = os.path.join(utils.safe_makedir(log_dir),
                                   "{0}.log".format(LOG_NAME))
            handler = logging.FileHandler(logfile)
            handler.setFormatter(formatter)
            logger.addHandler(handler)

import logbook
logger2 = logbook.Logger(LOG_NAME)


def create_log_handler(config, batch_records=False):
    log_dir = config.get("log_dir", None)
    email = config.get("email", None)

    if log_dir:
        utils.safe_makedir(log_dir)
        handler = logbook.FileHandler(os.path.join(log_dir, "%s.log" % LOG_NAME))
    else:
        handler = logbook.StreamHandler(sys.stdout)

    if email:        
        smtp_host = config.get("smtp_host",None)
        smtp_port = config.get("smtp_port",25)
        if smtp_host is not None: smtp_host = [smtp_host,smtp_port]
        
        email = email.split(",")
        handler = logbook.MailHandler(email[0], email, server_addr=smtp_host,
                                      format_string=u'''Subject: [BCBB pipeline] {record.extra[run]} \n\n {record.message}''',
                                      level='INFO', bubble=True)
    if batch_records:
        handler = logbook.handlers.GroupHandler(handler)
        
    return handler

def utc_time():
    """
    Make an utc_time with appended 'Z'
    Borrowed from scilifelab.utils.timestamp
    """
    return str(datetime.utcnow()) + 'Z'
    
