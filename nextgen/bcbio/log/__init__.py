"""Utility functionality for logging.
"""
import os
import sys

import logbook
from bcbio import utils

def create_log_handler(config, log_name):
    log_dir = config.get("log_dir", None)
    email = config.get("email", None)
    
    if log_dir:
        utils.safe_makedir(log_dir)
        handler = logbook.FileHandler(os.path.join(log_dir, "%s.log" % log_name))
    else:
        handler = logbook.StreamHandler(sys.stdout)
        
    if email:
        smtp_host = config.get("smtp_host",None)
        smtp_port = config.get("smtp_port",25)
        if smtp_host is not None: smtp_host = [smtp_host,smtp_port]
        
        email = email.split(",")
        handler = logbook.MailHandler(email[0], email, server_addr=smtp_host,  
                                      format_string=u'''Subject: [BCBB pipeline] {record.extra[run]} \n\n {record.message}''',
                                      level='INFO', bubble = True)
    return handler
