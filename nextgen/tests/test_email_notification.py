"""
A test for sending out email notifications based on a logger
"""
import os
import unittest
import time
from bcbio.log import (create_log_handler,logger2)
from bcbio.pipeline.config_loader import load_config
import logbook

class TestEmailNotification(unittest.TestCase):
    
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data", "automated")
        self.config_file = os.path.join(self.data_dir,"post_process.yaml")
        assert os.path.exists(self.config_file), "Could not locate required configuration file %s" % self.config_file
 
    def _get_log_handler(self,config):
        log_handler = create_log_handler(config,True)
        return log_handler

    def _log_messages(self, log_handler, subject="Test email"):
        try:
            with log_handler.applicationbound():
                with logbook.Processor(lambda record: record.extra.__setitem__('run', subject)):
                    logger2.debug("DEBUG record test generated @ %s" % time.strftime("%x - %X"))
                    logger2.info("INFO record test generated @ %s" % time.strftime("%x - %X"))
                    logger2.notice("NOTICE record test generated @ %s" % time.strftime("%x - %X"))
                    logger2.warning("WARNING record test generated @ %s" % time.strftime("%x - %X"))
                    logger2.error("ERROR record test generated @ %s" % time.strftime("%x - %X"))
                    logger2.critical("CRITICAL record test generated @ %s" % time.strftime("%x - %X"))
        except Exception as e:
            return e
        return None
              
    def test_1_notify(self):
        config = load_config(self.config_file)
        if not "email" in config:
            print "No email configured, skipping test!"
            return
        log_handler = self._get_log_handler(config)
        result = self._log_messages(log_handler,"Pipeline notification test email @ %s" % time.strftime("%x - %X"))
        assert result is None, "%s" % result
        
    def test_2_report_notification(self):
        config = load_config(self.config_file)
        if not "gdocs_upload" in config or not "gdocs_email_notification" in config["gdocs_upload"]:
            print "Google docs email reporting not configured, skipping test!"
            return
        
        config["email"] = config["gdocs_upload"]["gdocs_email_notification"]
        
        log_handler = self._get_log_handler(config)
        result = self._log_messages(log_handler,"Google Docs report notification test email @ %s" % time.strftime("%x - %X"))
        assert result is None, "%s" % result
