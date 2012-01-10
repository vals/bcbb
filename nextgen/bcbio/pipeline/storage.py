"""Transfer raw files from finished NGS runs for backup and storage.
"""
import os

import yaml

from bcbio.log import logger
from bcbio.pipeline import log
from bcbio.log import create_log_handler
from bcbio.pipeline.config_loader import load_config
from bcbio.pipeline.transfer import remote_copy


def long_term_storage(remote_info, config_file):
    """Securely copy files from remote directory to the storage server.

    This requires ssh public keys to be setup so that no password entry
    is necessary, Fabric is used to manage setting up copies on the remote
    storage server.
    """
    import fabric.api as fabric
    import fabric.contrib.files as fabric_files
 
    logger.info("Copying run data over to remote storage: %s" % config["store_host"])
    logger.debug("The contents from AMQP for this dataset are:\n %s" % remote_info)
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    base_dir = config["store_dir"]
    try:
        protocol = config["transfer_protocol"]
    except KeyError:
        protocol = None
        pass

    fabric.env.host_string = "%s@%s" % \
    (config["store_user"], config["store_host"])
    remote_copy(remote_info, base_dir, protocol)
