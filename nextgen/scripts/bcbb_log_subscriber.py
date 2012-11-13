"""Catch log messages from a RabbitMQ queue as configured in post_process.yaml

\t%prog post_process.yaml [options]

Redirects the messages to handlers as specified by the options
"""
import sys
from optparse import OptionParser

import yaml

from logbook.queues import RabbitMQSubscriber
from logbook import NullHandler
from logbook import StderrHandler
from logbook import FileHandler
from logbook import NestedSetup
from logbook.more import DatabaseHandler
from logbook.more import CouchDBBackend


def main(config_file, **kwargs):
    with open(config_file) as fh:
        config = yaml.load(fh)

    try:
        rmq_settings = config["rabbitmq_logging"]
    except KeyError:
        print("RabbitMQ logging not configured in {}".format(config_file))
        sys.exit()

    handlers = [NullHandler()]
    if not kwargs["quiet"]:
        handlers.append(StderrHandler(bubble=True))

    if kwargs["filename"]:
        handlers.append(FileHandler(kwargs["filename"], bubble=True))

    if kwargs["log_db"]:
        try:
            cdb_settings = config["couchdb_logging"]
        except KeyError:
            print("CouchDB logging not configured in {}".format(config_file))
            sys.exit()

        db_handler = DatabaseHandler(cdb_settings["couchdb_url"],
                                     backend=CouchDBBackend,
                                     db=cdb_settings["database"],
                                     bubble=True)
        handlers.append(db_handler)

    setup = NestedSetup(handlers)

    print("Now waiting for log messages")
    with setup:
        subscriber = RabbitMQSubscriber(rmq_settings["url"], queue=rmq_settings["log_queue"])
        try:
            subscriber.dispatch_forever()

        except KeyboardInterrupt:
            print("\nLog subscriber shutting down")
            subscriber.close()

        except Exception:
            print("Log subscriber quit (unexpectedly)")

if __name__ == "__main__":
    parser = OptionParser(usage=__doc__)
    parser.add_option("-q", "--quiet",
                      action="store_true", dest="quiet", default=False,
                      help="Don't forward log records to StdErr")
    parser.add_option("-f", "--log-file", dest="filename", default=None,
                      help="Save log records in supplied file",
                      metavar="bcbb.log")
    parser.add_option("--db", "--couchdb", action="store_true", dest="log_db",
                      default=False,
                      help="Log to CouchDB database, configuration need to be \
                      supplied in the post_process.yaml")

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        sys.exit()

    config_file = args[0]

    kwargs = dict(quiet=options.quiet,
                  filename=options.filename,
                  log_db=options.log_db)

    main(config_file, **kwargs)
