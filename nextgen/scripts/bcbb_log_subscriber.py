"""Catch log messages from a RabbitMQ queue as configured in post_process.yaml

\t%prog post_process_file.yaml [options]

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


def main(config_file):
    with open(config_file) as fh:
        config = yaml.load(fh)

    try:
        rmq_settings = config["rabbitmq_logging"]
    except KeyError:
        print("RabbitMQ logging not configured in {}".format(config_file))
        sys.exit()

    handlers = [NullHandler()]
    handlers.append(StderrHandler(bubble=True))
    handlers.append(FileHandler("bcbb.log", bubble=True))

    cdb_settings = config["couchdb_logging"]
    db_handler = DatabaseHandler(cdb_settings["couchdb_url"],
                                 backend=CouchDBBackend,
                                 db="logging_test",
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

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        sys.exit()

    config_file = args[0]

    main(config_file)
