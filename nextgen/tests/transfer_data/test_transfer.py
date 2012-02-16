"""Testing using config files for transfer settings.
"""
import os
import fabric.api as fabric
import fabric.contrib.files as fabric_files
import time
from bcbio.log import logger
from bcbio.pipeline.config_loader import load_config

from bcbio.pipeline.toplevel import _copy_from_sequencer
from bcbio.pipeline.storage import _copy_for_storage

from illumina_finished_msg import _files_to_copy

test_dir_structure = \
{"999999_XX999_9999_AC9999ACXX": [ \
    {"Data": [ \
        {"Intensities": [ \
            {"BaseCalls": [ \
                {"fastq": [ \
                    "1_999999_AC9999ACXX_1_fastq.txt", \
                    "1_999999_AC9999ACXX_2_fastq.txt", \
                    "2_999999_AC9999ACXX_1_fastq.txt"
                ]}, \
                {"Plots": [ \
                    "s_1_1101_all.png", \
                    "s_1_1102_all.png", \
                    "s_1_1103_all.png"
                ]}, \
                "All.htm", \
                "BustasrdSummary.xsl", \
                "BustasrdSummary.xml" \
            ]} \
        ]}, \
        {"reports": [ \
            {"level1": [ \
                {"level2": [ \
                    "file21", \
                    "file22" \
                ]}, \
                "file11", \
                "file12" \
            ]}, \
            "NumClusters_Chart.png", \
            "NumClusters_Chart.xml", \
            "NumPassedFilter25_Chart.png" \
        ]}, \
        {"Status_Files": [ \
            "ByCycleFrame.htm", \
            "ByCycle.htm", \
            "ByCycle.js" \
        ]}, \
        "Status.htm" \
    ]}, \
    {"InterOp": [ \
        "ControlMetricsOut.bin", \
        "CorrectedIntMetricsOut.bin", \
        "ErrorMetricsOut.bin" \
    ]}, \
    "RunInfo.xml", \
    "run_info.yaml", \
    "runParameters.xml" \
]}


def _remove_transferred_files(remote_info, config):
    """Remove the files transferred in a previous test.
    """
    copy_to = os.path.realpath("../transfer_data/copy_to")
    with fabric.settings(host_string="%s@%s" % \
         (config["store_user"], config["store_host"])):
        rm_str = "rm -r %s/%s" % \
         (copy_to, os.path.split(remote_info["directory"])[1])
        logger.debug(rm_str)
        fabric.run(rm_str)


def make_dirs_or_files(argument, residue=""):
    """Function for making the directories and files as defined in
    a nesting of dictionaries, lists and strings.
    """
    if type(argument) == dict:
        residue += argument.keys()[0] + "/"
        a_list = argument.values()[0]
        make_dirs_or_files(a_list, residue)
    if type(argument) == list:
        for elt in argument:
            make_dirs_or_files(elt, residue)
    if type(argument) == str:
        try:
            os.makedirs(residue)
        except OSError as e:
            if e.errno == 17:  # The directory already exists
                pass
            else:
                raise e

        open(residue + argument, 'w').close()


def print_file_list():
    for root, dirs, files in os.walk("999999_XX999_9999_AC9999ACXX"):
        for f in files:
            print(root + "/" + f)


def perform_transfer(transfer_function, protocol_config, \
        remove_before_copy=True, should_overwrite=False):
    """Sets up dictionaries simulating loaded remote_info and config
    from various sources. Then test transferring files with the function
    using the standard setting.

    Note that there need to be passphrase free SSH between this machine
    and what is set up as store_host in post_process.yaml.

    The file post_process.yaml need the information "store_user" and
    "store_host".

    transfer_function : function - The function to use for transferring
    the files.
        This function should accept the two dictionaries config and
        remote_info as parameters.
    protocol_config : dictionary - containing the key "transfer_protocol"
    which give the string
    of a transfer protocol.
    """
    config = load_config("transfer_test_post_process.yaml")

    for file_dir in "copy_to", "to_copy":
        if not os.path.isdir(file_dir):
            os.mkdir(file_dir)

    store_dir = os.path.realpath("copy_to")
    config["store_dir"] = store_dir
    config.update(protocol_config)
    copy_dir = os.path.realpath("to_copy")

    # Generate test files
    test_files = []
    make_dirs_or_files(test_dir_structure, "to_copy/")
    base = test_dir_structure.keys()[0]
    for root, dirs, files in os.walk("to_copy" + "/" + base):
        for f in files:
            test_files.append(root[8:] + "/" + f)

    test_data = {}
    for test_file in test_files:
        with open("to_copy/" + test_file, 'w+') as file_to_write:
            # We just use the current processor time as test data, the
            # important part is that it will be different enough between
            # test just so we know we are not comparing with files copied
            # in a previous test during the assertion.
            test_data[test_file] = str(time.clock())
            file_to_write.write(test_data[test_file])

    remote_info = {}
    remote_info["directory"] = copy_dir
    remote_info["user"] = config["store_user"]
    remote_info["hostname"] = config["store_host"]

    files_to_copy = _files_to_copy("to_copy/" + base)
    files_to_copy = list(set(sum(files_to_copy, [])))
    files_to_copy = map(lambda fpath: base + "/" + fpath, files_to_copy)

    # Trim out the unavailable files to copy
    unavailable_files = []
    for cfl in files_to_copy:
        available_for_copy = False
        for tfl in test_files:
            if cfl in tfl:
                available_for_copy = True

        if not available_for_copy:
            unavailable_files.append(cfl)

    for ufl in unavailable_files:
        files_to_copy.remove(ufl)

    remote_info["to_copy"] = files_to_copy

    # Perform copy with settings
    with fabric.settings(host_string="%s@%s" % \
         (remote_info["user"], remote_info["hostname"])):

        if fabric_files.exists("%s/%s" % \
           (store_dir, os.path.split(copy_dir)[1])) and remove_before_copy:
            _remove_transferred_files(remote_info, config)

        # Copy
        transfer_function(remote_info, config)

    # Check if the copy succeeded
    for test_file, test_value in test_data.items():
        test_file_path = "%s/%s" % (store_dir, "to_copy/" + test_file)
        # Did the files get copied correctly
        assert os.path.isfile(test_file_path), "File not copied: %s" % test_file
        if os.path.isfile(test_file_path):
            with open(test_file_path, "r") as copied_file:
                read_data = copied_file.read()
                fail_string = "File copy failed: %s " % read_data + \
                "doesn't equal %s (for %s). " % \
                (test_data[test_file], test_file_path) + \
                "Remove is %s and Overwrite is %s." % \
                (str(remove_before_copy), str(should_overwrite))
                # Assertion that passes when:
                #  - The files got copied if we removed the old files in the
                # target directory before the copy.
                #  - The new files did not replace the old files in the
                # target directory if they where not supposed to overwrite.
                #  - The new files did replace the old files in the target
                # directory of we specified that this should happen.
                assert (read_data == test_data[test_file]) == \
                (remove_before_copy or should_overwrite), fail_string


def test__copy_for_storage():
    """Test using the copy function without any specification
    as to how to do it.
    """
    config = {}
    perform_transfer(_copy_for_storage, config)
    perform_transfer(_copy_for_storage, config,
    remove_before_copy=False)


def test__copy_for_storage_scp():
    """Test using the copy function with scp.
    """
    config = {"transfer_protocol": "scp"}
    perform_transfer(_copy_for_storage, config)
    perform_transfer(_copy_for_storage, config,
    remove_before_copy=False)


def test__copy_for_storage_rsync():
    """Test using the copy function with rsync.
    """
    config = {"transfer_protocol": "rsync"}
    perform_transfer(_copy_for_storage, config)
    perform_transfer(_copy_for_storage, config,
    remove_before_copy=False, should_overwrite=True)

# NOTE: rdiff-backup not available
# def test__copy_for_storage_rdiff_backup():
#     """Test using the copy function with rdiff-backup.
#     """
#     config = {"transfer_protocol": "rdiff-backup"}
#     perform_transfer(_copy_for_storage, config)
#     perform_transfer(_copy_for_storage, config,
#     remove_before_copy=False, should_overwrite=True)


def test__copy_from_sequencer():
    """Test using the copy function without any specification
    as to how to do it.
    """
    config = {}
    perform_transfer(_copy_from_sequencer, config)
    perform_transfer(_copy_from_sequencer, config,
    remove_before_copy=False)


def test__copy_from_sequencer_scp():
    """Test using the copy function with scp.
    """
    config = {"transfer_protocol": "scp"}
    perform_transfer(_copy_from_sequencer, config)
    perform_transfer(_copy_from_sequencer, config,
    remove_before_copy=False)


def test__copy_from_sequencer_rsync():
    """Test using the copy function with rsync.
    """
    config = {"transfer_protocol": "rsync"}
    perform_transfer(_copy_from_sequencer, config)
    perform_transfer(_copy_from_sequencer, config,
    remove_before_copy=False, should_overwrite=True)

# NOTE: rdiff-backup not available
# def test__copy_from_sequencer_rdiff_backup():
#     """Test using the copy function with rdiff-backup.
#     """
#     config = {"transfer_protocol": "rdiff-backup"}
#     perform_transfer(_copy_from_sequencer, config)
#     perform_transfer(_copy_from_sequencer, config,
#     remove_before_copy=False, should_overwrite=True)
