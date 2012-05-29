#!/usr/bin/env python
"""Functions for getting barcode statistics from demultiplexing.
"""

import copy
import bcbio.google.connection
import bcbio.google.document
from bcbio.google import spreadsheet as g_spreadsheet
from bcbio.google import _to_unicode
from bcbio.log import logger2
import bcbio.solexa.flowcell
import bcbio.pipeline.flowcell

# The structure of the demultiplex result
BARCODE_STATS_HEADER = [
         ['Project name', 'project_name'],
         ['Lane', 'lane'],
         ['Lane description', 'description'],
         ['Sample name', 'sample_name'],
         ['bcbb internal barcode index', 'bcbb_barcode_id'],
         ['Barcode name', 'barcode_name'],
         ['Barcode sequence', 'barcode_sequence'],
         ['Barcode type', 'barcode_type'],
         ['Demultiplexed read (pair) count', 'read_count'],
         ['Demultiplexed read (pair) count (millions)', 'rounded_read_count'],
         ['Comment', 'comment']
        ]

# The structure of the sequencing result
SEQUENCING_RESULT_HEADER = [
                 ['Sample name', 'sample_name'],
                 ['Run', 'run'],
                 ['Lane', 'lane'],
                 ['Read (pair) count', 'read_count'],
                 ['Read (pair) count (millions)', 'rounded_read_count'],
                 ['Barcode sequence'],
                 ['Comment', 'comment'],
                 ['Pass', 'pass']
                ]


def _create_header(header, columns):
    header = copy.deepcopy(header)
    names = []
    for column in columns:
        for i, head in enumerate(header):
            if (len(head) > 1 and head[1] == column):
                names.append(head[0])
                header.pop(i)
    for head in header:
        names.append(head[0])
    return names


def get_spreadsheet(ssheet_title, encoded_credentials):
    """Connect to Google docs and get a spreadsheet"""

    # Convert the spreadsheet title to unicode
    ssheet_title = _to_unicode(ssheet_title)

    # Create a client class which will make HTTP requests with Google Docs server.
    client = g_spreadsheet.get_client()
    bcbio.google.connection.authenticate(client, encoded_credentials)

    # Locate the spreadsheet
    ssheet = g_spreadsheet.get_spreadsheet(client, ssheet_title)

    # Check that we got a result back
    if not ssheet:
        logger2.warn("No document with specified title '%s' found in \
                      GoogleDocs repository" % ssheet_title)
        return (None, None)

    return (client, ssheet)


def _write_project_report_to_gdocs(client, ssheet, flowcell):
    """Writes report data to corresponding worksheets in a google docs
    spreadsheet.
    """
    # Flatten the project_data structure into a list
    samples = {}
    for sample in flowcell.get_samples():
        if sample.get_name() in samples:
            samples[sample.get_name()].add_sample(sample)
        else:
            samples[sample.get_name()] = sample

    column_headers = [col_header[0] for col_header in SEQUENCING_RESULT_HEADER[:6]]

    success = True
    for sample in samples.values():
        wsheet_title = sample.get_name()
        run_name = "{}_{}".format(flowcell.fc_date, flowcell.fc_name)
        row = (sample.get_name(), \
               run_name, \
               sample.get_lane(), \
               sample.get_read_count(), \
               sample.get_rounded_read_count(), \
               sample.barcode_sequence)

        success &= _write_to_worksheet(client, \
                                       ssheet, \
                                       wsheet_title, \
                                       [row], \
                                       column_headers, \
                                       append=True)

    return success


def write_project_report_summary_to_gdocs(client, ssheet):
    """Summarize the data from the worksheets and write them to a "Summary"
    worksheet.
    """
    # Summary data
    flowcells = {}
    samples = {}
    # Get the list of worksheets in the spreadsheet
    wsheet_feed = g_spreadsheet.get_worksheets_feed(client, ssheet)
    # Loop over the worksheets and parse the data from the ones that contain
    # flowcell data.

    for wsheet in wsheet_feed.entry:
        wsheet_title = wsheet.title.text
        if wsheet_title.endswith("_QC"):
            continue

        if wsheet_title == "Summary":
            continue

        wsheet_data = g_spreadsheet.get_cell_content(client, ssheet, wsheet, '2')
        delim = ';'

        # Add the results from the worksheet to the summarized data.
        for row in wsheet_data:
            sample_name = row[0]
            run_name = row[1]
            lane_name = row[2]
            read_count = row[3]
            try:
                barcode_sequence = row[5]
            except IndexError:
                # For compatibility with old worksheets.
                barcode_sequence = None

            data = {"name": sample_name,
                    "read_count": read_count,
                    "sequence": barcode_sequence}

            lane = bcbio.pipeline.flowcell.Lane({"lane": lane_name})
            sample = bcbio.pipeline.flowcell.BarcodedSample(data, lane)

            logger2.debug("Comment in Sample object: {}".format(sample.get_comment()))

            if sample_name in samples:
                samples[sample_name]["object"].add_sample(sample, delim)
                samples[sample_name]["flowcells"] += "{}{}".format(delim, run_name)
                if not samples[sample_name]["object"].barcode_sequence and barcode_sequence:
                    samples[sample_name]["object"].barcode_sequence = barcode_sequence

            else:
                samples[sample_name] = {"object": sample, "flowcells": run_name}

    wsheet_title = "Summary"

    # Try getting already existing 'comment' and 'pass' values.
    existing_summary_wsheet = g_spreadsheet.get_worksheet(client, ssheet, wsheet_title)
    num_rows = g_spreadsheet.row_count(existing_summary_wsheet)
    name_data = {}
    for row_num in range(2, num_rows + 1):
        content = g_spreadsheet.get_row(client, ssheet, existing_summary_wsheet, row_num)
        sample_name = content[0]
        comment = content[5]
        pass_field = content[6]

        name_data[sample_name] = [comment, pass_field]

    # Flatten the project_data structure into a list
    rows = []
    for sample_data in samples.values():
        sample = sample_data["object"]
        flowcells = sample_data["flowcells"]

        sample_name = sample.get_name()
        comment, pass_field = name_data.get(sample_name, [None, ""])

        logger2.debug("Comment passed in to 'rows': {}".format(comment))

        row = [sample_name, \
               flowcells, \
               sample.get_lane(), \
               sample.get_read_count(), \
               sample.get_rounded_read_count(), \
               sample.barcode_sequence, \
               comment, \
               pass_field]

        rows.append(row)

    # Write the data to the worksheet
    column_headers = [col_header[0] for col_header in SEQUENCING_RESULT_HEADER]
    return _write_to_worksheet(client, \
                               ssheet, \
                               wsheet_title, \
                               rows, \
                               column_headers, \
                               False)


def write_run_report_to_gdocs(fc, fc_date, fc_name, ssheet_title, \
    encoded_credentials, wsheet_title=None, append=False, split_project=False):
    """Upload the barcode read distribution for a run to google docs.
    """

    # Connect to google and get the spreadsheet
    client, ssheet = get_spreadsheet(ssheet_title, encoded_credentials)
    if not client or not ssheet:
        return False

    # Get the projects in the run
    projects = fc.get_project_names()
    logger2.info("Will write data from the run {}_{} for " \
        "projects: {!r}".format(fc_date, fc_name, "', '".join(projects)))

    # If we will split the worksheet by project, use the project
    # names as worksheet titles.
    success = True
    header = _create_header(BARCODE_STATS_HEADER, fc.columns())
    if split_project:
        # Filter away the irrelevent project entries and write the
        # remaining to the appropriate worksheet.
        for project in projects:
            pruned_fc = fc.prune_to_project(project)
            success &= _write_to_worksheet(client, ssheet, project, \
                pruned_fc.to_rows(), header, append)

    # Otherwise, set the default title of the worksheet to be a string of
    # concatenated date and flowcell id.
    else:
        if wsheet_title is None:
            wsheet_title = "{}_{}".format(fc_date, fc_name)

        success &= _write_to_worksheet(client, ssheet, wsheet_title, \
            fc.to_rows(), header, append)

    return success


def _write_to_worksheet(client, ssheet, wsheet_title, rows, header, append):
    """Generic method to write a set of rows to a worksheet on google docs.
    """
    # Convert the worksheet title to unicode
    wsheet_title = _to_unicode(wsheet_title)

    # Add a new worksheet, possibly appending or replacing a pre-existing
    # worksheet according to the append-flag.
    wsheet = g_spreadsheet.add_worksheet(client, \
                                         ssheet, \
                                         wsheet_title, \
                                         len(rows) + 1, \
                                         len(header), \
                                         append)
    if wsheet is None:
        logger2.error("ERROR: Could not add a worksheet {!r} to " \
            "spreadsheet {!r}".format(wsheet_title, ssheet.title.text))
        return False

    # Write the data to the worksheet
    success = g_spreadsheet.write_rows(client, ssheet, wsheet, header, rows)
    if success:
        logger2.info("Wrote data to the {!r}:{!r} " \
                     "worksheet".format(ssheet.title.text, wsheet_title))
    else:
        logger2.error("ERROR: Could not write data to the {!r}:{!r} " \
                      "worksheet".format(ssheet.title.text, wsheet_title))
    return success
