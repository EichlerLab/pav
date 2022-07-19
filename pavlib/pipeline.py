"""
PAV pipeline utilities.

Functions for finding data.
"""

import numpy as np
import os
import pandas as pd
import re

import svpoplib

from Bio import SeqIO
import Bio.bgzf



def get_asm_config(asm_name, hap, asm_table, config):
    """
    Get a dictionary of parameters and paths for one assembly.

    :param asm_name: Assembly name.
    :param hap: Haplotype (h1, h2).
    :param asm_table: Assembly table (assemblies.tsv) as a Pandas DataFrame.
    :param config: Pipeline config dictionary.
    """

    config = config.copy()  # Altered by overridden configuration options

    # Check values
    if hap is None or hap.strip() == '':
        raise RuntimeError('Cannot get assembly config: "hap" is missing')

    if not re.match('^h\d+$', hap):
        raise RuntimeError('Cannot get assembly config: "hap" is malformed: Expected "h" followed by an integer: {hap}'.format(hap=hap))

    if asm_name is None or asm_name.strip() == '':
        raise RuntimeError('Cannot get assembly config: "asm_name" is missing')

    asm_name = asm_name.strip()
    hap = hap.strip()

    # Get assembly table entry
    if asm_name in asm_table.index:
        asm_source = 'table'
        asm_table_entry = asm_table.loc[asm_name]
    else:
        asm_source = 'config'
        asm_table_entry = None

    # Get config override
    config_string = None

    if asm_table_entry is not None and 'CONFIG' in asm_table_entry:
        config_string = asm_table_entry['CONFIG']

        if not pd.isnull(config_string) and config_string is not None:
            config_string = config_string.strip()

            if not config_string:
                config_string = None
        else:
            config_string = None

    config_override = get_config_override_dict(config_string)
    config = get_config_with_override(config, config_override)

    # Get input filename pattern (table takes precedence over config)
    if asm_source == 'table':

        # Get haplotype
        haplotype_col = re.sub('^h', 'HAP', hap)

        if haplotype_col not in asm_table_entry:
            raise RuntimeError(
                'Cannot get assembly config: No haplotype column "{}" in assembly table for haplotype "{}"'.format(
                    haplotype_col, hap
                )
            )

        filename_pattern = asm_table_entry[haplotype_col]

        if not pd.isnull(filename_pattern):
            filename_pattern = filename_pattern.strip()

            if not filename_pattern:
                filename_pattern = None
        else:
            filename_pattern = None

    else:
        # Get from config pattern
        filename_pattern = config.get('asm_pattern', None)

        if filename_pattern is None:
            raise RuntimeError('Cannot get assembly config: {}: No assembly table entry and no "asm_pattern" in config and no matching assembly table entry'.format(asm_name))

        if '{asm_name}' not in filename_pattern:
            raise RuntimeError('Cannot get contig filter BED: {}: Missing {{asm_name}} wildcard in config["asm_pattern"]: {}'.format(asm_name, filename_pattern))

        if '{hap}' not in filename_pattern:
            raise RuntimeError('Cannot get contig filter BED: {}: Missing {{hap}} wildcard in config["asm_pattern"]: {}'.format(asm_name, filename_pattern))

    # Get filter
    tig_filter_pattern = None
    tig_filter_source = None

    filter_col = re.sub('^h', 'FILTER_HAP', hap)

    if asm_table_entry is not None and filter_col in asm_table_entry:
        tig_filter_pattern = asm_table_entry[filter_col]
        tig_filter_source = 'table'

    elif 'tig_filter_pattern' in config:
        tig_filter_pattern = config['tig_filter_pattern']
        tig_filter_source = 'config'

        if not pd.isnull(tig_filter_pattern) and tig_filter_pattern is not None:

            if '{asm_name}' not in tig_filter_pattern:
                raise RuntimeError('Cannot get contig filter BED: {}: Missing {{asm_name}} wildcard in config["tig_filter_pattern"]: {}'.format(asm_name, tig_filter_pattern))

            if '{hap}' not in tig_filter_pattern:
                raise RuntimeError('Cannot get contig filter BED: {}: Missing {{hap}} wildcard in config["tig_filter_pattern"]: {}'.format(asm_name, tig_filter_pattern))

        else:
            tig_filter_pattern = None

    if pd.isnull(tig_filter_pattern):
        tig_filter_pattern = None

    if tig_filter_pattern is not None:
        tig_filter_pattern = tig_filter_pattern.strip()

        if not tig_filter_pattern:
            tig_filter_pattern = None

    if tig_filter_pattern is not None:
        try:
            tig_filter_bed = tig_filter_pattern.format(asm_name=asm_name, hap=hap)
        except KeyError as e:
            raise RuntimeError('Cannot get contig filter BED: {}: Unknown wildcard in contig filter filename pattern: {}: {}'.format(asm_name, e, tig_filter_pattern))
    else:
        tig_filter_bed = None

    # Return config dictionary
    return {
        'asm_name': asm_name,
        'hap': hap,
        'asm_source': asm_source,
        'filename_pattern': filename_pattern,
        'tig_filter_bed': tig_filter_bed,
        'tig_filter_source': tig_filter_source,
        'config_override_string': config_string,
        'config_override': config_override
    }


def get_asm_input_list(asm_name, hap, asm_table, config):
    """
    Get a list of input files.

    :param asm_name: Assembly name.
    :param hap: Haplotype.
    :param asm_table: Assembly table (assemblies.tsv).
    :param config: Pipeline config.

    :return: A list of input files.
    """

    # Get config
    asm_config = get_asm_config(asm_name, hap, asm_table, config)
    config = get_config_with_override(config, asm_config['config_override'])

    record_source = asm_config['asm_source']

    filename_pattern = asm_config['filename_pattern']

    # Return empty list if there are no paths
    if filename_pattern is None:
        return list()

    # Separate semi-colon-delimited paths and fit wildcards
    path_list = list()

    for file_name in filename_pattern.split(';'):

        file_name = file_name.strip()

        if not file_name:
            continue

        # Check for "parent" wildcard (hifiasm-trio)
        if '{parent}' in file_name:

            if hap == 'h1':
                parent = 'mat'
            elif hap == 'h2':
                parent = 'pat'
            else:
                raise RuntimeError('"Cannot get input file for {asm_name} with "{{parent}}" pattern: Haplotype must be "h1" or "h2": Found {hap}: {file_name}'.format(asm_name=wildcards.asm_name, file_name=file_name, hap=wildcards.hap))

        else:
            parent = None

        # Check for required fields if the record comes from a config entry (must have asm_name and hap or parent)
        if record_source == 'config':

            if '{asm_name}' not in file_name:
                raise RuntimeError('Cannot get input file for {asm_name}: Path from config "asm_name" must have {{asm_name}} wildcard: {file_name}'.format(asm_name=wildcards.asm_name, file_name=file_name))

            if '{hap}' not in file_name and parent is None:
                raise RuntimeError('Cannot get input file for {asm_name}: Path from config "asm_pattern" must have {{hap}} or {{parent}} wildcard: {file_name}'.format(asm_name=wildcards.asm_name, file_name=file_name))

        if '{sample}' in file_name:
            delim = config.get('sample_delimiter', '_')

            if delim == '.':
                delim = '\.'

            elif delim not in {'_', '-', '+', '#'}:
                raise RuntimeError('Sample delimiter in config ("sample_delimiter") must be one of "_", ".", "-", "+", "#": {}'.format(delim))

            re_match = re.match('^([^{delim}]+){delim}.*'.format(delim=delim), asm_name)

            if re_match is None:
                raise RuntimeError('"Cannot get input fasta for {asm_name}: "asm_pattern" contains wildcard {{sample}}, but sample cannot be extracted from the assembly name (expected to be at the start of the assembly name and separated by an underscore): {file_name}'.format(asm_name=wildcards.asm_name, file_name=file_name))

            sample = re_match[1]
        else:
            sample = None

        # Create filename
        file_name_parsed = file_name.format(
            sample=sample,
            asm_name=asm_name,
            hap=hap,
            parent=parent
        )

        if not os.path.isfile(file_name_parsed):
            # Try substituting "hap1" for "h1" and "hap2" for "h2" if h1/h2 is not found.
            # Return this version of the file if it exists, let the pipeline fail on the h1/h2 filename if neither exists.

            if re.match('^h\d+$', hap):
                alt_hap = re.sub('^h', 'hap', hap)

                file_name_parsed_alt = file_name.format(
                    sample=sample,
                    asm_name=asm_name,
                    hap=alt_hap,
                    parent=parent
                )

                if os.path.isfile(file_name_parsed_alt):
                    file_name_parsed = file_name_parsed_alt

        if file_name_parsed is not None:
            path_list.append(file_name_parsed)

        # Return paths
        return path_list


def expand_input(file_name_list):
    """
    Expand input to a list of tuples:
    * [0]: File name.
    * [1]: File type ("fasta", "fastq", or "gfa")

    File type does not change if the file is gzipped, the downstream input function will resolve that.

    This function traverses FOFN files (file of file names) recursively until FASTA, FASTQ, or GFA files are found.

    :param file_name_list: List of input file names.
    """

    # Check arguments
    if file_name_list is None:
        raise RuntimeError('Cannot create input FASTA: Input name list is None')

    # Check files
    if issubclass(file_name_list.__class__, str):
        file_name_list = [file_name_list]

    elif issubclass(file_name_list.__class__, tuple):
        file_name_list = list(file_name_list)

    elif issubclass(file_name_list.__class__, set):
        file_name_list = sorted(file_name_list)

    elif not issubclass(file_name_list.__class__, list):
        raise RuntimeError(f'Unrecognized type for input file name list: {file_name_list.__class__}')


    # Generate a list of files traversing into FOFN files
    file_name_tuples = list()
    fofn_list = list()  # Set of visited FOFN files (prevent recursive traversal)

    while len(file_name_list) > 0:

        # Next file name
        file_name = file_name_list[0].strip()
        file_name_list = file_name_list[1:]

        if not file_name:
            continue

        # Get file extension
        file_name_lower = file_name.lower()

        if file_name_lower.endswith('.gz'):  # Strip GZ, the downstream input functions will detect file type
            file_name_lower = file_name_lower.rsplit('.', 1)[0]

        if '.' not in file_name_lower:
            raise RuntimeError(f'No recognizable extension in file name: {file_name}')

        file_name_ext = file_name_lower.rsplit('.', 1)[1]

        # Expand FOFN files
        if file_name_ext == 'fofn':

            # Check for recursive FOFN traversal
            file_name_real = os.path.realpath(file_name)

            if file_name_real in fofn_list:
                raise RuntimeWarning(f'Detected recursive FOFN traversal, ignoring redundant entry: {file_name}')
                continue

            fofn_list.append(file_name_real)

            # Append FOFN entries to the input file list
            with svpoplib.seq.PlainOrGzReader(file_name) as in_file:
                for line in in_file:
                    line = line.strip()

                    if line:
                        file_name_list.append(line)

        elif file_name_ext in {'fasta', 'fa', 'fn'}:
            file_name_tuples.append((file_name, 'fasta'))

        elif file_name_ext in {'fastq', 'fq'}:
            file_name_tuples.append((file_name, 'fastq'))

        elif file_name_ext == 'gfa':
            file_name_tuples.append((file_name, 'gfa'))

        else:
            raise RuntimeError(f'Unrecognized file extension {file_name_ext}: {file_name}')

    # Return tuples
    return file_name_tuples, fofn_list


def get_rule_input_list(asm_name, hap, asm_table, config):
    """
    Get a full list of input files.

    :param asm_name: Assembly name.
    :param hap: Haplotype.
    :param asm_table: Assembly table (assemblies.tsv).
    :param config: Pipeline config.

    :return: A list of all files that may affect input. This includes both sequence data files and FOFN files. This
        rule is designed to be called by a Snakemake input rule.
    """

    input_list = get_asm_input_list(asm_name, hap, asm_table, config)

    file_name_tuples, fofn_list = expand_input(input_list)

    return [
        filename for filename in fofn_list
    ] + [
        filename for filename, filetype in file_name_tuples
    ]


def input_tuples_to_fasta(file_name_tuples, out_file_name):
    """
    Convert a list of input files to a single FASTA entry. Input files may be FASTA, FASTQ, or GFA.

    :param file_name_tuples: List of tuples for each input entry ([0]: File name, [1]: File format). The file format
        must be "fasta", "fastq", or "gfa" (case sensitive).
    :param out_file_name: Output file. If `file_name_tuples` or contains only empty files, the output file is also
        empty (0 bytes) indicating to PAV that data for this haplotype are missing.
    """

    # Check input files, fail early
    has_data = False

    if file_name_tuples is None:
        file_name_tuples = []

    for index in range(len(file_name_tuples)):
        file_name, file_format = file_name_tuples[index]

        if file_format not in {'fasta', 'fastq', 'gfa'}:
            raise RuntimeError(f'Unrecognized file format "{file_format}": {file_name}')

        if not os.path.isfile(file_name):
            raise RuntimeError(f'Input file does not exist or is not a regular file: {file_name}')

        if os.stat(file_name).st_size > 0:
            has_data = True
        else:
            file_name_tuples[index] = (file_name_tuples[index][0], 'empty')  # Mark file as empty

    # Stop if there are no records. An empty file to signal downstream steps that there is no data for this haplotype.
    if not has_data:
        with open(out_file_name, 'w') as out_file:
            pass

        return

    # Record iterator
    def input_record_iter():

        record_id_set = set()

        for file_name, file_format in file_name_tuples:

            if file_format in {'fasta', 'fastq'}:
                for record in svpoplib.seq.fa_to_record_iter(file_name, input_format=file_format):

                    if record.id in record_id_set:
                        raise RuntimeError(f'Duplicate record ID in input: {record.id}')

                    record_id_set.add(record.id)

                    yield record

            elif file_format == 'gfa':
                for record in svpoplib.seq.gfa_to_record_iter(file_name):

                    if record.id in record_id_set:
                        raise RuntimeError(f'Duplicate record ID in input: {record_id}')

                    record_id_set.add(record.id)

                    yield record

            elif file_format not in {'skip', 'empty'}:
                raise RuntimeError(f'Program bug: Unrecognized file type "{file_format}" after checking input.')

    # Traverse entries
    with Bio.bgzf.open(out_file_name, 'wb') as out_file:
        Bio.SeqIO.write(input_record_iter(), out_file, 'fasta')

    return


def get_config_override_dict(config_string):
    """
    Get a dictionary of overridden parameters using the CONFIG column of the assembly table.

    :param config_string: Config override string (e.g. attr1=val1;attr2=val2). Must be colon separated and each
        element must have an equal sign. Whitespace around semi-colons and equal signs is ignored.

    :return: Dict of overridden parameters or an empty dict if no parameters were overridden.
    """

    config_override = dict()

    # Check string
    if config_string is None or pd.isnull(config_string) or not config_string.strip():
        return config_override

    config_string = config_string.strip()

    # Process each config directive
    tok_list = config_string.split(';')

    for tok in tok_list:

        # Check tok
        tok = tok.strip()

        if not tok:
            continue

        if '=' not in tok:
            raise RuntimeError('Cannot get assembly config: Missing "=" in CONFIG token {}: {}'.format(tok, config_string))

        # Get attribute and value
        key, val = tok.split('=', 1)

        key = key.strip()
        val = val.strip()

        if not key:
            raise RuntimeError('Cannot get assembly config: Missing key (key=value) in CONFIG token {}: {}'.format(tok, config_string))

        if not val:
            raise RuntimeError('Cannot get assembly config: Missing value (key=value) in CONFIG token {}: {}'.format(tok, config_string))

        # Set
        config_override[key] = val

    return config_override


def get_config_with_override(config, override_config):
    """
    Get a config dict with values replaced by overridden values. The dict in parameter `config` is copied if it is
    modified. The original (unmodified) config or a modified copy is returned.

    :param config: Existing config. Original object will not be modified.
    :param override_config: A defined set of values that will override entries in `config`.

    :return: A config object.
    """

    if override_config is None:
        return config

    if config is None:
        config = dict()

    config = config.copy()

    for key, val in override_config.items():
        if key in {'reference'}:
            raise RuntimeError('The reference configuration parameter cannot be defined per sample.')

        config[key] = val

    return config


def get_override_config(config, asm_name, asm_table):
    """
    Get a config dict with values replaced by overridden values. The dict in parameter `config` is copied if it is
    modified. The original (unmodified) config or a modified copy is returned.

    :param config: Existing config. Original object will not be modified.
    :param asm_name: Name of the assembly.

    :return: A config object.
    """

    if asm_name is None:
        raise RuntimeError('Cannot get config overide for assembly: None')

    if asm_table is None:
        raise RuntimeError('Cannot get config overide for assembly table: None')

    # Get table entry
    if asm_name not in asm_table.index:
        return config

    asm_table_entry = asm_table.loc[asm_name]

    if 'CONFIG' not in asm_table_entry:
        return config

    return get_config_with_override(config, get_config_override_dict(asm_table_entry['CONFIG']))

