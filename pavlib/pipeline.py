"""
PAV pipeline utilities.

Functions for finding data.
"""

import collections
import itertools
import numpy as np
import os
import pandas as pd
import re

import svpoplib

from Bio import SeqIO
import Bio.bgzf


def expand_pattern(pattern, asm_table, **kwargs):

    # Check kwargs
    if kwargs is None or len(kwargs) == 0:
        # Empty kwargs so the kwargs product iterates over one item
        kwargs = {None: (None,)}

    for key, val in kwargs.items():
        if isinstance(val, str):
            # Single string value, iterate over one item, not each string character
            kwargs[key] = (val,)

    kwargs_keys = sorted(kwargs.keys())
    kwargs_n = len(kwargs_keys)

    sub_dict = dict()  # String substitution dict

    if 'asm_name' in kwargs.keys() and kwargs['asm_name'] is not None:
        asm_list = set(kwargs['asm_name'])
    else:
        asm_list = asm_table.index

    if 'hap' in kwargs.keys() and kwargs['hap'] is not None:
        hap_set = set(kwargs['hap'])
    else:
        hap_set = None

    # Process each assembly
    for asm_name in asm_list:
        sub_dict['asm_name'] = asm_name

        for hap in get_hap_list(asm_name, asm_table):

            if hap_set is not None and hap not in hap_set:
                continue

            sub_dict['hap'] = hap

            for kw_prod in itertools.product(*[kwargs[key] for key in kwargs_keys]):
                for i in range(kwargs_n):
                    sub_dict[kwargs_keys[i]] = kw_prod[i]

                yield pattern.format(**sub_dict)


def get_hap_list(asm_name, asm_table):
    """
    Get a list of haplotypes for an assembly name.

    :param asm_name: Assembly name.
    :param asm_table: Assembly table (assemblies.tsv) as a Pandas DataFrame.
    :param config: Pipeline config dictionary.
    """

    # Check values
    if asm_name is None or asm_name.strip() == '':
        raise RuntimeError('Cannot get assembly config: "asm_name" is missing')

    asm_name = asm_name.strip()

    # Get assembly table entry
    if asm_name not in asm_table.index:
        raise RuntimeError(f'Assembly name not found in the sample config: {asm_name}')

    asm_table_entry = asm_table.loc[asm_name]

    # Find haplotypes for a table entry
    hap_list = [
        col[len('HAP_'):]
            for col in asm_table_entry[~ pd.isnull(asm_table_entry)].index
            if col.startswith('HAP_')
    ]

    if len(hap_list) == 0:
        raise RuntimeError(f'No haplotypes found for assembly {asm_name}: All hap columns are missing or empty')

    return hap_list


def get_asm_config(asm_name, hap, asm_table, config):
    """
    Get a dictionary of parameters and paths for one assembly.

    :param asm_name: Assembly name.
    :param hap: Haplotype name (e.g. "h1", "h2").
    :param asm_table: Assembly table.
    :param config: Pipeline config dictionary.
    """

    config = config.copy()  # Altered by overridden configuration options

    # Check values
    if hap is None or hap.strip() == '':
        raise RuntimeError('Cannot get assembly config: "hap" is missing')

    if asm_name is None or asm_name.strip() == '':
        raise RuntimeError('Cannot get assembly config: "asm_name" is missing')

    asm_name = asm_name.strip()
    hap = hap.strip()

    if asm_name not in asm_table.index:
        raise RuntimeError(f'No assembly table entry: {asm_name}')

    if f'HAP_{hap}' not in asm_table.columns:
        raise RuntimeError(f'No haplotype in assembly table columns: {hap}')

    # Get assembly table entry
    asm_table_entry = asm_table.loc[asm_name]

    # Get config override
    config_string = asm_table_entry['CONFIG']

    if not pd.isnull(config_string):
        config_string = config_string.strip()
    else:
        config_string = ''

    if not config_string:
        config_string = None

    config_override = get_config_override_dict(config_string)

    # Get filename pattern
    assembly_input = asm_table_entry[f'HAP_{hap}']

    if not pd.isnull(assembly_input):
        assembly_input = assembly_input.strip().format(asm_name=asm_name, sample=asm_name, hap=hap)

        if not assembly_input:
            assembly_input = None
    else:
        assembly_input = None

    if assembly_input is not None:
        assembly_input = [val.strip() for val in assembly_input.split(';') if val.strip()]
    else:
        assembly_input = list()

    # Get filter
    filter_col = f'FILTER_{hap}'

    if filter_col in asm_table_entry and not pd.isnull(asm_table_entry[filter_col]):
        filter_input = asm_table_entry[filter_col].strip().format(asm_name=asm_name, sample=asm_name, hap=hap)
        filter_input = [val.strip() for val in filter_input.split(';') if val.strip()]

    else:
        filter_input = list()

    # Return config dictionary
    return {
        'asm_name': asm_name,
        'hap': hap,
        'assembly_input': assembly_input,
        'filter_input': filter_input,
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
    :param config_param: Configuration parameter (returned by get_asm_config()) to get the list of files from.

    :return: A list of input files.
    """

    # Get config
    asm_config = get_asm_config(asm_name, hap, asm_table, config)

    assembly_input = asm_config['assembly_input']

    # Return empty list if there are no paths
    if assembly_input is None or len(assembly_input) == 0:
        return list()

    # Separate semi-colon-delimited paths and fit wildcards
    path_list = list()

    for file_name in assembly_input:

        # Create filename
        if os.stat(file_name).st_size > 0:
            path_list.append(file_name)

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

        elif file_name_ext in {'fasta', 'fa', 'fn', 'fna'}:
            file_name_tuples.append((file_name, 'fasta'))

        elif file_name_ext in {'fastq', 'fq', 'fnq'}:
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
                        raise RuntimeError(f'Duplicate record ID in input: {record.id}')

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


def read_assembly_table(asm_table_filename, config):
    """
    Read assembly table.

    Returns a tuple of 3 elements:
        0) Assembly input table with asm_name (sample) as the key and haplotype names as the columns [pandas.DataFrame].
        1) The configuration column from the input (all NaN/NA if the config column was not present) [pandas.Series].
        3) A map from haplotype names to original column names for generating errors and information related to the
           original input table [dict].

    The table returned by this function renames the input columns (e.g. "HAP1" or "HAP_h1") to the haplotype name (e.g.
    "h1") and keeps them in the order they were found in the input table. All other columns are removed.

    :param asm_table_filename: Input filename to read. If None, produce an empty table.
    :param config: Pipeline configuration.

    :return: A tuple of the haplotype input table, configuration column (Series), and hap to column name map.
    """

    if asm_table_filename is None:
        raise RuntimeError('Cannot read assembly table: None')

    asm_table_filename = asm_table_filename.strip()

    if not os.path.isfile(asm_table_filename):
        raise RuntimeError(f'Assembly table file missing or is not a regular file: {asm_table_filename}')

    if config is None:
        config = dict()

    ignore_cols = set(config.get('ignore_cols', set())) | {'CONFIG'}

    # Read table
    if not os.path.isfile(asm_table_filename):
        raise RuntimeError('Missing assembly table: {}'.format(asm_table_filename))

    filename_lower = asm_table_filename.lower()

    if filename_lower.endswith(('.tsv', '.tsv.gz', '.tsv.txt', 'tsv.txt.gz')):
        df = pd.read_csv(asm_table_filename, sep='\t', header=0, dtype=str)
    elif filename_lower.endswith('.xlsx'):
        df = pd.read_excel(asm_table_filename, header=0, dtype=str)
    elif filename_lower.endswith(('.csv', '.csv.gz', '.csv.txt', '.csv.txt.gz')):
        df = pd.read_csv(asm_table_filename, header=0, dtype=str)
    else:
        raise RuntimeError(f'Unrecoginized table file type (expected ".tsv", ".tsv.gz", ".xlsx", ".csv", ".csv.gz"): {asm_table_filename}')

    # Check for required columns
    if 'NAME' not in df.columns:
        raise RuntimeError('Missing assembly table column: NAME')

    # Check assembly names
    if np.any(pd.isnull(df['NAME'])):
        raise RuntimeError(f'Found {sum(pd.isnull(df["NAME"]))} entries in the assembly table with empty NAME values: {asm_table_filename}')

    bad_name = {name for name in df['NAME'] if pd.isnull(name) or re.search(r'^[a-zA-Z0-9_-]+$', name) is None}

    if bad_name:
        bad_name_str = '"' + '", "'.join(sorted(set(bad_name))[:3]) + '"' + ('...' if len(bad_name) > 3 else '')
        raise RuntimeError(f'Found {len(bad_name)} column names with non-alphanumeric/underscore/dash characters in the name: {bad_name_str}: {asm_table_filename}')

    dup_asm_list = [assembly_name for assembly_name, count in collections.Counter(df.index).items() if count > 1]

    if dup_asm_list:
        raise RuntimeError(f'Found {len(dup_asm_list)} duplicate assembly names "{", ".join(dup_asm_list)}": {asm_table_filename}')


    # Set index and config
    df.set_index('NAME', inplace=True)

    if 'CONFIG' not in df.columns:
        df['CONFIG'] = np.nan

    # Map haplotype names to column names
    hap_list = list()
    hap_col_map = dict()
    filter_list = list()
    unknown_cols = list()

    col_set = set()

    for col in df.columns:

        if col in ignore_cols:
            continue

        if col in col_set:
            raise RuntimeError(f'Duplicate column name "{col}" found in assembly table: {asm_table_filename}')

        match_hap_named = re.search(r'^HAP_([a-zA-Z0-9-+.]+)$', col)
        match_hap_num = re.search(r'^HAP([0-9]+)$', col)
        match_filter = re.search(r'^FILTER_([a-zA-Z0-9-+.]+)$', col)

        if match_hap_named:
            hap = match_hap_named[1]
        elif match_hap_num:
            hap = f'h{match_hap_num[1]}'
        elif match_filter:
            filter_list.append(col)
            continue
        else:
            unknown_cols.append(col)
            continue

        hap_list.append(hap)

        if hap in hap_col_map:
            if col != hap_col_map[hap]:
                dup_source = f'(duplicate column {col})'
            else:
                dup_source = f'(derived from columns {hap_col_map[hap]} and {col})'

            raise RuntimeError(f'Duplicate haplotype name "{hap}" found in assembly table {dup_source}: {asm_table_filename}')

        hap_col_map[hap] = col

    if unknown_cols:
        col_list = ', '.join(unknown_cols[:5]) + '...' if len(unknown_cols) > 5 else ''
        raise RuntimeError(f'Unknown columns in assembly table: {col_list}: {asm_table_filename}')

    # Set index and column names
    df_hap = df[[hap_col_map[hap] for hap in hap_list]]
    df_hap.columns = [f'HAP_{hap}' for hap in hap_list]

    df_filter = df[filter_list]

    filter_hap_map = {
        'FILTER_' + (val[len('HAP_'):] if val.startswith('HAP_') else val): 'FILTER_' + key
            for key, val in hap_col_map.items()
    }

    missing_set = {col for col in df_filter.columns if col not in filter_hap_map}

    if missing_set:
        missing_str = ', '.join(sorted(missing_set)[:3]) + ('...' if len(missing_set) > 3 else '')
        raise RuntimeError(f'Found {len(missing_set)} filter columns that do not match haplotype input columns: {missing_str}: {asm_table_filename}')

    df_filter.columns = [filter_hap_map[col] for col in df_filter.columns]

    return pd.concat([df_hap, df_filter, df[['CONFIG']]], axis=1)


def get_config(wildcards, config, asm_table, key=None, default_val=None, default_none=False):
    """
    Get a config object that might be modified by CONFIG parameters in the assembly table.

    If "key" is None, the full config dictionary is returned. If "key" is defined, then the value of config for
    that key is returned with an optional default value.

    :param wildcards: Rule wildcards.
    :param key: Key of the value to get from config.
    :param default_val: Default value.
    :param default_none: If True and default_val is None, return None instead of throwing an exception.

    :return: Config object. Original global config, if unmodified, or a modified copy of it.
    """

    local_config = get_override_config(config, wildcards.asm_name, asm_table)

    if key is None:
        return local_config

    if key not in local_config:
        if default_val is None and not default_none:
            raise RuntimeError('Configuration does not include key "{}" for asm_name "{}" with no default specified'.format(key, wildcards.asm_name))

        return default_val

    return local_config[key]
