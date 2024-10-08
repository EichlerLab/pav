"""
Manage pipeline configurations
"""

import pandas as pd
import sys
import textwrap

import pavlib

DEFAULT_ALIGNER = 'minimap2'

DEFAULT_ALIGNER_PARAMS = {
    'minimap2': '-x asm5',
    'lra': '',
}

NAMED_ALIGNER_PARAMS = {
    'pav2': {
        'minimap2': '-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5',
        'lra': ''
    }
}

KNOWN_ALIGNERS = ['minimap2', 'lra']


#
# Configuration parameters
#

class ConfigParam(object):
    """
    A configuration parameter object.

    Minimum and maximum values can be either a single value to check or a tuple of (value, inclusive/exclusive). If
    it is a single value, then the min/max is inclusive. If it is a tuple, then min/max is inclusive if the second
    element of the tuple is `True` and exclusive if it is `False`.

    Fields:
    * name: Parameter name.
    * val_type: Type of parameter as a string.
    * default: Default value.
    * min: Minimum value if not `None`.
    * max: Maximum value if not `None`.
    * allowed: Set of allowed values if not `None`.
    * to_lower: String value is converted to lower case if `True`. Only valid if `val_type` is `str`.
    * fail_none: If `True`, fail if a parameter value is `None`, otherwise, return the default value.
    * description: Description of the parameter.
    * is_null: `True` if the parameter is included for documentation purposes, but parameter processing is handled
        outside this class. Alignment parameters must be adjusted for the aligner, and so it is not processed here,
        however, the alignment ConfigParam objects are included to simplify parameter documentation.
    * advanced: `True` if the parameter is an advanced option and should not be shown in brief documentation.
    """

    def __init__(self,
                 name, val_type,
                 default=None,
                 min=None, max=None, allowed=None,
                 to_lower=False, fail_none=False,
                 description=None,
                 is_null=False,
                 advanced=False
                 ):
        self.name = name
        self.val_type = val_type
        self.default = default
        self.min = min
        self.max = max
        self.allowed = allowed
        self.to_lower = to_lower
        self.fail_none = fail_none
        self.description = description
        self.is_null = is_null
        self.advanced = advanced

        # Check name
        if self.name is None or not isinstance(self.name, str) or not self.val_type.strip():
            raise RuntimeError('name is missing or empty')

        self.name = self.name.strip().lower()

        # Check type
        if self.val_type is None or not isinstance(self.val_type, str) or not self.val_type.strip():
            raise RuntimeError('Type is missing or empty')

        self.val_type = self.val_type.strip().lower()

        if self.val_type not in {'int', 'float', 'bool', 'str'}:
            raise RuntimeError(f'Unrecognized parameter type: {self.val_type}')

        # Check remaining parameters
        if allowed is not None and not isinstance(allowed, set):
            raise RuntimeError(f'Allowed values must be a set: {type(allowed)}')

        if to_lower not in {True, False}:
            raise RuntimeError(f'to_lower must be True or False (bool): {type(to_lower)}')

        if fail_none not in {True, False}:
            raise RuntimeError(f'fail_none must be True or False (bool): {type(fail_none)}')

    def get_value(self, val):
        """
        Check and get value.

        :param val: Value to check.

        :return: Value after checking and type conversion.
        """

        # Check default.
        if val is None:
            if self.fail_none:
                raise RuntimeError(f'Missing value for parameter {self.name}: Receieved None')

            val = self.default

        # Check and cast type
        if self.val_type == 'int':
            try:
                val = int(val)
            except ValueError:
                raise RuntimeError(f'Failed casting {self.name} to int: {str(val)}')

        elif self.val_type == 'float':
            try:
                val = float(val)
            except ValueError:
                raise RuntimeError(f'Failed casting {self.name} to float: {str(val)}')
        elif self.val_type == 'bool':
            val = pavlib.util.as_bool(val, fail_to_none=True)

            if val is None:
                raise RuntimeError(f'Failed casting {self.name} to bool: {str(val)}')

        elif self.val_type == 'str':
            val = str(val)

        else:
            raise RuntimeError(f'Unrecognized parameter type (PROGRAM BUG) for {self.name}: {self.val_type}')

        # Convert to lower case
        if self.to_lower:
            if not self.val_type == 'str':
                raise RuntimeError(f'Cannot specify `to_lower=True` for non-string type {self.name}: {self.val_type}')

            val = val.lower()

        # Check allowed values
        if self.allowed is not None and val not in self.allowed:
            raise RuntimeError(f'Illegal value for {self.name}: {val} (allowed values: {self.allowed})')

        # Enforce min/max
        if self.min is not None:
            if isinstance(self.min, tuple):
                min_val, min_inclusive = self.min
            else:
                min_val, min_inclusive = self.min, True

            if val < min_val or (val == min_val and not min_inclusive):
                raise RuntimeError(f'Illegal range for {self.name}: Minimum allowed value is {min_val} ({"inclusive" if min_inclusive else "exclusive"})')

        if self.max is not None:
            if isinstance(self.max, tuple):
                max_val, max_inclusive = self.max
            else:
                max_val, max_inclusive = self.max, True

            if val > max_val or (val == max_val and not max_inclusive):
                raise RuntimeError(f'Illegal range for {self.name}: Maximum allowed value is {max_val} ({"inclusive" if max_inclusive else "exclusive"})')

        # Done converting and checking
        return val


CONFIG_PARAM_LIST = [

    # Alignments
    ConfigParam('aligner', 'str', allowed={'minimap2', 'lra'}, is_null=True,
                description='Alignment program to use.'),
    ConfigParam('align_params', 'str', is_null=True,
                description='Parameters for the aligner. Default depends on aligner (minimap2: "-x asm5", lra: ""). '
                            'Keyword "pav2" reverts to legacy parameters used by PAV versions 1 & 2'),
    ConfigParam('align_score_model', 'str', pavlib.align.score.DEFAULT_ALIGN_SCORE_MODEL,
                description='Default alignment score model as a string argument to pavlib.align.score.get_score_model(). '
                            'These parameters are also used for scoring large variants.',
                advanced=True),
    ConfigParam('redundant_callset', 'bool', False,
                description='Per haplotype assembly, callset is nonredundant per assembled sequence instead of globally '
                            'across all assembly sequences. Allows for multiple representations of the same locus '
                            'assembled in different sequences. May be useful for somatic variation, but requires more '
                            'significant downstream work, but will increase false-positive calls and requires more '
                            'downstream processind and QC to obtain a good-quality callset.',
                advanced=True),

    # Variant calling
    ConfigParam('merge_partitions', 'int', 20, min=1,
                description='Split variants into this many partitions to merge.',
                advanced=True),
    ConfigParam('cigar_partitions', 'int', 10, min=1,
                description='For intra-alignment (not alignment truncating), split chromosomes into this many '
                            'partitions and search for INS/DEL/SNV inside alignment records for each partition',
                advanced=True),
    ConfigParam('query_filter', 'str', None,
                description='Query filter BED file. May be multiple file names separated by semicolons (";"). Each BED '
                            'file contains regions in query-coordinates (assembled sequence) matching sequences in the '
                            'input FASTA file. Any variants intersecting these loci are dropped from the callset. May '
                            'be used to  apply quality filtering for known mis-assembled loci.'),
    ConfigParam('min_anchor_score', 'str', pavlib.const.DEFAULT_MIN_ANCHOR_SCORE,
                description='Minimum score of an aligned segement to allow it to be used as an anchor. This value may '
                            'be the absolute score value or a relative value adjusted for the score of a perfectly '
                            'aligned segment of some length (e.g. "1000bp" would be the score of 1000 aligned bases '
                            'with no gaps or mismatches, i.e. 2000 with default alignment parameters with match=2). '
                            'Any alignment record with a score of at least this value may be used as an anchor for '
                            'alignment-truncating variants.'),
    ConfigParam('lg_dot_graph', 'bool', True,
                description='Write a "dot" file describing the graph structure of each query sequence. The dot file '
                            'be transformed to an image with "graphviz" (https://graphviz.org/). The file can be found '
                            'in "results/{asm_name}/lg_sv/lgsv_graph_{hap}.tar"'),

    # Inversion site flagging from variant call clusters
    ConfigParam('inv_sig_cluster_win', 'int', 200,
                description='Cluster variants within this many bases'),
    ConfigParam('inv_sig_cluster_win_min', 'int', 500,
                description='Window must reach this size'),
    ConfigParam('inv_sig_cluster_snv_min', 'int', 20,
                description='Minimum number if SNVs in window'),
    ConfigParam('inv_sig_cluster_indel_min', 'int', 10,
                description='Minimum number of indels in window'),
    ConfigParam('inv_sig_insdel_cluster_flank', 'int', 2,
                description='For each insertion, multiply SVLEN by this to get the distance to the nearest deletion it may intersect'),
    ConfigParam('inv_sig_insdel_merge_flank', 'int', 2000,
                description='Merge clusters within this distance (bp)'),
    ConfigParam('inv_sig_cluster_svlen_min', 'int', 4,
                description='Discard indels less than this size'),
    ConfigParam('inv_sig_merge_flank', 'int', 500,
                description='Merge windows within this many bp'),
    ConfigParam('inv_sig_part_count', 'int', 10,
                description='Partition signature regions into this many partitions for the caller. Marked here so that this file can be cross-referenced with the inversion caller log'),
    ConfigParam('inv_sig_filter', 'str', 'svindel',
                description='Filter flagged regions'),
    ConfigParam('inv_max_overlap', 'float', 0.2,
                min=0.0, max=1.0,
                description='Maximum allowed reciprocal overlap between inversions in the same haplotype'),

    # Inversions
    ConfigParam('inv_min', 'int', 0, min=0, description='Minimum inversion size'),
    ConfigParam('inv_max', 'int', None, min=0, description='Maximum inversion size'),
    ConfigParam('inv_inner', 'str', 'core',
                allowed={'core', 'none', 'full'}, to_lower=True,
                description='Filter variants inside the inverted core (no flanking repeats) inversions if "core", '
                    'full inversion including repeats if "full", and do not filter if "none"',
                advanced=True),

    ConfigParam('inv_region_limit', 'int', pavlib.const.INV_REGION_LIMIT,
                description='maximum region size when searching for inversions. Value 0 ignores limits and allows regions to be any size',
                advanced=True),
    ConfigParam('inv_min_expand', 'int', pavlib.const.INV_MIN_EXPAND_COUNT,
                description='The default number of region expansions to try (including the initial expansion) and '
                            'finding only fwd k-mer states after smoothing before giving up on the region.',
                advanced=True),
    ConfigParam('inv_init_expand', 'int', pavlib.const.INV_INIT_EXPAND,
                description='Expand the flagged region by this (bp) before starting.',
                advanced=True),
    ConfigParam('inv_min_kmers', 'int', pavlib.const.INV_MIN_KMERS,
                description='Minimum number of k-mers with a distinct state (sum of FWD, FWDREV, and REV). Stop if the '
                            'number of k-mers is less after filtering uninformative and high-count k-mers.',
                advanced=True),
    ConfigParam('inv_max_ref_kmer_count', 'int', pavlib.const.INV_MAX_REF_KMER_COUNT,
                description='If canonical reference k-mers have a higher count than this, they are discarded',
                advanced=True),
    ConfigParam('inv_repeat_match_prop', 'float', pavlib.const.INV_REPEAT_MATCH_PROP,
                description='When scoring INV structures, give a bonus to inverted repeats that are similar in size '
                            'scaled by this factor',
                advanced=True),
    ConfigParam('inv_min_inv_kmer_run', 'int', pavlib.const.INV_MIN_INV_KMER_RUN,
                description='Minimum continuous run of strictly inverted k-mers',
                advanced=True),
    ConfigParam('inv_min_qry_ref_prop', 'float', pavlib.const.INV_MIN_QRY_REF_PROP,
                description='Minimum query and reference region size proportion',
                advanced=True),

    ConfigParam('inv_k_size', 'int', pavlib.const.INV_K_SIZE, description='K-mer size', advanced=True),

    ConfigParam('inv_kde_bandwidth', 'float', pavlib.const.INV_KDE_BANDWIDTH,
                description='Convolution KDE bandwidth',
                advanced=True),
    ConfigParam('inv_kde_trunc_z', 'float', pavlib.const.INV_KDE_TRUNC_Z,
                description='Convolution KDE truncated normal Z-score based on a standard normal (N(0,1)) distribution',
                advanced=True),
    ConfigParam('inv_kde_func', 'str', pavlib.const.INV_KDE_FUNC, allowed={'auto', 'fft', 'conv'}, to_lower=True,
                description='Convolution method. "fft" uses a Fast-Fourier Transform, "conv" is a standard truncated '
                            'normal distribution. "auto" defaults to "fft" if scipy.signal is available and "conv" '
                            'otherwise.',
                advanced=True)
]

CONFIG_PARAM_DICT = {
    param.name: param for param in CONFIG_PARAM_LIST
}


def format_config_md(out_file=sys.stdout, width=80, advanced=True):
    """
    Write markdown-formatted help for configuration options.

    :param out_file: Output file.
    :param width: Line-wrap length.
    :param advanced: Include advanced options.
    """

    for param in CONFIG_PARAM_LIST:

        if not advanced and param.advanced:
            continue

        first_line = f'* {param.name} [{param.val_type}'

        if param.default is not None:
            if param.val_type == 'str':
                first_line += f', "{param.default}"'
            else:
                first_line += f', {param.default}'

        if param.min is not None or param.max is not None:
            if param.min is not None:
                range = ('[' if isinstance(param.min, tuple) and param.min[1] else '(') + str(param.min) + ':'
            else:
                range = '(-inf : '

            if param.max is not None:
                range += str(param.max) + (']' if isinstance(param.max, tuple) and param.max[1] else ')')
            else:
                range += 'inf)'

            first_line += f', {range}'

        if param.allowed is not None:
            first_line += f', {param.allowed}'

        first_line += ']: '

        out_file.write(
            '\n'.join(textwrap.wrap(param.description, initial_indent=first_line, subsequent_indent='  ', width=width))
        )

        out_file.write('\n')


#
# Configuration retrieval
#

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
            raise RuntimeError('The reference configuration parameter cannot be overridden or defined per sample.')

        config[key] = val

    return config


def get_override_config(asm_name, config,  asm_table):
    """
    Get a config dict with values replaced by overridden values. The dict in parameter `config` is copied if it is
    modified. The original (unmodified) config or a modified copy is returned.

    :param asm_name: Name of the assembly.
    :param config: Existing config. Original object will not be modified.
    :param asm_table: Assembly table (DataFrame).

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


def get_config(key, wildcards, config, asm_table):
    """
    Get a config object that might be modified by CONFIG parameters in the assembly table.

    If "key" is None, the full config dictionary is returned. If "key" is defined, then the value of config for
    that key is returned with an optional default value.

    :param key: Key of the value to get from config.
    :param wildcards: Rule wildcards.
    :param config: Global pipeline config (dict).
    :param asm_table: Assembly table (DataFrame).

    :return: Validated value of a configuration parameter cast to the correct type for that value.
    """

    if 'asm_name' not in wildcards.keys():
        local_config = get_override_config(wildcards.asm_name, config, asm_table)
    else:
        local_config = config

    if key is None:
        raise RuntimeError(f'key is missing')

    if key not in CONFIG_PARAM_DICT:
        raise RuntimeError(f'Unknown parameter "{key}" for assembly "{wildcards.asm_name}"')

    return CONFIG_PARAM_DICT[key].get_value(local_config[key] if key in local_config else None)


#
# Aligner parameters
#

def get_aligner(wildcards, config, asm_table):
    """
    Get aligner.

    :param wildcards: Rule wildcards. If this is a string object, it is assumed to be the assembly name.
    :param config: Pipeline config.
    :param asm_table: Assembly table.

    :return: Name of aligner. Will pull from default values if not overridden by config.
    """

    if isinstance(wildcards, str):

        asm_name = wildcards.strip()

        if not asm_name:
            raise RuntimeError(f'Cannot get aligner for assembly: wildcards is and empty string')

    else:

        if 'asm_name' not in wildcards.keys():
            raise RuntimeError(f'Cannot get aligner for assembly: wildcards.asm_name is missing')

        asm_name = wildcards.asm_name

    config = get_override_config(asm_name, config, asm_table)

    # Get aligner
    if 'aligner' in config:
        aligner = config.get('aligner').strip()

        if aligner not in KNOWN_ALIGNERS:
            raise RuntimeError(f'Unknown aligner in parameter "aligner": {aligner}: Known aligners are {", ".join(KNOWN_ALIGNERS)}')

    else:
        aligner = DEFAULT_ALIGNER

    return aligner


def get_aligner_params(wildcards, config, asm_table, aligner=None):
    """
    Get alignment parameters.

    :param wildcards: Rule wildcards. If this is a string object, it is assumed to be the assembly name.
    :param config: Pipeline config.
    :param asm_table: Assembly table.
    :param aligner: Alignment program name (minimap2, lra).

    :return: A string of parameters for the aligner. Will pull from default values if not overridden by config.
    """

    # Check arguments
    if isinstance(wildcards, str):

        asm_name = wildcards.strip()

        if not asm_name:
            raise RuntimeError(f'Cannot get aligner for assembly: wildcards is and empty string')

    else:

        if 'asm_name' not in wildcards.keys():
            raise RuntimeError(f'Cannot get aligner for assembly: wildcards.asm_name is missing')

        asm_name = wildcards.asm_name

    config = get_override_config(asm_name, config, asm_table)

    if aligner is None:
        aligner = get_aligner(wildcards, config, asm_table)

    # Get aligner config parameter
    align_param_key = None

    if f'{aligner}_params' in config:
        if 'align_params' in config and str(config[f'{aligner}_params']).strip() != str(config['align_params']).strip():
            raise RuntimeError(f'Conflicting alignment parameters: Both "{aligner}_params", and "align_params" are defined and are not the same parameter string.')

        align_param_key = f'{aligner}_params'

    elif 'align_params' in config:
        align_param_key = 'align_params'

    if align_param_key is None:
        return DEFAULT_ALIGNER_PARAMS.get(aligner, None)

    align_params = config[align_param_key].strip()

    if align_params.lower() in NAMED_ALIGNER_PARAMS:
        if aligner not in NAMED_ALIGNER_PARAMS[align_params.lower()]:
            raise RuntimeError(f'Named alignment parameters are not defined for this aligner: {align_params}, aligner={aligner}')

        return NAMED_ALIGNER_PARAMS[align_params.lower()][aligner]

    return align_params


def get_aligner_input(wildcards, config, asm_table, aligner=None):
    """

    :param wildcards: Rule wildcards.
    :param config: Pipeline config.
    :param asm_table: Assembly table.
    :param aligner: Alignment program name (minimap2, lra).

    :return: A list of input files for this aligner and sample.
    """

    # Check parameters
    if 'asm_name' not in wildcards.keys():
        raise RuntimeError(f'Cannot get alignment parameters for assembly: wildcards.asm_name is missing')

    config = get_override_config(wildcards.asm_name, config, asm_table)

    if aligner is None:
        aligner = get_aligner(wildcards, config, asm_table)

    # Return list of input file (FASTA file is first)
    if aligner == 'minimap2':
        return [
            'data/query/{asm_name}/query_{hap}.fa.gz',
            'data/query/{asm_name}/query_{hap}.fa.gz.fai',
            'data/query/{asm_name}/query_{hap}.fa.gz.gzi'
        ]

    if aligner == 'lra':
        return [
            'temp/{asm_name}/align/query/query_{hap}.fa',
            'data/ref/ref.fa.gz.gli',
            'data/ref/ref.fa.gz.mms'
        ]

    raise RuntimeError(f'Unknown aligner: {aligner}')

