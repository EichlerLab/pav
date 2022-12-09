"""
Variant caller functions.
"""


import collections
import intervaltree
import numpy as np
import pandas as pd
import re

import svpoplib
import pavlib


def get_gt(row, hap, map_tree):
    """
    Get variant genotype based on haplotype and mappability.

    :param row: Variant call row (post h1/h2 merge).
    :param hap: Haplotype being called.
    :param map_tree: Tree of mappable regions.

    :return: '1' if the variant call is in this haplotype, '0' if it was not called but is in mappable regions,
        and '.' if it is not in a mappable region.
    """

    if hap in row['HAP'].split(';'):
        return '1'

    for interval in map_tree[row['#CHROM']].overlap(row['POS'], row['END']):

        if (interval.begin <= row['POS']) and (interval.end >= row['END']):
            return '0'

    return '.'


def val_per_hap(df, df_h1, df_h2, col_name, delim=';'):
    """
    Construct a field from a merged variant DataFrame (`df`) by pulling values from each pre-merged haplotype. Matches
    the merged variant IDs from df to the correct original ID in each haplotype and returns a comma-separated string
    of values. `df_h1` and `df_h2` must be indexed with the original variant ID from the corresponding haplotype.
    Function returns a Pandas Series keyed by the `df` index.

    :param df: Merged DataFrame of variants.
    :param df_h1: Pre-merged DataFrame of variants from h1. Index must be variant IDs for the original unmerged calls.
    :param df_h2: Pre-merged DataFrame of variants from h2. Index must be variant IDs for the original unmerged calls.
    :param col_name: Get this column name from `df_h1` and `df_h2`.
    :param delim: Separate values by this delimiter.

    :return: A Pandas series keyed by indices from `df` with comma-separated values extracted from each haplotype
        DataFrame.
    """

    df_dict = {'h1': df_h1, 'h2': df_h2}

    # Generate a Series keyed by variant ID with tuples of (hap, id) for the haplotype and variant ID in that haplotype.
    # Then, apply a function to join the target values from each dataframe in the correct order.
    return df.apply(lambda row:
         tuple(zip(
             row['HAP'].split(';'), row['HAP_VARIANTS'].split(';')
         )),
         axis=1
    ).apply(
        lambda val_list: delim.join(df_dict[val[0]].loc[val[1], col_name] for val in val_list)
        ## Get from original table, correct back to original ID (take off everything after ".')
        #lambda val_list: delim.join(df_dict[val[0]].loc[val[1].rsplit('.', 1)[0], col_name] for val in val_list)
    )


def filter_by_ref_tree(df, filter_tree, match_tig=False):
    """
    Filter DataFrame by a dict (keyed by chromosome) of interval trees.

    :param df: DataFrame to filter.
    :param filter_tree: Dict of trees.
    :param match_tig: Match TIG_REGION contig name to the value in the filter tree for intersected intervals (filter
        tree values should be contig names for the interval). When set, only filters variants inside regions with a
        matching contig. This causes PAV to treat each contig as a separate haplotype.

    :return: Filtered DataFrame.
    """

    if df.shape[0] == 0:
        return df, df.copy()

    if match_tig:
        filter_pass = df.apply(
            lambda row: not np.any(
                [row['TIG_REGION'].split(':', 1)[0] == region.data for region in filter_tree[row['#CHROM']][row['POS']:row['END']]]
            ),
            axis=1
        )


    else:
        filter_pass = df.apply(
            lambda row: len(filter_tree[row['#CHROM']][row['POS']:row['END']]) == 0,
            axis=1
        )

    return df.loc[filter_pass], df.loc[~ filter_pass]


def filter_by_tig_tree(df, tig_filter_tree):
    """
    Filter records from a callset DataFrame by matching "TIG_REGION" with regions in an IntervalTree.

    :param df: DataFrame to filter. Must contain field "TIG_REGION".
    :param tig_filter_tree: A `collections.defaultdict` of `intervaltree.IntervalTree` (indexed by contig name) of
        no-call regions. Variants with a tig region intersecting these records will be removed (any intersect). If
        `None`, then `df` is not filtered.

    :return: A tuple of dataframes: (passed, filtered).
    """

    if tig_filter_tree is None:
        return df, pd.DataFrame([], columns=df.columns)

    if df.shape[0] == 0:
        return df, df.copy()

    rm_index_set = set()

    for index, row in df.iterrows():
        match_obj = re.match('^([^:]+):(\d+)-(\d+)$', row['TIG_REGION'])

        if match_obj is None:
            raise RuntimeError('Unrecognized TIG_REGION format for record {}: {}'.format(index, row['TIG_REGION']))

        if tig_filter_tree[match_obj[1]][int(match_obj[2]) - 1:int(match_obj[3])]:
            rm_index_set.add(index)

    # Return
    return df.loc[[val not in rm_index_set for val in df.index]].copy(), df.loc[[val in rm_index_set for val in df.index]].copy()


def left_homology(pos_tig, seq_tig, seq_sv):
    """
    Determine the number of perfect-homology bp upstream of an SV/indel using the SV/indel sequence (seq_sv), a contig
    or reference sequence (seq_tig) and the position of the first base upstream of the SV/indel (pos_tig) in 0-based
    coordinates. Both the contig and SV/indel sequence must be in the same orientation (reverse-complement if needed).
    Generally, the SV/indel sequence is in reference orientation and the contig sequence is the reference or an
    aligned contig in reference orientation (reverse-complemented if needed to get to the + strand).

    This function traverses from `pos_tig` to upstream bases in `seq_tig` using bases from the end of `seq_sv` until
    a mismatch between `seq_sv` and `seq_tig` is found. Search will wrap through `seq_sv` if homology is longer than
    the SV/indel.

    WARNING: This function assumes upper-case for the sequences. Differing case will break the homology search. If any
    sequence is None, 0 is returned.

    :param pos_tig: Contig/reference position (0-based) in reference orientation (may have been reverse-complemented by an
        alignment) where the homology search begins.
    :param seq_tig: Contig sequence as an upper-case string and in reference orientation (may have been reverse-
        complemented by the alignment).
    :param seq_sv: SV/indel sequence as an upper-case string.

    :return: Number of perfect-homology bases between `seq_sv` and `seq_tig` immediately upstream of `pos_tig`. If any
        of the sequneces are None, 0 is returned.
    """

    if seq_sv is None or seq_tig is None:
        return 0

    svlen = len(seq_sv)

    hom_len = 0

    while hom_len <= pos_tig:  # Do not shift off the edge of a contig.
        seq_tig_base = seq_tig[pos_tig - hom_len]

        # Do not match ambiguous bases
        if seq_tig_base not in {'A', 'C', 'G', 'T'}:
            break

        # Match the SV sequence (dowstream SV sequence with upstream reference/contig)
        if seq_sv[-((hom_len + 1) % svlen)] != seq_tig_base:
            # Circular index through seq in reverse from last base to the first, then back to the first
            # if it wraps around. If the downstream end of the SV/indel matches the reference upstream of
            # the SV/indel, shift left. For tandem repeats where the SV was placed in the middle of a
            # repeat array, shift through multiple perfect copies (% oplen loops through seq).
            break

        hom_len += 1

    # Return shifted amount
    return hom_len


def right_homology(pos_tig, seq_tig, seq_sv):
    """
    Determine the number of perfect-homology bp downstream of an SV/indel using the SV/indel sequence (seq_sv), a contig
    or reference sequence (seq_tig) and the position of the first base downstream of the SV/indel (pos_tig) in 0-based
    coordinates. Both the contig and SV/indel sequence must be in the same orientation (reverse-complement if needed).
    Generally, the SV/indel sequence is in reference orientation and the contig sequence is the reference or an
    aligned contig in reference orientation (reverse-complemented if needed to get to the + strand).

    This function traverses from `pos_tig` to downstream bases in `seq_tig` using bases from the beginning of `seq_sv` until
    a mismatch between `seq_sv` and `seq_tig` is found. Search will wrap through `seq_sv` if homology is longer than
    the SV/indel.

    WARNING: This function assumes upper-case for the sequences. Differing case will break the homology search. If any
    sequence is None, 0 is returned.

    :param pos_tig: Contig/reference position (0-based) in reference orientation (may have been reverse-complemented by an
        alignment) where the homology search begins.
    :param seq_tig: Contig sequence as an upper-case string and in reference orientation (may have been reverse-
        complemented by the alignment).
    :param seq_sv: SV/indel sequence as an upper-case string.

    :return: Number of perfect-homology bases between `seq_sv` and `seq_tig` immediately downstream of `pos_tig`. If any
        of the sequences are None, 0 is returned.
    """

    if seq_sv is None or seq_tig is None:
        return 0

    svlen = len(seq_sv)
    tig_len = len(seq_tig)

    hom_len = 0
    pos_tig_limit = tig_len - pos_tig

    while hom_len < pos_tig_limit:  # Do not shift off the edge of a contig.
        seq_tig_base = seq_tig[pos_tig + hom_len]

        # Do not match ambiguous bases
        if seq_tig_base not in {'A', 'C', 'G', 'T'}:
            break

        # Match the SV sequence (dowstream SV sequence with upstream reference/contig)
        if seq_sv[hom_len % svlen] != seq_tig_base:
            # Circular index through seq in reverse from last base to the first, then back to the first
            # if it wraps around. If the downstream end of the SV/indel matches the reference upstream of
            # the SV/indel, shift left. For tandem repeats where the SV was placed in the middle of a
            # repeat array, shift through multiple perfect copies (% oplen loops through seq).
            break

        hom_len += 1

    # Return shifted amount
    return hom_len


def merge_haplotypes(h1_file_name, h2_file_name, h1_callable, h2_callable, config_def, threads=1, chrom=None, is_inv=None):
    """
    Merge haplotypes for one variant type.

    :param h1_file_name: h1 variant call BED file name.
    :param h2_file_name: h2 variant call BED file name.
    :param h1_callable: h1 callable region BED file name.
    :param h2_callable: h2 callable region BED file name.
    :param config_def: Merge definition.
    :param threads: Number of threads for each merge.
    :param chrom: Chromosome to merge or `None` to merge all chromosomes in one step.
    :param is_inv: Add inversion columns if `True`, autodetect if `None`.

    :return: A dataframe of variant calls.
    """

    # Merge
    df = svpoplib.svmerge.merge_variants(
        bed_list=[h1_file_name, h2_file_name],
        sample_names=['h1', 'h2'],
        strategy=config_def,
        threads=threads,
        subset_chrom=chrom
    )

    #df.set_index(df['ID'].apply(lambda val: val.rsplit('.', 1)[0]), inplace=True, drop=False)
    df.set_index('ID', inplace=True, drop=False)
    df.index.name = 'INDEX'

    # Set is_inv
    if is_inv is None:
        is_inv = np.any(df['SVTYPE'] == 'INV')

    if is_inv and not np.all(df['SVTYPE'] == 'INV'):
        raise RuntimeError('Detected inversions in merge, but not all variants are inversions ({} of {})'.format(
            np.sum(df['SVTYPE'] == 'INV'), df.shape[0]
        ))

    # Restructure columns
    if 'HAP' in df.columns:
        del (df['HAP'])

    if 'DISC_CLASS' in df.columns:
        del (df['DISC_CLASS'])

    if 'HAP_AC' in df.columns:
        del (df['HAP_AC'])

    if 'HAP_AF' in df.columns:
        del (df['HAP_AF'])

    df.columns = [re.sub('^MERGE_', 'HAP_', val) for val in df.columns]

    del (df['HAP_SRC'])
    del (df['HAP_SRC_ID'])

    df.columns = ['HAP' if val == 'HAP_SAMPLES' else val for val in df.columns]

    if df.shape[0] > 0:
        # Change , to ; from merger
        df['HAP'] = df['HAP'].apply(lambda val: ';'.join(val.split(',')))
        df['HAP_VARIANTS'] = df['HAP_VARIANTS'].apply(lambda val: ';'.join(val.split(',')))

        if 'HAP_RO' in df.columns:
            df['HAP_RO'] = df['HAP_RO'].apply(lambda val: ';'.join(val.split(',')))

        if 'HAP_SZRO' in df.columns:
            df['HAP_SZRO'] = df['HAP_SZRO'].apply(lambda val: ';'.join(val.split(',')))

        if 'HAP_OFFSET' in df.columns:
            df['HAP_OFFSET'] = df['HAP_OFFSET'].apply(lambda val: ';'.join(val.split(',')))

        if 'HAP_OFFSZ' in df.columns:
            df['HAP_OFFSZ'] = df['HAP_OFFSZ'].apply(lambda val: ';'.join(val.split(',')))

        if 'HAP_MATCH' in df.columns:
            df['HAP_MATCH'] = df['HAP_MATCH'].apply(lambda val: ';'.join(val.split(',')))

        # Add h1 and h2 to columns
        df_h1 = svpoplib.pd.read_csv_chrom(h1_file_name, chrom=chrom, sep='\t', low_memory=False)
        df_h1.set_index('ID', inplace=True, drop=False)
        df_h1['CLUSTER_MATCH'].fillna('NA', inplace=True)
        df_h1 = df_h1.astype(str)

        df_h2 = svpoplib.pd.read_csv_chrom(h2_file_name, chrom=chrom, sep='\t', low_memory=False)
        df_h2.set_index('ID', inplace=True, drop=False)
        df_h2['CLUSTER_MATCH'].fillna('NA', inplace=True)
        df_h2 = df_h2.astype(str)

        df['TIG_REGION'] = pavlib.call.val_per_hap(df, df_h1, df_h2, 'TIG_REGION')
        df['QUERY_STRAND'] = pavlib.call.val_per_hap(df, df_h1, df_h2, 'QUERY_STRAND')
        df['CI'] = pavlib.call.val_per_hap(df, df_h1, df_h2, 'CI')
        df['ALIGN_INDEX'] = pavlib.call.val_per_hap(df, df_h1, df_h2, 'ALIGN_INDEX')
        df['CLUSTER_MATCH'] = pavlib.call.val_per_hap(df, df_h1, df_h2, 'CLUSTER_MATCH')
        df['CALL_SOURCE'] = pavlib.call.val_per_hap(df, df_h1, df_h2, 'CALL_SOURCE')

        # Set inversion columns
        if is_inv:
            del(df['RGN_REF_DISC'])
            del(df['RGN_TIG_DISC'])
            del(df['FLAG_ID'])
            del(df['FLAG_TYPE'])

            df['RGN_REF_INNER'] = pavlib.call.val_per_hap(df, df_h1, df_h2, 'RGN_REF_INNER')
            df['RGN_TIG_INNER'] = pavlib.call.val_per_hap(df, df_h1, df_h2, 'RGN_TIG_INNER')

        # Load mapped regions
        map_tree_h1 = collections.defaultdict(intervaltree.IntervalTree)
        map_tree_h2 = collections.defaultdict(intervaltree.IntervalTree)

        df_map_h1 = pd.read_csv(h1_callable, sep='\t')
        df_map_h2 = pd.read_csv(h2_callable, sep='\t')

        for index, row in df_map_h1.iterrows():
            map_tree_h1[row['#CHROM']][row['POS']:row['END']] = True

        for index, row in df_map_h2.iterrows():
            map_tree_h2[row['#CHROM']][row['POS']:row['END']] = True

        # Get genotypes setting no-call for non-mappable regions
        df['GT_H1'] = df.apply(pavlib.call.get_gt, hap='h1', map_tree=map_tree_h1, axis=1)
        df['GT_H2'] = df.apply(pavlib.call.get_gt, hap='h2', map_tree=map_tree_h2, axis=1)

        df['GT'] = df.apply(lambda row: '{}|{}'.format(row['GT_H1'], row['GT_H2']), axis=1)

        if np.any(df['GT'].apply(lambda val: val == '0|0')):
            raise RuntimeError('Program bug: Found 0|0 genotypes after merging haplotypes')

        del df['GT_H1']
        del df['GT_H2']

    else:

        df['TIG_REGION'] = np.nan
        df['QUERY_STRAND'] = np.nan
        df['CI'] = np.nan
        df['ALIGN_INDEX'] = np.nan
        df['CLUSTER_MATCH'] = np.nan
        df['CALL_SOURCE'] = np.nan

        if is_inv:
            del(df['RGN_REF_DISC'])
            del(df['RGN_TIG_DISC'])
            del(df['FLAG_ID'])
            del(df['FLAG_TYPE'])

            df['RGN_REF_INNER'] = np.nan
            df['RGN_TIG_INNER'] = np.nan

        df['GT'] = np.nan

    # Return merged BED
    return df


def get_merge_params(wildcards, config):
    """
    Get merging parameters.

    :param wildcards: Rule wildcards.
    :param config: Config parameters.

    :return: An SV-Pop merge definition string describing how variants should be intersected.
    """

    # Get svtype
    vartype, svtype = wildcards.vartype_svtype.split('_')

    # Get merge parameters
    config_def = None

    if svtype in {'ins', 'del', 'inv'}:  # Search config for a matching override
        if f'merge_{svtype}' in config:
            config_def = config[f'merge_{svtype}']
        elif 'merge_insdel' in config:
            config_def = config['merge_insdel']
        elif 'merge_insdelinv' in config:
            config_def = config['merge_insdelinv']

    elif svtype == 'snv' and 'merge_snv' in config:
        config_def = config['merge_snv']

    if config_def is None:  # Get default
        config_def = pavlib.constants.MERGE_PARAM_DEFAULT.get(svtype, None)

    if config_def is None:
        raise RuntimeError(f'No merge parameters for svtype: {svtype}')

    # Return config definition
    return config_def
