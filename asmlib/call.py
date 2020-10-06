"""
Variant caller functions.
"""

import re

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
    )

def filter_by_tig_tree(df, tig_filter_tree):
    """
    Filter records from a callset DataFrame by matching "TIG_REGION" with regions in an IntervalTree.

    :param df: DataFrame to filter. Must contain field "TIG_REGION".
    :param tig_filter_tree: A `collections.defaultdict` of `intervaltree.IntervalTree` (indexed by contig name) of
        no-call regions. Variants with a tig region intersecting these records will be removed (any intersect). If
        `None`, then `df` is not filtered.

    :return: Filtered `df`.
    """

    if tig_filter_tree is None:
        return df

    rm_index_set = set()

    for index, row in df.iterrows():
        match_obj = re.match('^([^:]+):(\d+)-(\d+)$', row['TIG_REGION'])

        if match_obj is None:
            raise RuntimeError('Unrecognized TIG_REGION format for record {}: {}'.format(index, row['TIG_REGION']))

        if tig_filter_tree[match_obj[1]][int(match_obj[2]) - 1:int(match_obj[3])]:
            rm_index_set.add(index)

    # Return
    if not rm_index_set:
        return df

    return df.loc[[val not in rm_index_set for val in df.index]]

