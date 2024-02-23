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

# Explanations for filter codes
FILTER_REASON = {
    'PASS': 'Variant passed filters',
    'QRY_FILTER': 'Query filter region',
    'COMPOUND': 'Inside larger variant',
    'SVLEN': 'Variant size out of bounds',
    'TRIM': 'Alignment trimming removed variant region'
}


def version_variant_bed_id(df, re_version=False):
    """
    Version IDs in a variant call BED table (pd.Series). Table should have "ID", "FILTER", and "QRY_REGION" fields. If
    ID is missing, it will be generated from the table. Variants are prioritized by FILTER=PASS with "dot" versions
    (e.g. ID.1, ID.2) being preferrentially added to the variant call that do not have FILTER=PASS.

    :param df: Variant call DataFrame.
    :param re_version: If True, eliminate existing versions from IDs before re-versioning.

    :return: A Series object with the new IDs.
    """

    # Determine which columns are available for sorting.
    get_cols = [col for col in ['ID', 'FILTER', 'QRY_REGION'] if col in df.columns]

    if 'ID' not in get_cols:
        id_col = svpoplib.variant.get_variant_id(df, apply_version=False)
    else:
        id_col = None

    # Subset columns
    if len(get_cols) > 0:
        # Make a copy (order is changed)
        df_re = df[get_cols].copy().reset_index()

        if id_col is not None:
            df_re['ID'] = list(id_col)
            get_cols = ['ID'] + get_cols
    else:
        # Create with a single ID column
        df_re = pd.DataFrame(id_col)
        df_re.columns = ['ID']

        get_cols = ['ID']

    if 'FILTER' not in get_cols:
        df_re['FILTER'] = 'PASS'

    if 'QRY_REGION' not in get_cols:
        df_re['QRY_REGION'] = 'chrUn:0-0'

    # Remove versions (.1, .2, etc) if present
    if re_version:
        df_re['ID'] = df_re['ID'].apply(lambda val: val.rsplit('.', 1)[0])

    # Sort by ID, FILTER, then query region
    df_re['FILTER'].fillna('')
    df_re['FILTER'] = df_re['FILTER'].apply(lambda val: ('a' if val == 'PASS' else 'b') + val)  # Force PASS to sort

    qry_region = df_re['QRY_REGION'].apply(pavlib.seq.region_from_string)
    df_re['QRY_ID'] = qry_region.apply(lambda val: val.chrom)
    df_re['QRY_POS'] = qry_region.apply(lambda val: val.pos)

    del(qry_region)

    df_re.sort_values(['ID', 'FILTER', 'QRY_ID', 'QRY_POS'], inplace=True)

    # Re-ID, first PASS, then non-PASS
    df_re.loc[df_re['FILTER'] == 'aPASS', 'ID'] = svpoplib.variant.version_id(
        df_re.loc[df_re['FILTER'] == 'aPASS', 'ID']
    )

    df_re.loc[df_re['FILTER'] != 'aPASS', 'ID'] = svpoplib.variant.version_id(
        df_re.loc[df_re['FILTER'] != 'aPASS', 'ID'],
        existing_id_set=set(df_re.loc[df_re['FILTER'] == 'aPASS', 'ID'])
    )

    # Return IDs
    df_re.sort_index(inplace=True)
    df_re.index = df.index

    return df_re['ID']


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


def val_per_hap(df, df_dict, col_name, delim=';'):
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

    # Generate a Series keyed by variant ID with tuples of (hap, id) for the haplotype and variant ID in that haplotype.
    # Then, apply a function to join the target values from each dataframe in the correct order.
    return df.apply(lambda row:
         tuple(zip(
             row['HAP'].split(';'), row['HAP_VARIANTS'].split(';')
         )),
         axis=1
    ).apply(
        lambda val_list: delim.join(str(df_dict[val[0]].loc[val[1], col_name]) for val in val_list)
    )


def filter_by_ref_tree(row, filter_tree, reason, match_tig=False):
    """
    Filter DataFrame by a dict (keyed by chromosome) of interval trees.

    `filter_tree` is a dictionary with one key per chromosome and an IntervalTree object as the value for each
    chromosome. IntervalTrees have a start, end, and data. Elements in the filter tree are usually constructed from
    the genomic locations of variant calls. The data should be a tuple of two elements: 1) The contig where the variant
    was called from, and 2) The ID of the variant. The first element is needed for `match_tig`, and the second is needed
    to annotate why variants are being dropped if they are inside another.

    :param row: Variant call row.
    :param filter_tree: Dict of trees, one entry per chromosome.
    :param reason: Filter reason. "FILTER" is set to `reason` if it is "PASS", otherwise, it is appended to the
        filter column with a comma. (e.g. "REASON1,REASON2").
    :param match_tig: Match QRY_REGION contig name to the value in the filter tree for intersected intervals (filter
        tree values should be contig names for the interval). When set, only filters variants inside regions with a
        matching contig. This causes PAV to treat each contig as a separate haplotype.

    :return: Filtered DataFrame.
    """

    intersect_set = filter_tree[row['#CHROM']][row['POS']:row['END']]

    if intersect_set and match_tig:
        tig_name = row['QRY_REGION'].split(':', 1)[0]
        intersect_set = {record for record in intersect_set if record.data[0] != tig_name}

    if intersect_set:
        row['FILTER'] = 'COMPOUND'
        row['COMPOUND'] = ','.join(val.data[1] for val in intersect_set)

        if row['FILTER'] == 'PASS':
            return reason

        return row['FILTER'] + f',{reason}'

    return row['REASON']  # Leave unaltered


def filter_by_tig_tree(row, filter_tree, reason):
    """
    Filter records from a callset DataFrame by matching "QRY_REGION" with regions in an IntervalTree.

    :param row: DataFrame to filter. Must contain field "QRY_REGION" and "FILTER" initialized to "PASS" (may contain
        other non-PASS values).
    :param filter_tree: A `collections.defaultdict` of `intervaltree.IntervalTree` (indexed by contig name) of
        no-call regions. Variants with a tig region intersecting these records will be removed (any intersect). If
        `None`, then `df` is not filtered.
    :param reason: Filter reason. "FILTER" is set to `reason` if it is "PASS", otherwise, it is appended to the
        filter column with a comma. (e.g. "REASON1,REASON2").

    :return: Filter column value.
    """

    match_obj = re.match(r'^([^:]+):(\d+)-(\d+)$', row['QRY_REGION'])

    if match_obj is None:
        raise RuntimeError('Unrecognized QRY_REGION format for record {}: {}'.format(row.name, row['QRY_REGION']))

    if filter_tree[match_obj[1]][int(match_obj[2]) - 1:int(match_obj[3])]:
        if row['FILTER'] == 'PASS':
            return reason

        return row['FILTER'] + f',{reason}'

    return row['FILTER']  # Leave unaltered


def read_variant_table(filename_list, drop_dup_id=False, version_id=True):
    """
    Read a variant table and prepare for callset integration. Returns the table and two `collections.defaultdict(set)`
    objects (one for filters, and one for compound variant IDs). Callset integration uses these objects to track
    variant filters (filters) and IDs of variants that cover this variant (compound) (i.e. the ID of a large DEL a
    small variant appears in will appear in the small variant's COMPOUND column and "COMPOUND" will be added to
    filters).

    :param filename_list: List of filenames (list) or a single filename (str) to read.
    :param drop_dup_id: If `True`, drop duplicates by variant ID keeping the first variant in the table. If `False`,
        duplicate IDs are versioned (ID.1, ID.2).
    :param version_id: If `True`, version variants with duplicate IDs. Has no effect if `drop_dup_id` is `True`.

    :return: A tuple of [0] variant table, [1] filter dict, [2] compound dict.
    """

    if isinstance(filename_list, str):
        filename_list = [filename_list]

    df_list = [
        pd.read_csv(
            filename, sep='\t',
            low_memory=False, keep_default_na=False,
            dtype={'#CHROM': str, 'QRY_ID': str}
        )
        for filename in filename_list
    ]

    if len(filename_list) > 1:
        df = pd.concat(df_list, axis=0)
    elif len(filename_list) == 1:
        df = df_list[0]
    else:
        raise RuntimeError('Filename list length is 0 or is not a list')

    df = df.sort_values(
        ['#CHROM', 'POS', 'END', 'ID']
    ).reset_index(drop=True)

    if 'FILTER' not in df.columns:
        df['FILTER'] = 'PASS'

    if drop_dup_id:
        df.drop_duplicates('ID', keep='first', inplace=True)

    if version_id and not drop_dup_id:
        df['ID'] = svpoplib.variant.version_id(df['ID'])

    df.set_index('ID', inplace=True, drop=False)
    df.index.name = 'INDEX'

    filter_dict = collections.defaultdict(set)
    compound_dict = collections.defaultdict(set)

    if 'COMPOUND' in df.columns:
        for index, row in df.loc[~ pd.isnull(df['COMPOUND'])].iterrows():
            compound_dict[index] |= {val.strip() for val in row['COMPOUND'].split(',') if len(val.strip()) > 0}

        del(df['COMPOUND'])

    for index, row in df.iterrows():
        if row['FILTER'] != 'PASS':
            filter_dict[index].add(row['FILTER'])

    return df, filter_dict, compound_dict


class DepthContainer:
    """
    Computes average alignment depth and coverage at SV sites from a coverage BED file.
    """

    def __init__(self, df_cov):
        self.df_cov = df_cov

        if df_cov is None or df_cov.shape[0] == 0:
            raise RuntimeError('Coverage table is missing or empty')

        # Make a table of indexes into df_cov
        chrom = None

        last_end = 0
        first_index = 0

        self.index_dict = dict()

        for index in range(df_cov.shape[0]):
            row = df_cov.iloc[index]

            if row['#CHROM'] != chrom:

                if chrom is not None:
                    self.index_dict[chrom] = (first_index, index)

                if row['#CHROM'] in self.index_dict:
                    raise RuntimeError(f'Discontiguous chromosome order: Found {row["#CHROM"]} in multiple blocks')

                if row['POS'] != 0:
                    raise RuntimeError(f'First record for chromosome {row["#CHROM"]} is not 0: {row["POS"]} (record location {index})')

                chrom = row['#CHROM']
                last_end = row['END']
                first_index = index

            else:
                if row['POS'] != last_end:
                    raise RuntimeError(f'Discontiguous or out of order record in {row["#CHROM"]} (record location {index}): POS={row["POS"]}, expected POS={last_end}')

                last_end = row['END']

        assert chrom is not None, 'Missing chromosome at the end of the coverage table'

        self.index_dict[chrom] = (first_index, df_cov.shape[0])

        # Set state to the first chromosome
        row = self.df_cov.iloc[0]

        self.chrom = row['#CHROM']
        self.index, self.last_index = self.index_dict[self.chrom]

        self.pos = row['POS']
        self.end = row['END']
        self.depth = row['DEPTH']
        self.qry_id = set(row['QRY_ID'].split(',')) if not pd.isnull(row['QRY_ID']) else set()

    def get_depth(self, row):

        # Switch chromosomes
        if row['#CHROM'] != self.chrom:

            if row['#CHROM'] not in self.index_dict:
                sv_id = row['ID'] if 'ID' in row.index else '<UNKNOWN>'
                raise RuntimeError(f'Variant "{sv_id}" (variant row index {row.name}) assigned to chromosome that is not in the depth table: {row["#CHROM"]}')

            self.chrom = row['#CHROM']
            self.index, self.last_index = self.index_dict[self.chrom]

            self.pos = self.df_cov.iloc[self.index]['POS']
            self.end = self.df_cov.iloc[self.index]['END']
            self.depth = self.df_cov.iloc[self.index]['DEPTH']
            self.qry_id = set(self.df_cov.iloc[self.index]['QRY_ID'].split(',')) if not pd.isnull(self.df_cov.iloc[self.index]['QRY_ID']) else set()

        # Catch up to the coverage record where this record begins
        #assert False, 'WORKING ON DEPTH CONTAINER'

        is_end_ins = False

        while row['POS'] >= self.end:
            self.index += 1

            if self.index >= self.last_index:

                # Rescue insertions added to the end of the chromosome. This can happen if the contig aligns up to
                # the end of the reference chromosome without clipping and with extra sequence added as an insertion
                # operation to the end of the alignment (CIGAR string). In this case, the "depth" is the number of
                # query records reaching the end of the reference chromosome that this insertion could have been
                # appended to.

                # self.pos, self.end, self.depth, and self.qry_id are already set

                if not (self.index == self.last_index and row['SVTYPE'] == 'INS' and row['END'] == row['POS'] + 1):
                    # Variant is not in the bounds of this reference chromosome
                    sv_id = row['ID'] if 'ID' in row.index else '<UNKNOWN>'
                    raise RuntimeError(f'Ran out of depth records on "{self.chrom}" to the beginning of variant record {sv_id} (variant row index {row.name})')

                self.index -= 1
                is_end_ins = True
                break

            self.pos = self.df_cov.iloc[self.index]['POS']
            self.end = self.df_cov.iloc[self.index]['END']
            self.depth = self.df_cov.iloc[self.index]['DEPTH']

            self.qry_id = set(self.df_cov.iloc[self.index]['QRY_ID'].split(',')) if not pd.isnull(self.df_cov.iloc[self.index]['QRY_ID']) else set()

        # If the variant fully contained within this coverage record, return stats from this region
        if row['END'] < self.end or is_end_ins:
            return self.depth, 1 if self.depth > 0 else 0, ','.join(sorted(self.qry_id))

        # Get coverage from the variant position to the end of this coverage record
        # sum_depth = self.depth * (self.end - row['POS'])
        # sum_align = (1 if self.depth > 0 else 0) * (self.end - row['POS'])
        # qry_id = set(self.qry_id.split(',')) if not pd.isnull(self.qry_id) else set()

        sum_depth = 0
        sum_align = 0
        qry_id = set()

        step_index = self.index

        svlen = 0

        last_end = self.df_cov.iloc[step_index]['POS']

        # Get coverage from all coverage records fully contained within the variant record
        while row['END'] > last_end:
            if step_index >= self.last_index:
                sv_id = row['ID'] if 'ID' in row.index else '<UNKNOWN>'
                raise RuntimeError(f'Ran out of depth records on "{self.chrom}" to the end of variant record {sv_id} (variant row index {row.name})')

            record_len = min(
                [row['END'], self.df_cov.iloc[step_index]['END']]
            ) - max(
                [row['POS'], self.df_cov.iloc[step_index]['POS']]
            )

            last_end = self.df_cov.iloc[step_index]['END']

            svlen += record_len

            sum_depth += self.df_cov.iloc[step_index]['DEPTH'] * record_len

            sum_align += record_len if self.df_cov.iloc[step_index]['DEPTH'] > 0 else 0

            if not pd.isnull(self.df_cov.iloc[step_index]['QRY_ID']):
                qry_id |= set(self.df_cov.iloc[step_index]['QRY_ID'].split(',')) if not pd.isnull(self.df_cov.iloc[step_index]['QRY_ID']) else {}

            step_index += 1

        assert svlen == row['END'] - row['POS'], f'Record length for "{row["ID"] if "ID" in row.index else "<UNKNOWN>"}" after scanning a variant spanning multiple records does not match the expected variant length: Found {record_len}, expected {svlen}: Scanned depth table from index {self.index} to {step_index - 1} (inclusive)'

        return (
            sum_depth / svlen,
            sum_align / svlen,
            ','.join(sorted(qry_id)) if len(qry_id) > 0 else np.nan
        )

def update_filter_compound_fields(df, filter_dict, compound_dict):
    """
    Update FILTER and COMPOUND fields in a variant call DataFrame.

    :param df: Variant call DataFrame.
    :param filter_dict: Dict keyed by variant IDs containing sets of FILTER strings.
    :param compound_dict: Dict keyed by variant IDs containing IDs of variants intersecting this variant.

    :return: `df` with FILTER updated and COMPOUND added after the filter column.
    """

    df['FILTER'] = pd.Series(filter_dict).apply(
        lambda vals: ','.join(sorted(vals))
    ).reindex(df.index, fill_value='PASS')

    col_compound = pd.Series(compound_dict).apply(
        lambda vals: ','.join(sorted(vals))
    ).reindex(df.index, fill_value='')

    if 'COMPOUND' in df.columns:
        df['COMPOUND'] = col_compound
    else:
        df.insert(
            [index for col, index in zip(df.columns, range(df.shape[1])) if col == 'FILTER'][0] + 1,
            'COMPOUND',
            col_compound
        )


def apply_compound_filter(df, compound_filter_tree, filter_dict, compound_dict, update=True):
    """
    Apply the compound-variant filter. If a variant intersects one already seen, then "COMPOUND" is added to the filter
    (`filter_dict`) for that variant and the variant it intersected is added to `compound_dict` for that variant. This
    tracks which variants are being removed because they are in a larger event (such as small variants in a deletion)

    :param df: DataFrame to apply filter for.
    :param compound_filter_tree: Tree of existing variants to test intersects against.
    :param filter_dict: Filter dict matching variants in `df` keyed by variant IDs. Contains sets of filters for each
        variant.
    :param compound_dict: Compound dict matching variants in `df` keyed by variant IDs. Contains sets of variant IDs
        "covering" this event (i.e. large variants already in the callset a variant intersected). The result of this
        can be used to reclaim smaller variants filtered by COMPOUND if all the variant IDs in it's `compound_dict`
        entry were removed from the callset.
    :param update: If `True`, update `compound_filter_tree` if a variant is not filtered and does not intersect another
        variant already in `compound_filter_tree`. Set to `False` if variants in `df` should not be intersected with
        future variants to test for compound filtering.
    """

    for index, row in df.sort_values(['SVLEN', 'POS'], ascending=(False, True)).iterrows():
        intersect_set = compound_filter_tree[row['#CHROM']][row['POS']:row['END']]

        if len(intersect_set) == 0:
            if update and index not in filter_dict.keys():
                compound_filter_tree[row['#CHROM']][row['POS']:row['END']] = row['ID']
        else:
            filter_dict[index].add('COMPOUND')
            compound_dict[index] |= {val.data for val in intersect_set}


def apply_qry_filter_tree(df, qry_filter_tree, filter_dict):
    """
    Match QRY_REGION from variant calls to a filter tree in query coordinates.

    :param df: Variant call dataframe.
    :param qry_filter_tree: Filter (dict keyed by query IDs with intervaltree objects covering filtered loci). If
        None, no filtering is applied (skips call).
    :param filter_dict: A dictionary (defaultdict) keyed by variant call IDs with sets of filter names. Filter
        "QRY_FILTER" is added to the set for all IDs intersecting a region in `qry_filter_tree`.
    """

    if qry_filter_tree is not None:

        filter_set = df.apply(lambda row: pavlib.seq.region_from_string(row['QRY_REGION']), axis=1).apply(
            lambda region: len(qry_filter_tree[region.chrom][region.pos:region.end]) > 0
        )

        for index in filter_set[filter_set].index:
            filter_dict[index].add('QRY_FILTER')


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


def merge_haplotypes(bed_list, callable_list, hap_list, config_def, threads=1, subset_chrom=None):
    """
    Merge haplotypes for one variant type.

    :param bed_list: List of input BED files.
    :param callable_list: List of callable reference loci in BED files.
    :param hap_list: List of haplotypes.
    :param config_def: Merge definition.
    :param threads: Number of threads for each merge.
    :param subset_chrom: Chromosome or set of chromosome names to merge, or `None` to merge all chromosomes in one step.
    :param is_inv: Add inversion columns if `True`, autodetect if `None`.

    :return: A dataframe of variant calls.
    """

    # Check input
    n_hap = len(hap_list)

    if len(bed_list) != n_hap:
        raise RuntimeError(f'Input variant BED list length ({len(bed_list)}) does not match the haplotype name list length: ({n_hap})')

    if len(callable_list) != n_hap:
        raise RuntimeError(f'Input callable BED list length ({len(callable_list)}) does not match the haplotype name list length: ({n_hap})')

    # Merge
    df = svpoplib.svmerge.merge_variants(
        bed_list=bed_list,
        sample_names=hap_list,
        strategy=config_def,
        threads=threads,
        subset_chrom=subset_chrom
    )

    df.set_index('ID', inplace=True, drop=False)
    df.index.name = 'INDEX'

    # Restructure columns
    for col in ('HAP', 'RGN_REF_DISC', 'RGN_QRY_DISC', 'FLAG_ID', 'FLAG_TYPE', 'MERGE_SRC', 'MERGE_SRC_ID'):
        if col in df.columns:
            del(df[col])

    df.columns = [re.sub('^MERGE_', 'HAP_', val) for val in df.columns]
    df.columns = ['HAP' if val == 'HAP_SAMPLES' else val for val in df.columns]

    # Change , to ; from merger
    for col in ('HAP', 'HAP_VARIANTS', 'HAP_RO', 'HAP_SZRO', 'HAP_OFFSET', 'HAP_OFFSZ', 'HAP_MATCH'):
        if col in df.columns:
            df[col] = df[col].apply(lambda val: ';'.join(val.split(',')))

    # Pack fields with values from all haplotypes
    df_dict = {
        key: val for key, val in zip(hap_list, bed_list)
    }

    if df.shape[0] > 0:
        for col in ('QRY_REGION', 'QRY_STRAND', 'CI', 'ALIGN_INDEX', 'CALL_SOURCE', 'RGN_REF_INNER', 'RGN_QRY_INNER', 'COV_MEAN', 'COV_PROP', 'COV_QRY'):
            if col in df.columns:
                df[col] = pavlib.call.val_per_hap(df, df_dict, col)

    # Load mapped regions
    map_tree_list = list()

    for index in range(n_hap):
        map_tree = collections.defaultdict(intervaltree.IntervalTree)

        for index, row in pd.read_csv(callable_list[index], sep='\t').iterrows():
            map_tree[row['#CHROM']][row['POS']:row['END']] = True

        map_tree_list.append(map_tree)

    # Get genotypes setting no-call for non-mappable regions
    if df.shape[0] > 0:
        df['GT'] = pd.concat(
            [
                df.apply(pavlib.call.get_gt, hap=hap_list[index], map_tree=map_tree_list[index], axis=1)
                    for index in range(n_hap)
            ],
            axis=1
        ).apply(lambda vals: '|'.join(vals), axis=1)
    else:
        df['GT'] = ''

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
