"""
Call inversions from aligned query sequences.

Inversion calling has two key steps:
1) Flag: Find signatures of inversions from alignments and/or variant calls.
2) Call: Use flagged loci to call inversions. Calling starts with the flagged locus, then it expands until it
   finds the full inversion including unique reference sequence on each flank.
"""

import collections
import intervaltree
import numpy as np
import os
import pandas as pd

import pavlib
import svpoplib
import kanapy

global expand
global temp
global get_config
global REF_FA


#
# Definitions
#

# Column names for the inversion call table
INV_CALL_COLUMNS = [
    '#CHROM', 'POS', 'END',
    'ID', 'SVTYPE', 'SVLEN',
    'HAP',
    'QRY_REGION', 'QRY_STRAND',
    'CI',
    'RGN_REF_INNER', 'RGN_QRY_INNER',
    'RGN_REF_DISC', 'RGN_QRY_DISC',
    'FLAG_ID', 'FLAG_TYPE',
    'ALIGN_INDEX',
    'CALL_SOURCE',
    'FILTER',
    'SEQ'
]


def _input_call_inv_cluster(wildcards):
    """
    Get input for flagging by clustered variant calls. Return both ins & del indels or snvs.

    :param wildcards: Wildcards.

    :return: A list of input files.
    """

    if wildcards.vartype == 'indel':
        return [
            'temp/{asm_name}/cigar/merged/svindel_insdel_{hap}.bed.gz'.format(**wildcards),
        ]

    elif wildcards.vartype == 'snv':
        return [
            'temp/{asm_name}/cigar/merged/snv_snv_{hap}.bed.gz'.format(**wildcards),
        ]


    raise RuntimeError('Unknown variant type {} for input into call_inv_cluster: Expected "indel" or "snv"'.format(
        wildcards.vartype
    ))

def _call_inv_accept_flagged_region(row, allow_single_cluster=False, match_any=set()):
    """
    Annotate which flagged regions are "accepted" (will try to call INV).

    If `allow_single_cluster` is `False` and `match_any` is an empty set, then the only signatures accepted are
    matched SV events or matched indel events.

    :param row: Row from merged flagged regions.
    :param allow_single_cluster: Try to resolve inversions if True for loci with only signatures of clustered SNVs
        and/or clustered indels. Will try many regions and may increase false-positives.
    :param match_any: If defined, contains a set of signatures where at least one must match. Can be used to
        restrict accepted regions to SV-supported signatures only, but inversions with a small uniquely-inverted
        region will likely be missed.

    :return: `True` if the inversion caller should try the region.
    """

    if not allow_single_cluster and (row['TYPE'] == {'CLUSTER_SNV'} or row['TYPE'] == {'CLUSTER_INDEL'}):
        return False

    if match_any and not row['TYPE'] & match_any:
        return False

    return True

BATCH_COUNT_DEFAULT = 60


#############
### Rules ###
#############


#
# Call inversions
#

# Merge batches.
# noinspection PyTypeChecker
rule call_inv_batch_merge:
    input:
        bed=lambda wildcards: expand('temp/{{asm_name}}/inv_caller/batch/{{hap}}/inv_call_{batch}.bed.gz', batch=range(int(get_config(wildcards, 'inv_sig_batch_count', BATCH_COUNT_DEFAULT))))
    output:
        bed=temp('temp/{asm_name}/inv_caller/sv_inv_{hap}.bed.gz')
    run:

        df = pd.concat(
            [pd.read_csv(file_name, sep='\t') for file_name in input.bed],
            axis=0
        )

        df.drop_duplicates('ID', inplace=True)

        df.sort_values(
            ['#CHROM', 'POS', 'END', 'ID']
        ).to_csv(
            output.bed, sep='\t', index=False, compression='gzip'
        )

# Call inversions in batches of flagged regions.
rule call_inv_batch:
    input:
        bed_flag='results/{asm_name}/inv_caller/flagged_regions_{hap}.bed.gz',
        bed_aln='results/{asm_name}/align/trim-qryref/aligned_query_{hap}.bed.gz',
        qry_fa='temp/{asm_name}/align/query_{hap}.fa.gz',
        qry_fai='temp/{asm_name}/align/query_{hap}.fa.gz.fai'
    output:
        bed=temp('temp/{asm_name}/inv_caller/batch/{hap}/inv_call_{batch}.bed.gz')
    log:
        log='log/{asm_name}/inv_caller/log/{hap}/inv_call_{batch}.log'
    params:
        align_score=lambda wildcards: get_config(wildcards, 'align_score_model', pavlib.align.score.DEFAULT_ALIGN_SCORE_MODEL)
    threads: 1
    run:

        call_source = 'FLAG-DEN'

        # Get params
        k_size = get_config(wildcards, 'inv_k_size', pavlib.const.INV_K_SIZE)

        region_limit = get_config(wildcards, 'inv_region_limit', pavlib.const.INV_REGION_LIMIT)
        min_expand = get_config(wildcards, 'inv_min_expand', pavlib.const.INV_MIN_EXPAND_COUNT)
        init_expand = get_config(wildcards, 'inv_init_expand', pavlib.const.INV_INIT_EXPAND)
        min_kmers = get_config(wildcards, 'inv_min_kmers', pavlib.const.INV_MIN_KMERS)
        max_ref_kmer_count = get_config(wildcards, 'inv_max_ref_kmer_count', pavlib.const.INV_MAX_REF_KMER_COUNT)

        kde_bandwidth = get_config(wildcards, 'inv_kde_bandwidth', pavlib.const.INV_KDE_BANDWIDTH)
        kde_trunc_z = get_config(wildcards, 'inv_kde_trunc_z', pavlib.const.INV_KDE_TRUNC_Z)
        kde_func = get_config(wildcards, 'inv_kde_func', pavlib.const.INV_KDE_FUNC)

        batch = int(wildcards.batch)

        density_out_dir = 'results/{asm_name}/inv_caller/density_table'.format(**wildcards)

        os.makedirs(density_out_dir, exist_ok=True)

        # Read and subset table to records in this batch
        df_flag = pd.read_csv(input.bed_flag, sep='\t', header=0)
        df_flag = df_flag.loc[df_flag['BATCH'] == batch]

        if df_flag.shape[0] == 0:
            # No records in batch

            df_bed = pd.DataFrame(
                [], columns=INV_CALL_COLUMNS
            )

        else:

            # Init
            k_util = kanapy.util.kmer.KmerUtil(k_size)

            align_lift = pavlib.align.lift.AlignLift(
                pd.read_csv(input.bed_aln, sep='\t'),
                svpoplib.ref.get_df_fai(input.qry_fai)
            )

            id_set = set()  # The caller can find the same inversion through different flagged regions, track and drop duplicates

            kde = pavlib.kde.KdeTruncNorm(
                kde_bandwidth, kde_trunc_z, kde_func
            )

            # Call inversions
            call_list = list()

            with open(log.log, 'w') as log_file:
                for index, row in df_flag.iterrows():

                    # Scan for inversions
                    region_flag = pavlib.seq.Region(row['#CHROM'], row['POS'], row['END'])

                    try:
                        inv_call = pavlib.inv.scan_for_inv(
                            region_flag, REF_FA, input.qry_fa,
                            align_lift=align_lift,
                            k_util=k_util,
                            nc_ref=None,
                            nc_qry=None,
                            region_limit=region_limit,
                            min_expand=min_expand,
                            init_expand=init_expand,
                            min_kmers=min_kmers,
                            max_ref_kmer_count=max_ref_kmer_count,
                            kde=kde,
                            log=log_file,
                        )


                    except RuntimeError as ex:
                        log_file.write('RuntimeError in scan_for_inv(): {}\n'.format(ex))
                        inv_call = None

                    # Save inversion call
                    if inv_call is not None and inv_call.id not in id_set:

                        # Get seq
                        seq = pavlib.seq.region_seq_fasta(
                            inv_call.region_qry_outer,
                            input.qry_fa,
                            rev_compl=inv_call.region_qry_outer.is_rev
                        )

                        # Get alignment record data
                        align_index = ','.join(sorted(
                            pavlib.util.collapse_to_set(
                                (
                                    inv_call.region_ref_outer.pos_aln_index,
                                    inv_call.region_ref_outer.end_aln_index,
                                    inv_call.region_ref_inner.pos_aln_index,
                                    inv_call.region_ref_inner.end_aln_index
                                ),
                                to_type=str
                            )
                        ))

                        # Save call
                        call_list.append(pd.Series(
                            [
                                inv_call.region_ref_outer.chrom,
                                inv_call.region_ref_outer.pos,
                                inv_call.region_ref_outer.end,

                                inv_call.id,
                                'INV',
                                inv_call.svlen,

                                wildcards.hap,

                                inv_call.region_qry_outer.to_base1_string(),
                                '-' if inv_call.region_qry_outer.is_rev else '+',

                                0,

                                inv_call.region_ref_inner.to_base1_string(),
                                inv_call.region_qry_inner.to_base1_string(),

                                inv_call.region_ref_discovery.to_base1_string(),
                                inv_call.region_qry_discovery.to_base1_string(),

                                inv_call.region_flag.region_id(),
                                row['TYPE'],

                                align_index,

                                pavlib.inv.CALL_SOURCE,

                                'PASS',

                                seq
                            ],
                            index=[
                                '#CHROM', 'POS', 'END',
                                'ID', 'SVTYPE', 'SVLEN',
                                'HAP',
                                'QRY_REGION', 'QRY_STRAND',
                                'CI',
                                'RGN_REF_INNER', 'RGN_QRY_INNER',
                                'RGN_REF_DISC', 'RGN_QRY_DISC',
                                'FLAG_ID', 'FLAG_TYPE',
                                'ALIGN_INDEX',
                                'CALL_SOURCE',
                                'FILTER',
                                'SEQ'
                            ]
                        ))

                        id_set.add(inv_call.id)

                        # Save density table
                        inv_call.df.to_csv(
                            os.path.join(density_out_dir, 'density_{}_{}.tsv.gz'.format(inv_call.id, wildcards.hap)),
                            sep='\t', index=False, compression='gzip'
                        )

            # Merge records
            if len(call_list) > 0:
                df_bed = pd.concat(call_list, axis=1).T.sort_values(['#CHROM', 'POS', 'END', 'ID'])

            else:
                # Create emtpy data frame
                df_bed = pd.DataFrame(
                    [],
                    columns=[
                        '#CHROM', 'POS', 'END',
                        'ID', 'SVTYPE', 'SVLEN',
                        'HAP',
                        'QRY_REGION', 'QRY_STRAND',
                        'CI',
                        'RGN_REF_INNER', 'RGN_QRY_INNER',
                        'RGN_REF_DISC', 'RGN_QRY_DISC',
                        'FLAG_ID', 'FLAG_TYPE',
                        'ALIGN_INDEX',
                        'CALL_SOURCE',
                        'SEQ'
                    ]
                )

        # Write
        df_bed.to_csv(output.bed, sep='\t', index=False, compression='gzip')


#
# Flag regions to scan for inversions
#

# Merge flagged overlapping regions (SV ins/del match, indel ins/del match, indel cluster, SNV cluster) into single
# flagged record regions and annotate each. Column "TRY_INV" will be used by the inversion caller to determine which
# regions to try calling.
rule call_inv_merge_flagged_loci:
    input:
        bed_insdel_sv='temp/{asm_name}/inv_caller/flag/insdel_sv_{hap}.bed.gz',
        bed_insdel_indel='temp/{asm_name}/inv_caller/flag/insdel_indel_{hap}.bed.gz',
        bed_cluster_indel='temp/{asm_name}/inv_caller/flag/cluster_indel_{hap}.bed.gz',
        bed_cluster_snv='temp/{asm_name}/inv_caller/flag/cluster_snv_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/inv_caller/flagged_regions_{hap}.bed.gz'
    run:

        # Parameters
        flank = int(get_config(wildcards, 'inv_sig_merge_flank', 500))  # Merge windows within this many bp
        batch_count = int(get_config(wildcards, 'inv_sig_batch_count', BATCH_COUNT_DEFAULT))  # Batch signature regions into this many batches for the caller. Marked here so that this file can be cross-referenced with the inversion caller log
        inv_sig_filter = get_config(wildcards, 'inv_sig_filter', 'svindel')   # Filter flagged regions

        # Get region filter parameters
        allow_single_cluster = False
        match_any = set()

        if inv_sig_filter is not None:
            if inv_sig_filter == 'single_cluster':
                allow_single_cluster = True

            elif inv_sig_filter == 'svindel':
                match_any.add('MATCH_SV')
                match_any.add('MATCH_INDEL')

            elif inv_sig_filter == 'sv':
                match_any.add('MATCH_SV')

            else:
                raise RuntimeError(f'Unrecognized region filter: {inv_sig_filter} (must be "single_cluster", "svindel", or "sv")')

        # Read
        df_insdel_sv = pd.read_csv(input.bed_insdel_sv, sep='\t')
        df_insdel_indel = pd.read_csv(input.bed_insdel_indel, sep='\t')
        df_cluster_indel = pd.read_csv(input.bed_cluster_indel, sep='\t')
        df_cluster_snv = pd.read_csv(input.bed_cluster_snv, sep='\t')

        # Annotate counts
        df_insdel_sv['COUNT_INDEL'] = 0
        df_insdel_sv['COUNT_SNV'] = 0
        df_insdel_sv['TYPE'] = [{'MATCH_SV'} for index in range(df_insdel_sv.shape[0])]

        df_insdel_indel['COUNT_INDEL'] = 0
        df_insdel_indel['COUNT_SNV'] = 0
        df_insdel_indel['TYPE'] = [{'MATCH_INDEL'} for index in range(df_insdel_indel.shape[0])]

        df_cluster_indel['COUNT_INDEL'] = df_cluster_indel['COUNT']
        df_cluster_indel['COUNT_SNV'] = 0
        df_cluster_indel['TYPE'] = [{'CLUSTER_INDEL'} for index in range(df_cluster_indel.shape[0])]
        del(df_cluster_indel['COUNT'])

        df_cluster_snv['COUNT_INDEL'] = 0
        df_cluster_snv['COUNT_SNV'] = df_cluster_snv['COUNT']
        df_cluster_snv['TYPE'] = [{'CLUSTER_SNV'} for index in range(df_cluster_snv.shape[0])]
        del(df_cluster_snv['COUNT'])

        # Merge
        df = pd.concat([
            df_insdel_sv,
            df_insdel_indel,
            df_cluster_indel,
            df_cluster_snv
        ], axis=0).sort_values(['#CHROM', 'POS'])

        # Merge flagged regions
        region_list = list()

        chrom = None
        pos = 0
        end = 0

        indel_count = 0
        snv_count = 0

        type_set = set()

        for index, row in df.iterrows():

            if (row['POS'] < end + flank) and (row['#CHROM'] == chrom):

                # Add to existing region
                type_set |= row['TYPE']
                end = row['END']

                indel_count += row['COUNT_INDEL']
                snv_count += row['COUNT_SNV']

            else:

                # Write region
                if type_set:
                    region_list.append(pd.Series(
                        [
                            chrom, pos, end,
                            '{}-{}-RGN-{}'.format(chrom, pos, end - pos),
                            'RGN', end - pos,
                            type_set,
                            indel_count, snv_count
                        ],
                        index=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'TYPE', 'COUNT_INDEL', 'COUNT_SNV']
                    ))

                # Start new region
                type_set = row['TYPE'].copy()
                pos = row['POS']
                end = row['END']
                chrom = row['#CHROM']

                indel_count = row['COUNT_INDEL']
                snv_count = row['COUNT_SNV']

        # Final region
        if type_set:
            region_list.append(pd.Series(
                [
                    chrom, pos, end,
                    '{}-{}-RGN-{}'.format(chrom, pos, end - pos),
                    'RGN', end - pos,
                    type_set,
                    #','.join(sorted(type_set)),
                    indel_count, snv_count
                ],
                index=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'TYPE', 'COUNT_INDEL', 'COUNT_SNV']
            ))

        # Merge
        if len(region_list) > 0:
            df_merged = pd.concat(region_list, axis=1).T.sort_values(['#CHROM', 'POS'])

            # Annotate accepted regions
            df_merged['TRY_INV'] = df_merged.apply(
                _call_inv_accept_flagged_region, allow_single_cluster=allow_single_cluster, match_any=match_any,
                axis=1
            )

            # Group into batches
            df_merged['BATCH'] = -1

            batch = 0

            for index, row in df_merged.iterrows():
                if row['TRY_INV']:
                    df_merged.loc[index, 'BATCH'] = batch
                    batch = (batch + 1) % batch_count

        else:
            df_merged = pd.DataFrame([], columns=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'TYPE', 'COUNT_INDEL', 'COUNT_SNV', 'TRY_INV', 'BATCH'])

        df_merged['TYPE'] = df_merged['TYPE'].apply(lambda vals: ','.join(sorted(vals)))

        # Write
        df_merged.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# Flag inversion regions by matched INS/DEL calls, which often occur inside inversions. Aligners will generate a
# deletion over the inversion with an insertion of the same size flanking it where the inversion is the inverted
# sequence.
rule call_inv_flag_insdel_cluster:
    input:
        bed='temp/{asm_name}/cigar/merged/svindel_insdel_{hap}.bed.gz'
    output:
        bed=temp('temp/{asm_name}/inv_caller/flag/insdel_{vartype}_{hap}.bed.gz')
    params:
        flank_cluster=lambda wildcards: int(get_config(wildcards, 'inv_sig_insdel_cluster_flank', 2)), # For each INS, multiply SVLEN by this to get the distance to the nearest DEL it may intersect
        flank_merge=lambda wildcards: int(get_config(wildcards, 'inv_sig_insdel_merge_flank', 2000)),  # Merge clusters within this distance
        cluster_min_svlen=lambda wildcards: int(get_config(wildcards, 'inv_sig_cluster_svlen_min', 4))    # Discard indels less than this size
    wildcard_constraints:
        vartype='sv|indel'
    run:

        # Params
        flank_cluster = params.flank_cluster
        flank_merge = params.flank_merge
        svlen_min = params.cluster_min_svlen if wildcards.vartype == 'indel' else 50

        # Input
        df = pd.read_csv(input.bed, sep='\t', header=0, low_memory=False)
        df = df.loc[df['FILTER'] == 'PASS']
        df = df.loc[df['SVLEN'] >= svlen_min]

        if wildcards.vartype == 'indel':
            df = df.loc[df['SVLEN'] < 50]

        # Stop if variants are empty
        if df.shape[0] == 0:
            pd.DataFrame(
                [],
                columns=['#CHROM', 'POS', 'END']
            ).to_csv(output.bed, sep='\t', index=False, compression='gzip')

            return

        # Subset by SVTYPE
        df_ins = df.loc[df['SVTYPE'] == 'INS']
        df_del = df.loc[df['SVTYPE'] == 'DEL']

        # Load deletion intervals into a tree
        deltree = collections.defaultdict(intervaltree.IntervalTree)

        for index, row in df_del.iterrows():
            deltree[row['#CHROM']][row['POS']:row['END']] = (row['ID'], row['SVLEN'], row['POS'], row['END'])

        # Flag matched INS/DELs
        match_list = list()

        for index, row in df_ins.iterrows():
            flank = row['SVLEN'] * flank_cluster

            match_set = deltree[row['#CHROM']][row['POS'] - flank : row['POS'] + flank]

            if match_set:
                match_list.append(pd.Series(
                    [
                        row['#CHROM'],
                        np.min([record.data[2] for record in match_set]),
                        np.max([record.data[3] for record in match_set])
                    ],
                    index=['#CHROM', 'POS', 'END']
                ))

        if len(match_list) == 0:
            pd.DataFrame(
                [],
                columns=['#CHROM', 'POS', 'END']
            ).to_csv(output.bed, sep='\t', index=False, compression='gzip')

            return

        # Merge overlapping intervals
        df_match = pd.concat(match_list, axis=1).T.sort_values(['#CHROM', 'POS'])

        match_merged_list = list()

        chrom = None
        pos = None
        end = None

        for index, row in df_match.iterrows():

            # Switch chromosomes
            if chrom != row['#CHROM']:

                # Record current record
                if chrom is not None:
                    match_merged_list.append(pd.Series(
                        [chrom, pos, end],
                        index=['#CHROM', 'POS', 'END']
                    ))

                chrom = row['#CHROM']
                pos = row['POS']
                end = row['END']

            # Same chromosome, check record
            if row['POS'] - flank_merge <= end:
                end = np.max([end, row['END']])

            else:
                match_merged_list.append(pd.Series(
                    [chrom, pos, end],
                    index=['#CHROM', 'POS', 'END']
                ))

                pos = row['POS']
                end = row['END']

        if len(match_merged_list) > 0:
            df_match = pd.concat(match_merged_list, axis=1).T.sort_values(['#CHROM', 'POS'])

        else:
            df_match = pd.DataFrame(
                [],
                columns=['#CHROM', 'POS', 'END']
            )

        # Write
        df_match.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# Detect clusters of indels and SNVs that common occur when query sequences are aligned through inversions in direct
# orientation
#
# noinspection PyTypeChecker
rule call_inv_cluster:
    input:
        bed=_input_call_inv_cluster
    output:
        bed=temp('temp/{asm_name}/inv_caller/flag/cluster_{vartype}_{hap}.bed.gz')
    params:
        cluster_win=lambda wildcards: get_config(wildcards, 'inv_sig_cluster_win', 200),            # Cluster variants within this many bases
        cluster_win_min=lambda wildcards: get_config(wildcards, 'inv_sig_cluster_win_min', 500),    # Window must reach this size
        cluster_min_snv=lambda wildcards: get_config(wildcards, 'inv_sig_cluster_snv_min', 20),     # Minimum number if SNVs in window (if vartype == snv)
        cluster_min_indel=lambda wildcards: get_config(wildcards, 'inv_sig_cluster_indel_min', 10)  # Minimum number of indels in window (if vartype == indel)
    wildcard_constraints:
        vartype='indel|snv'
    run:

        # Params
        cluster_win = params.cluster_win
        cluster_win_min = params.cluster_win

        if wildcards.vartype == 'indel':
            cluster_min = params.cluster_min_indel
        elif wildcards.vartype == 'snv':
            cluster_min = params.cluster_min_snv
        else:
            raise RuntimeError('Bad variant type {}: Expected "indel" or "snv"')

        # Read
        df = pd.concat(
            [
                pd.read_csv(
                    input_file_name, sep='\t', usecols=('#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'FILTER'),
                    low_memory=False, dtype={'#CHROM': str}
                ) for input_file_name in input.bed
            ],
            axis=0
        ).sort_values(['#CHROM', 'POS'])

        df = df.loc[df['FILTER'] == 'PASS']

        if wildcards.vartype == 'indel':
            df = df.loc[df['SVLEN'] < 50]

        # DEL to midpoint
        df['POS'] = (df['END'] + df['POS']) // 2

        # Find clusters
        cluster_list = list()

        chrom = None
        cluster_pos = 0
        cluster_end = 0
        cluster_count = 0

        for index, row in df.iterrows():

            if (row['POS'] < cluster_end + cluster_win) and (row['#CHROM'] == chrom):

                # Add to existing cluster
                cluster_count += 1
                cluster_end = row['POS']

            else:

                # Save last cluster
                if (cluster_count >= cluster_min) and (cluster_end - cluster_pos >= cluster_win_min):
                    cluster_list.append(pd.Series(
                        [chrom, cluster_pos, cluster_end, cluster_count],
                        index=['#CHROM', 'POS', 'END', 'COUNT']
                    ))

                # Start new cluster
                cluster_count = 1
                cluster_pos = row['POS']
                cluster_end = row['POS']
                chrom = row['#CHROM']

        # Save final cluster
        if (cluster_count >= cluster_min) and (cluster_end - cluster_pos >= cluster_win_min):
            cluster_list.append(pd.Series(
                [chrom, cluster_pos, cluster_end, cluster_count],
                index=['#CHROM', 'POS', 'END', 'COUNT']
            ))

        # Merge records
        if len(cluster_list) > 0:
            df_cluster = pd.concat(cluster_list, axis=1).T
        else:
            df_cluster = pd.DataFrame([], columns=['#CHROM', 'POS', 'END', 'COUNT'])

        # Write
        df_cluster.to_csv(output.bed, sep='\t', index=False, compression='gzip')
