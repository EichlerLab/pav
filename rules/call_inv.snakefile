"""
Call inversions from aligned contigs.

Inversion calling has two key steps:
1) Flag: Find signatures of inversions from alignments and/or variant calls.
2) Call: Use flagged loci to call inversions. Calling starts with the flagged locus, then it expands until it
   finds the full inversion including unique reference sequence on each flank.
"""

###################
### Definitions ###
###################


def _input_call_inv_cluster(wildcards):
    """
    Get input for flagging by clustered variant calls. Return both ins & del indels or snvs.

    :param wildcards: Wildcards.

    :return: A list of input files.
    """

    if wildcards.vartype == 'indel':
        return [
            'temp/{asm_name}/cigar/pre_inv/svindel_insdel_{hap}.bed.gz'.format(**wildcards),
        ]

    elif wildcards.vartype == 'snv':
        return [
            'temp/{asm_name}/cigar/pre_inv/snv_snv_{hap}.bed.gz'.format(**wildcards),
        ]


    raise RuntimeError('Unknown variant type {} for input into call_inv_cluster: Expected "indel" or "snv"'.format(
        wildcards.vartype
    ))

def _call_inv_accept_flagged_region(row):
    """
    Annotate which flagged regions are "accepted" (will try to call INV).

    :param row: Row from merged flagged regions.

    :return: `True` if the inversion caller should try the region.
    """

    if row['TYPE'] in {'CLUSTER_SNV', 'CLUSTER_INDEL'}:
        return False

    return True

BATCH_COUNT = config.get('inv_sig_batch_count', 60)


#############
### Rules ###
#############

#
# Inter-inversion variants
#

# call_inv_interinv
#
# Call inter-inversion variants.
# Implementation on hold.
# rule call_inv_interinv:
#     input:
#         bed_cigar='results/{asm_name}/inv_caller/sv_inv_{hap}.bed.gz',
#         tig_fa='results/{asm_name}/align/contigs_{hap}.fa.gz'
#     output:
#         bed='results/{asm_name}/inv_caller/sv_inv_{hap}.bed.gz'
#     run:
#         pass
#
#         # Read SVs
#         df_inv = pd.read_csv(input.bed_cigar, sep='\t')
#
#         # Get alignments
#         for index, row in df_inv.iterrows():
#
#             # Get regions
#             region_ref = asmlib.seq.Region(row['#CHROM'], row['POS'], row['END'])
#             region_tig = asmlib.seq.region_from_string(row['RGN_TIG_OUTER'], is_rev=(row['STRAND'] == '-'))
#
#             # Get sequence
#             seq_ref = asmlib.seq.region_seq_fasta(region_ref, REF_FA)
#             seq_tig = asmlib.seq.region_seq_fasta(region_tig, input.tig_fa)

#
# Call inversions
#

# call_inv_batch_merge
#
# Merge batches.
rule call_inv_batch_merge:
    input:
        bed=expand('temp/{{asm_name}}/inv_caller/batch/{{hap}}/inv_call_{batch}.bed.gz', batch=range(BATCH_COUNT))
    output:
        bed='results/{asm_name}/inv_caller/sv_inv_{hap}.bed.gz'
    run:

        pd.concat(
            [pd.read_csv(file_name, sep='\t') for file_name in input.bed],
            axis=0
        ).sort_values(
            ['#CHROM', 'POS']
        ).to_csv(
            output.bed, sep='\t', index=False, compression='gzip'
        )

# call_inv_batch
#
# Call inversions in batches of flagged regions.
rule call_inv_batch:
    input:
        bed_flag='results/{asm_name}/inv_caller/flagged_regions_{hap}.bed.gz',
        bed_aln='results/{asm_name}/align/aligned_tig_{hap}.bed.gz',
        tig_fa='results/{asm_name}/align/contigs_{hap}.fa.gz',
        fai='results/{asm_name}/align/contigs_{hap}.fa.gz.fai'
    output:
        bed=temp('temp/{asm_name}/inv_caller/batch/{hap}/inv_call_{batch}.bed.gz')
    log:
        log='results/{asm_name}/inv_caller/batch/{hap}/inv_call_{batch}.log'
    params:
        k_size=config.get('inv_k_size', 31),
        inv_threads=config.get('inv_threads', 12),
        inv_mem=config.get('inv_mem', '4G')
    run:

        # Get params
        batch = int(wildcards.batch)
        k_size = params.k_size

        density_out_dir = 'results/{asm_name}/inv_caller/density_table'.format(**wildcards)

        os.makedirs(density_out_dir, exist_ok=True)

        # Get SRS (state-run-smooth)
        srs_tree = asmlib.inv.get_srs_tree(config.get('srs_list'), None)  # If none, tree contains a default for all region sizes

        # Read and subset table to records in this batch
        df_flag = pd.read_csv(input.bed_flag, sep='\t', header=0)
        df_flag = df_flag.loc[df_flag['BATCH'] == batch]

        if df_flag.shape[0] == 0:
            # No records in batch

            df_bed = pd.DataFrame(
                [],
                columns=[
                    '#CHROM', 'POS', 'END',
                    'ID', 'SVTYPE', 'SVLEN',
                    'HAP',
                    'TIG_REGION', 'QUERY_STRAND',
                    'CI',
                    'RGN_REF_INNER', 'RGN_TIG_INNER',
                    'RGN_REF_DISC', 'RGN_TIG_DISC',
                    'FLAG_ID', 'FLAG_TYPE',
                    'ALIGN_INDEX', 'CLUSTER_MATCH',
                    'CALL_SOURCE',
                    'SEQ'
                ]
            )

        else:

            # Init
            k_util = kanapy.util.kmer.KmerUtil(k_size)

            align_lift = asmlib.align.AlignLift(
                pd.read_csv(input.bed_aln, sep='\t'),
                analib.ref.get_df_fai(input.fai)
            )

            # Read alignment BED
            df_aln = pd.read_csv(
                input.bed_aln,
                sep='\t',
                usecols=('INDEX', 'CLUSTER_MATCH'),
                index_col='INDEX',
                squeeze=False
            )


            # Call inversions
            call_list = list()

            with open(log.log, 'w') as log_file:
                for index, row in df_flag.iterrows():

                    # Scan for inversions
                    region_flag = asmlib.seq.Region(row['#CHROM'], row['POS'], row['END'])

                    try:
                        inv_call = asmlib.inv.scan_for_inv(
                            region_flag, REF_FA, input.tig_fa, align_lift, k_util,
                            threads=params.inv_threads, log=log_file, srs_tree=srs_tree
                        )

                    except RuntimeError as ex:
                        log_file.write('RuntimeError in scan_for_inv(): {}\n'.format(ex))
                        inv_call = None

                    # Save inversion call
                    if inv_call is not None:

                        # Get seq
                        seq = asmlib.seq.region_seq_fasta(
                            inv_call.region_tig_outer,
                            input.tig_fa,
                            rev_compl=inv_call.region_tig_outer.is_rev
                        )

                        # Get alignment record data
                        aln_index_set = {
                            val for val_list in [
                                inv_call.region_ref_outer.pos_aln_index,
                                inv_call.region_ref_outer.end_aln_index,
                                inv_call.region_ref_inner.pos_aln_index,
                                inv_call.region_ref_inner.end_aln_index

                            ] for val in val_list
                        }

                        cluster_match_set = {df_aln.loc[index, 'CLUSTER_MATCH'] for index_tuple in aln_index_set for index in index_tuple}

                        if False in cluster_match_set:
                            cluster_match = False
                        elif np.any(pd.isnull(list(cluster_match_set))):
                            cluster_match = np.nan
                        else:
                            if cluster_match_set != {True}:
                                raise RuntimeError(
                                    'Found unexpected values in cluster_match: Expected True, False, and np.nan: {}'.format(
                                        ', '.join([str(val) for val in cluster_match_set])
                                    )
                                )

                            cluster_match = True

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

                                inv_call.region_tig_outer.to_base1_string(),
                                '-' if inv_call.region_tig_outer.is_rev else '+',

                                0,

                                inv_call.region_ref_inner.to_base1_string(),
                                inv_call.region_tig_inner.to_base1_string(),

                                inv_call.region_ref_discovery.to_base1_string(),
                                inv_call.region_tig_discovery.to_base1_string(),

                                inv_call.region_flag.region_id(),
                                row['TYPE'],

                                ','.join([str(val) for val in sorted(aln_index_set)]), cluster_match,

                                asmlib.inv.CALL_SOURCE,

                                seq
                            ],
                            index=[
                                '#CHROM', 'POS', 'END',
                                'ID', 'SVTYPE', 'SVLEN',
                                'HAP',
                                'TIG_REGION', 'QUERY_STRAND',
                                'CI',
                                'RGN_REF_INNER', 'RGN_TIG_INNER',
                                'RGN_REF_DISC', 'RGN_TIG_DISC',
                                'FLAG_ID', 'FLAG_TYPE',
                                'ALIGN_INDEX', 'CLUSTER_MATCH',
                                'CALL_SOURCE',
                                'SEQ'
                            ]
                        ))

                        # Save density table
                        inv_call.df.to_csv(
                            os.path.join(density_out_dir, 'density_{}_{}.tsv.gz'.format(inv_call.id, wildcards.hap)),
                            sep='\t', index=False, compression='gzip'
                        )

            # Merge records
            if len(call_list) > 0:
                df_bed = pd.concat(call_list, axis=1).T

            else:
                # Create emtpy data frame
                df_bed = pd.DataFrame(
                    [],
                    columns=[
                        '#CHROM', 'POS', 'END',
                        'ID', 'SVTYPE', 'SVLEN',
                        'HAP',
                        'TIG_REGION', 'QUERY_STRAND',
                        'CI',
                        'RGN_REF_INNER', 'RGN_TIG_INNER',
                        'RGN_REF_DISC', 'RGN_TIG_DISC',
                        'FLAG_ID', 'FLAG_TYPE',
                        'ALIGN_INDEX', 'CLUSTER_MATCH',
                        'CALL_SOURCE',
                        'SEQ'
                    ]
                )

        # Write
        df_bed.to_csv(output.bed, sep='\t', index=False, compression='gzip')


#
# Flag regions to scan for inversions
#

# call_inv_merge_flagged_loci
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
    params:
        flank=config.get('inv_sig_merge_flank', 500) , # Merge windows within this many bp
        batch_count=config.get('inv_sig_batch_count', 40)  # Batch signature regions into this many batches for the caller. Marked here so that this file can be cross-referenced with the inversion caller log
    run:
        # Parameters
        flank = params.flank

        # Read
        df_indsdel_sv = pd.read_csv(input.bed_insdel_sv, sep='\t')
        df_indsdel_indel = pd.read_csv(input.bed_insdel_indel, sep='\t')
        df_cluster_indel = pd.read_csv(input.bed_cluster_indel, sep='\t')
        df_cluster_snv = pd.read_csv(input.bed_cluster_snv, sep='\t')


        # Annotate counts
        df_indsdel_sv['COUNT_INDEL'] = 0
        df_indsdel_sv['COUNT_SNV'] = 0
        df_indsdel_sv['TYPE'] = 'MATCH_SV'

        df_indsdel_indel['COUNT_INDEL'] = 0
        df_indsdel_indel['COUNT_SNV'] = 0
        df_indsdel_indel['TYPE'] = 'MATCH_INDEL'

        df_cluster_indel['COUNT_INDEL'] = df_cluster_indel['COUNT']
        df_cluster_indel['COUNT_SNV'] = 0
        df_cluster_indel['TYPE'] = 'CLUSTER_INDEL'
        del(df_cluster_indel['COUNT'])

        df_cluster_snv['COUNT_INDEL'] = 0
        df_cluster_snv['COUNT_SNV'] = df_cluster_snv['COUNT']
        df_cluster_snv['TYPE'] = 'CLUSTER_SNV'
        del(df_cluster_snv['COUNT'])

        # Merge
        df = pd.concat([
            df_indsdel_sv,
            df_indsdel_indel,
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
                type_set.add(row['TYPE'])
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
                            ','.join(sorted(type_set)),
                            indel_count, snv_count
                        ],
                        index=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'TYPE', 'COUNT_INDEL', 'COUNT_SNV']
                    ))

                # Start new region
                type_set = {row['TYPE']}
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
                    ','.join(sorted(type_set)),
                    indel_count, snv_count
                ],
                index=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'TYPE', 'COUNT_INDEL', 'COUNT_SNV']
            ))

        # Merge
        if len(region_list) > 0:
            df_merged = pd.concat(region_list, axis=1).T.sort_values(['#CHROM', 'POS'])

            # Annotate accepted regions
            df_merged['TRY_INV'] = df_merged.apply(_call_inv_accept_flagged_region, axis=1)

            # Group into batches
            df_merged['BATCH'] = -1

            batch = 0

            for index, row in df_merged.iterrows():
                if row['TRY_INV']:
                    df_merged.loc[index, 'BATCH'] = batch
                    batch = (batch + 1) % BATCH_COUNT

        else:
            df_merged = pd.DataFrame([], columns=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'TYPE', 'COUNT_INDEL', 'COUNT_SNV', 'TRY_INV', 'BATCH'])

        # Write
        df_merged.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# call_inv_flag_insdel_cluster
#
# Flag inversion regions by matched INS/DEL calls, which often occur inside inversions. Aligners will generate a
# deletion over the inversion with an insertion of the same size flanking it where the inversion is the inverted
# sequence.
rule call_inv_flag_insdel_cluster:
    input:
        bed='temp/{asm_name}/cigar/pre_inv/svindel_insdel_{hap}.bed.gz'
    output:
        bed=temp('temp/{asm_name}/inv_caller/flag/insdel_{vartype}_{hap}.bed.gz')
    params:
        flank_cluster=int(config.get('inv_sig_insdel_cluster_flank', 2)), # For each INS, multiply SVLEN by this to get the distance to the nearest DEL it may intersect
        flank_merge=int(config.get('inv_sig_insdel_merge_flank', 2000)),  # Merge clusters within this distance
        cluster_min_svlen=int(config.get('inv_sig_cluster_svlen_min', 4))    # Discard indels less than this size
    wildcard_constraints:
        vartype='sv|indel'
    run:

        # Params
        flank_cluster = params.flank_cluster
        flank_merge = params.flank_merge
        svlen_min = params.cluster_min_svlen if wildcards.vartype == 'indel' else 50

        # Input
        df = pd.read_csv(input.bed, sep='\t', header=0, low_memory=False)
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

        df_match = pd.concat(match_merged_list, axis=1).T.sort_values(['#CHROM', 'POS'])

        # Write
        df_match.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# call_inv_cluster
#
# Detect clusters of indels and SNVs that common occur when contigs are aligned through inversions in direct orientation
rule call_inv_cluster:
    input:
        bed=_input_call_inv_cluster
    output:
        bed=temp('temp/{asm_name}/inv_caller/flag/cluster_{vartype}_{hap}.bed.gz')
    params:
        cluster_win=config.get('inv_sig_cluster_win', 200),            # Cluster variants within this many bases
        cluster_win_min=config.get('inv_sig_cluster_win_min', 500),    # Window must reach this size
        cluster_min_snv=config.get('inv_sig_cluster_snv_min', 20),     # Minimum number if SNVs in window (if vartype == snv)
        cluster_min_indel=config.get('inv_sig_cluster_indel_min', 10)  # Minimum number of indels in window (if vartype == indel)
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
            [pd.read_csv(input_file_name, sep='\t', usecols=('#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN')) for input_file_name in input.bed],
            axis=0
        ).sort_values(['#CHROM', 'POS'])

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

