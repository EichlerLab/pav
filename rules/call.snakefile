"""
Call variants from aligned contigs.
"""


#
# Definitions
#

CALL_CIGAR_BATCH_COUNT = 10


#
# Finalize variant calls
#

# call_final_bed
#
# Separate variants into final BED and FA files and write to results.
rule call_final_bed:
    input:
        bed_ins='temp/{asm_name}/bed/merged/svindel_ins.bed.gz',
        bed_del='temp/{asm_name}/bed/merged/svindel_del.bed.gz',
        bed_inv='temp/{asm_name}/bed/merged/sv_inv.bed.gz',
        bed_snv='temp/{asm_name}/bed/merged/snv_snv.bed.gz'
    output:
        bed_snv_snv='results/{asm_name}/bed/snv_snv.bed.gz',
        bed_indel_ins='results/{asm_name}/bed/indel_ins.bed.gz',
        bed_indel_del='results/{asm_name}/bed/indel_del.bed.gz',
        bed_sv_ins='results/{asm_name}/bed/sv_ins.bed.gz',
        bed_sv_del='results/{asm_name}/bed/sv_del.bed.gz',
        bed_sv_inv='results/{asm_name}/bed/sv_inv.bed.gz',
        fa_indel_ins='results/{asm_name}/bed/fa/indel_ins.fa.gz',
        fa_indel_del='results/{asm_name}/bed/fa/indel_del.fa.gz',
        fa_sv_ins='results/{asm_name}/bed/fa/sv_ins.fa.gz',
        fa_sv_del='results/{asm_name}/bed/fa/sv_del.fa.gz',
        fa_sv_inv='results/{asm_name}/bed/fa/sv_inv.fa.gz'
    run:

        # Process INS/DEL/INV (SV and indel)
        df_ins = pd.read_csv(input.bed_ins, sep='\t', low_memory=False)
        df_del = pd.read_csv(input.bed_del, sep='\t', low_memory=False)
        df_inv = pd.read_csv(input.bed_inv, sep='\t', low_memory=False)

        df_svtype_dict = {
            'ins': df_ins,
            'del': df_del,
            'inv': df_inv
        }

        for vartype, svtype in [('sv', 'ins'), ('sv', 'del'), ('sv', 'inv'), ('indel', 'ins'), ('indel', 'del')]:

            df = df_svtype_dict[svtype]

            # Subset
            if vartype == 'sv':
                df = df.loc[df['SVLEN'] >= 50]

            elif vartype == 'indel':
                df = df.loc[df['SVLEN'] < 50]

            else:
                raise RuntimeError('Program Bug: Unknown variant type: {}'.format(vartype))

            # Get output file names
            bed_file_name = output[f'bed_{vartype}_{svtype}']
            fa_file_name = output[f'fa_{vartype}_{svtype}']

            # Write FASTA
            with Bio.bgzf.BgzfWriter(fa_file_name, 'wb') as out_file:
                SeqIO.write(analib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')

            del(df['SEQ'])

            # Arrange columns
            df = analib.variant.order_variant_columns(
                df,
                head_cols=[
                    '#CHROM', 'POS', 'END', 'ID',
                    'SVTYPE', 'SVLEN', 'REF', 'ALT',
                    'HAP', 'GT', 'CLUSTER_MATCH', 'CALL_SOURCE',
                    'TIG_REGION', 'QUERY_STRAND', 'ALIGN_INDEX'
                    'CI',
                ],
                tail_cols=[
                    'CALL_SOURCE',
                    'HAP_VARIANTS',
                    'CI',
                    'HAP_RO', 'HAP_OFFSET', 'HAP_SZRO', 'HAP_OFFSZ'
                ],
                allow_missing=True
            )

            # Write BED
            df.to_csv(bed_file_name, sep='\t', index=False, compression='gzip')

        # Process SNVs
        df = pd.read_csv(input.bed_snv, sep='\t', low_memory=False)
        df.to_csv(output.bed_snv_snv, sep='\t', index=False, compression='gzip')


# pg_variant_bed
#
# Make variant BED and FASTA. Write all variant calls regardless of consensus loci.
rule call_merge_haplotypes:
    input:
        bed_var_h1='temp/{asm_name}/bed/integrated/h1/{vartype_svtype}.bed.gz',
        bed_var_h2='temp/{asm_name}/bed/integrated/h2/{vartype_svtype}.bed.gz',
        callable_h1='results/{asm_name}/callable/callable_regions_h1_500.bed.gz',
        callable_h2='results/{asm_name}/callable/callable_regions_h2_500.bed.gz'
    output:
        bed=temp('temp/{asm_name}/bed/merged/{vartype_svtype}.bed.gz')
    params:
        ro_min=float(config.get('ro_min', 0.5)),
        offset_max=int(config.get('offset_max', 200)),
        merge_threads=int(config.get('merge_threads', 12)),
        merge_mem=config.get('merge_mem', '2G')
    run:

        # Get configured merge definition
        if wildcards.vartype_svtype == 'snv_snv':
            config_def = 'nrid'
        else:
            config_def = 'nr:szro={}:offset={}'.format(int(params.ro_min * 100), params.offset_max)

        print('Merging with def: ' + config_def)
        sys.stdout.flush()

        # Merge
        df = analib.svmerge.merge_variants(
            bed_list=[input.bed_var_h1, input.bed_var_h2],
            sample_names=['h1', 'h2'],
            strategy=config_def,
            threads=params.merge_threads
        )

        df.set_index('ID', inplace=True, drop=False)

        # Restructure columns
        del(df['HAP'])
        del(df['DISC_CLASS'])

        df.columns = [re.sub('^MERGE_', 'HAP_', val) for val in df.columns]

        del(df['HAP_SRC'])
        del(df['HAP_SRC_ID'])
        del(df['HAP_AC'])
        del(df['HAP_AF'])

        df.columns = ['HAP' if val == 'HAP_SAMPLES' else val for val in df.columns]

        # Change , to ; from merger
        df['HAP'] = df['HAP'].apply(lambda val: ';'.join(val.split(',')))
        df['HAP_VARIANTS'] = df['HAP_VARIANTS'].apply(lambda val: ';'.join(val.split(',')))

        if 'HAP_RO' in df.columns:
            df['HAP_RO'] = df['HAP_RO'].apply(lambda val: ';'.join(val.split(',')))

        if 'HAP_OFFSET' in df.columns:
            df['HAP_OFFSET'] = df['HAP_OFFSET'].apply(lambda val: ';'.join(val.split(',')))

        if 'HAP_SZRO' in df.columns:
            df['HAP_SZRO'] = df['HAP_SZRO'].apply(lambda val: ';'.join(val.split(',')))

        if 'HAP_OFFSZ' in df.columns:
            df['HAP_OFFSZ'] = df['HAP_OFFSZ'].apply(lambda val: ';'.join(val.split(',')))

        # Add h1 and h2 to columns
        df_h1 = pd.read_csv(input.bed_var_h1, sep='\t', low_memory=False)
        df_h1.set_index('ID', inplace=True, drop=False)
        df_h1['CLUSTER_MATCH'].fillna('NA', inplace=True)
        df_h1 = df_h1.astype(str)

        df_h2 = pd.read_csv(input.bed_var_h2, sep='\t', low_memory=False)
        df_h2.set_index('ID', inplace=True, drop=False)
        df_h2['CLUSTER_MATCH'].fillna('NA', inplace=True)
        df_h2 = df_h2.astype(str)

        df['TIG_REGION'] = asmlib.call.val_per_hap(df, df_h1, df_h2, 'TIG_REGION')
        df['QUERY_STRAND'] = asmlib.call.val_per_hap(df, df_h1, df_h2, 'QUERY_STRAND')
        df['CI'] = asmlib.call.val_per_hap(df, df_h1, df_h2, 'CI')
        df['ALIGN_INDEX'] = asmlib.call.val_per_hap(df, df_h1, df_h2, 'ALIGN_INDEX')
        df['CLUSTER_MATCH'] = asmlib.call.val_per_hap(df, df_h1, df_h2, 'CLUSTER_MATCH')
        df['CALL_SOURCE'] = asmlib.call.val_per_hap(df, df_h1, df_h2, 'CALL_SOURCE')

        # Set inversion columns
        if wildcards.vartype_svtype == 'sv_inv':
            del(df['RGN_REF_DISC'])
            del(df['RGN_TIG_DISC'])
            del(df['FLAG_ID'])
            del(df['FLAG_TYPE'])

            df['RGN_REF_INNER'] = asmlib.call.val_per_hap(df, df_h1, df_h2, 'RGN_REF_INNER')
            df['RGN_TIG_INNER'] = asmlib.call.val_per_hap(df, df_h1, df_h2, 'RGN_TIG_INNER')

        # Load mapped regions
        map_tree_h1 = collections.defaultdict(intervaltree.IntervalTree)
        map_tree_h2 = collections.defaultdict(intervaltree.IntervalTree)

        df_map_h1 = pd.read_csv(input.callable_h1, sep='\t')
        df_map_h2 = pd.read_csv(input.callable_h2, sep='\t')

        for index, row in df_map_h1.iterrows():
            map_tree_h1[row['#CHROM']][row['POS']:row['END']] = True

        for index, row in df_map_h2.iterrows():
            map_tree_h2[row['#CHROM']][row['POS']:row['END']] = True

        # Get genotypes setting no-call for non-mappable regions
        df['GT_H1'] = df.apply(asmlib.call.get_gt, hap='h1', map_tree=map_tree_h1, axis=1)
        df['GT_H2'] = df.apply(asmlib.call.get_gt, hap='h2', map_tree=map_tree_h2, axis=1)

        df['GT'] = df.apply(lambda row: '{}|{}'.format(row['GT_H1'], row['GT_H2']), axis=1)

        if np.any(df['GT'].apply(lambda val: val == '0|0')):
            raise RuntimeError('Program bug: Found 0|0 genotypes after merging haplotypes')

        del df['GT_H1']
        del df['GT_H2']

        # Save BED
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# call_mappable_bed
#
# Make a table of mappable regions by merging aligned loci with loci covered by alignment-truncating events.
# "flank" parameter is an integer describing how far away records may be to merge (similar to the "bedtools merge"
# "slop" parameter). The flank is not added to the regions that are output.
rule call_mappable_bed:
    input:
        bed_align='results/{asm_name}/align/aligned_tig_{hap}.bed.gz',
        bed_lg_del='results/{asm_name}/lg_sv/sv_del_{hap}.bed.gz',
        bed_lg_ins='results/{asm_name}/lg_sv/sv_ins_{hap}.bed.gz',
        bed_lg_inv='results/{asm_name}/lg_sv/sv_inv_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/callable/callable_regions_{hap}_{flank}.bed.gz'
    run:

        # Get flank param
        try:
            flank = np.int32(wildcards.flank)

        except ValueError:
            raise RuntimeError('Flank parameter is not an integer: {flank}'.format(**wildcards))

        # Merge
        df = asmlib.util.region_merge(
            [
                input.bed_align,
                input.bed_lg_del,
                input.bed_lg_ins,
                input.bed_lg_inv
            ],
            pad=flank
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# call_correct_inter_inv
#
# Filter variants from inside inversions
rule call_integrate_sources:
    input:
        bed_cigar_insdel='temp/{asm_name}/cigar/pre_inv/svindel_insdel_{hap}.bed.gz',
        bed_cigar_snv='temp/{asm_name}/cigar/pre_inv/snv_snv_{hap}.bed.gz',
        bed_lg_ins='results/{asm_name}/lg_sv/sv_ins_{hap}.bed.gz',
        bed_lg_del='results/{asm_name}/lg_sv/sv_del_{hap}.bed.gz',
        bed_lg_inv='results/{asm_name}/lg_sv/sv_inv_{hap}.bed.gz',
        bed_inv='results/{asm_name}/inv_caller/sv_inv_{hap}.bed.gz'
    output:
        bed_ins=temp('temp/{asm_name}/bed/integrated/{hap}/svindel_ins.bed.gz'),
        bed_del=temp('temp/{asm_name}/bed/integrated/{hap}/svindel_del.bed.gz'),
        bed_inv=temp('temp/{asm_name}/bed/integrated/{hap}/sv_inv.bed.gz'),
        bed_snv=temp('temp/{asm_name}/bed/integrated/{hap}/snv_snv.bed.gz')
    params:
        min_inv=config.get('min_inv', 300),
        max_inv=config.get('max_inv', 2000000)
    run:

        # Set parameters
        if params.min_inv is not None and params.min_inv != 'unlimited':
            min_inv = int(params.min_inv)
        else:
            min_inv = None

        if params.max_inv is not None and params.max_inv != 'unlimited':
            max_inv = int(params.max_inv)
        else:
            max_inv = None

        # Read tig filter (if present)
        tig_filter_tree = None

        if 'tig_filter_pattern' in config:

            tig_filter_file = config['tig_filter_pattern'].format(**wildcards)

            if os.path.isfile(tig_filter_file):
                tig_filter_tree = collections.defaultdict(intervaltree.IntervalTree)
                df_filter = pd.read_csv(tig_filter_file, sep='\t', header=None, comment='#', usecols=(0, 1, 2))
                df_filter.columns = ['#CHROM', 'POS', 'END']

                for index, row in df_filter.iterrows():
                    tig_filter_tree[row['#CHROM']][row['POS']:row['END']] = True

        # Read INV calls
        df_inv = pd.concat(
            [
                pd.read_csv(input.bed_inv, sep='\t', low_memory=False),
                pd.read_csv(input.bed_lg_inv, sep='\t', low_memory=False)
            ],
            axis=0
        ).sort_values(
            ['#CHROM', 'POS']
        ).reset_index(drop=True)

        if min_inv is not None:
            df_inv = df_inv.loc[df_inv['SVLEN'] >= min_inv]

        if max_inv is not None:
            df_inv = df_inv.loc[df_inv['SVLEN'] <= max_inv]

        # Apply contig filter to INV
        df_inv = asmlib.call.filter_by_tig_tree(df_inv, tig_filter_tree)

        # Filter overlapping inversion calls
        inv_tree = collections.defaultdict(intervaltree.IntervalTree)
        inv_index_set = set()

        for index, row in df_inv.sort_values(['SVLEN', 'POS']).iterrows():
            if len(inv_tree[row['#CHROM']][row['POS']:row['END']]) == 0:
                inv_index_set.add(index)
                inv_tree[row['#CHROM']][row['POS']:row['END']] = row['ID']

        df_inv = df_inv.loc[inv_index_set].sort_values(['#CHROM', 'POS'])

        # Initialize filter with inversions
        filter_tree = collections.defaultdict(intervaltree.IntervalTree)

        for index, row in df_inv.iterrows():
            filter_tree[row['#CHROM']][row['POS']:row['END']] = row['ID']

        # Read large variants and filter by inversions
        df_lg_ins = pd.read_csv(input.bed_lg_ins, sep='\t', low_memory=False)
        df_lg_del = pd.read_csv(input.bed_lg_del, sep='\t', low_memory=False)

        if df_lg_ins.shape[0] > 0:
            df_lg_ins = df_lg_ins.loc[
                df_lg_ins.apply(
                    lambda row: len(filter_tree[row['#CHROM'][row['POS']:row['END']]]) == 0,
                    axis=1
                )
            ]

            # Apply contig filter to large INS
            df_lg_ins = asmlib.call.filter_by_tig_tree(df_lg_ins, tig_filter_tree)

        if df_lg_del.shape[0] > 0:
            df_lg_del = df_lg_del.loc[
                df_lg_del.apply(
                    lambda row: len(filter_tree[row['#CHROM'][row['POS']:row['END']]]) == 0,
                    axis=1
                )
            ]

            # Apply contig filter to large DEL
            df_lg_del = asmlib.call.filter_by_tig_tree(df_lg_del, tig_filter_tree)

            # Add large deletions to filter
            for index, row in df_lg_del.iterrows():
                filter_tree[row['#CHROM']][row['POS']:row['END']] = row['ID']

        # Read CIGAR calls
        df_cigar_insdel = pd.read_csv(input.bed_cigar_insdel, sep='\t', low_memory=False)
        df_snv = pd.read_csv(input.bed_cigar_snv, sep='\t', low_memory=False)

        # Check column conformance among INS/DEL callsets (required for merging)
        if list(df_cigar_insdel.columns) != list(df_lg_ins.columns):
            raise RuntimeError('Columns from CIGAR and large SV INS callsets do not match')

        if list(df_cigar_insdel.columns) != list(df_lg_del.columns):
            raise RuntimeError('Columns from CIGAR and large SV DEL callsets do not match')

        # Filter CIGAR calls
        if df_cigar_insdel.shape[0] > 0:
            df_cigar_insdel = df_cigar_insdel.loc[
                df_cigar_insdel.apply(
                    lambda row: len(filter_tree[row['#CHROM'][row['POS']:row['END']]]) == 0,
                    axis=1
                )
            ]

        if df_snv.shape[0] > 0:
            df_snv = df_snv.loc[
                df_snv.apply(
                    lambda row: len(filter_tree[row['#CHROM'][row['POS']:row['END']]]) == 0,
                    axis=1
                )
            ]

        # Apply contig filter to small variants
        df_cigar_insdel = asmlib.call.filter_by_tig_tree(df_cigar_insdel, tig_filter_tree)
        df_snv = asmlib.call.filter_by_tig_tree(df_snv, tig_filter_tree)

        # Merge insertion/deletion variants
        df_insdel = pd.concat(
            [
                df_lg_ins,
                df_lg_del,
                df_cigar_insdel
            ],
            axis=0
        ).sort_values(['#CHROM', 'POS'])

        # Write
        df_insdel.loc[
            df_insdel['SVTYPE'] == 'INS'
        ].to_csv(
            output.bed_ins, sep='\t', index=False, compression='gzip'
        )

        df_insdel.loc[
            df_insdel['SVTYPE'] == 'DEL'
        ].to_csv(
            output.bed_del, sep='\t', index=False, compression='gzip'
        )

        df_inv.to_csv(output.bed_inv, sep='\t', index=False, compression='gzip')
        df_snv.to_csv(output.bed_snv, sep='\t', index=False, compression='gzip')

# call_inv_bed
#
# Make inversion call BED.
rule call_inv_bed:
    input:
        bed='results/{asm_name}/inv_caller/sv_inv_{hap}.bed.gz'
    output:
        bed=temp('temp/{asm_name}/bed/pre_merge/sv_inv_{hap}.bed.gz'),
        bed_dropped='results/{asm_name}/bed/dropped/shortinv_sv_inv_{hap}.bed.gz'
    params:
        min_svlen=config.get('inv_min_svlen', 300)
    run:

        # Read inversions
        df = pd.read_csv(input.bed, sep='\t')

        # Filter
        df_drop = df.loc[df['SVLEN'] < params.min_svlen]
        df_drop.to_csv(output.bed_dropped, sep='\t', index=False, compression='gzip')

        df = df.loc[[index for index in df.index if index not in df_drop.index]]

        df.drop_duplicates('ID', inplace=True)

        # Write BED
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


#
# Call from CIGAR
#

# call_cigar_merge
#
# Merge discovery sets from each batch.
rule call_cigar_merge:
    input:
        bed_insdel=expand('temp/{{asm_name}}/cigar/batched/insdel_{{hap}}_{batch}.bed.gz', batch=range(CALL_CIGAR_BATCH_COUNT)),
        bed_snv=expand('temp/{{asm_name}}/cigar/batched/snv.bed_{{hap}}_{batch}.gz', batch=range(CALL_CIGAR_BATCH_COUNT))
    output:
        bed_insdel=temp('temp/{asm_name}/cigar/pre_inv/svindel_insdel_{hap}.bed.gz'),
        bed_snv=temp('temp/{asm_name}/cigar/pre_inv/snv_snv_{hap}.bed.gz')
    run:

        # Read, merge, sort, write
        pd.concat(
            [pd.read_csv(file_name, sep='\t') for file_name in input.bed_insdel],
            axis=0
        ).sort_values(
            ['#CHROM', 'POS']
        ).to_csv(
            output.bed_insdel, sep='\t', index=False, compression='gzip'
        )

        pd.concat(
            [pd.read_csv(file_name, sep='\t') for file_name in input.bed_snv],
            axis=0
        ).sort_values(
            ['#CHROM', 'POS']
        ).to_csv(
            output.bed_snv, sep='\t', index=False, compression='gzip'
        )


# call_cigar
#
# Call variants by alignment CIGAR parsing.
rule call_cigar:
    input:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz',
        tig_fa_name='results/{asm_name}/align/contigs_{hap}.fa.gz'
    output:
        bed_insdel=temp('temp/{asm_name}/cigar/batched/insdel_{hap}_{batch}.bed.gz'),
        bed_snv=temp('temp/{asm_name}/cigar/batched/snv.bed_{hap}_{batch}.gz')
    run:

        batch = int(wildcards.batch)

        # Read
        df_align = pd.read_csv(input.bed, sep='\t', dtype={'#CHROM': str})

        df_align = df_align.loc[df_align['INDEX'] % CALL_CIGAR_BATCH_COUNT == batch]

        # Call
        df_snv, df_insdel = asmlib.cigarcall.make_insdel_snv_calls(df_align, REF_FA, input.tig_fa_name, wildcards.hap)

        # Write
        df_insdel.to_csv(output.bed_insdel, sep='\t', index=False, compression='gzip')
        df_snv.to_csv(output.bed_snv, sep='\t', index=False, compression='gzip')
