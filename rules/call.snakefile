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
        bed_insdel='temp/{asm_name}/bed/svindel_insdel.bed.gz',
        bed_inv='temp/{asm_name}/bed/sv_inv.bed.gz',
        bed_snv='temp/{asm_name}/bed/snv_snv.bed.gz'
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
        df_insdel = pd.read_csv(input.bed_insdel, sep='\t', low_memory=False)
        df_inv = pd.read_csv(input.bed_inv, sep='\t', low_memory=False)

        for vartype, svtype in [('sv', 'ins'), ('sv', 'del'), ('sv', 'inv'), ('indel', 'ins'), ('indel', 'del')]:

            # Subset
            if vartype == 'sv':
                if svtype == 'inv':
                    df = df_inv
                else:
                    df = df_insdel.loc[df_insdel['SVLEN'] >= 50]

            elif vartype == 'indel':
                df = df_insdel.loc[df_insdel['SVLEN'] < 50]

            else:
                raise RuntimeError('Program Bug: Unknown variant type: {}'.format(vartype))

            df = df.loc[df['SVTYPE'] == svtype.upper()]

            # Get output file names
            bed_file_name = output[f'bed_{vartype}_{svtype}']
            fa_file_name = output[f'fa_{vartype}_{svtype}']

            # Write FASTA
            with Bio.bgzf.BgzfWriter(fa_file_name, 'wb') as out_file:
                SeqIO.write(analib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')

            del(df['SEQ'])

            # Write BED
            df.to_csv(bed_file_name, sep='\t', index=False, compression='gzip')

        # Process SNVs
        df = pd.read_csv(input.bed_snv_snv, sep='\t', low_memory=False)
        df.to_csv(output.bed_snv_snv, sep='\t', index=False, compression='gzip')


# pg_variant_bed
#
# Make variant BED and FASTA. Write all variant calls regardless of consensus loci.
rule call_merge_haplotypes:
    input:
        bed_var_h1='temp/{asm_name}/bed/h1/{vartype_svtype}.bed.gz',
        bed__var_h2='temp/{asm_name}/bed/h2/{vartype_svtype}.bed.gz',
        bed_align_h1='results/{asm_name}/align/aligned_tig_h1.bed.gz',
        bed_align_h2='results/{asm_name}/align/aligned_tig_h2.bed.gz',
        bed_lg_del_h1='results/{asm_name}/lg_sv/sv_del_h1.bed.gz',
        bed_lg_del_h2='results/{asm_name}/lg_sv/sv_del_h2.bed.gz'
    output:
        bed=temp('temp/{asm_name}/bed/{vartype_svtype}.bed.gz'),
        fa=temp('temp/{asm_name}/fa/{vartype_svtype}.fa.gz')
    params:
        ro_min=float(config.get('ro_min', 0.5)),
        offset_max=int(config.get('offset_max', 200)),
        merge_threads=int(config.get('merge_threads', 12)),
        merge_mem=config.get('merge_mem', '2G')
    run:

        # Get configured merge definition
        if wildcards.vartype == 'snv':
            config_def = 'nrid'
        else:
            config_def = 'nr:szro={}:offset={}:roor'.format(int(ro_min * 100), offset_max)

        print('Merging with def: ' + config_def)
        sys.stdout.flush()

        # Merge
        df = analib.svmerge.merge_variants(
            bed_list=[input.bed1, input.bed2],
            sample_names=['h1', 'h2'],
            strategy=config_def,
            threads=wildcards.merge_threads
        )

        # Restructure columns
        del(df['DISC_CLASS'])

        df.columns = [re.sub('^MERGE_', 'HAP_', val) for val in df.columns]

        del(df['HAP_SRC'])
        del(df['HAP_SRC_ID'])
        del(df['HAP_AC'])
        del(df['HAP_AF'])

        df.columns = ['HAP_SRC' if val == 'HAP_SAMPLES' else val for val in df.columns]

        # Restructure columns
        # TODO: Add alternate QUERY_ID:POS-END and alternate strand for h2 in homozygous calls. Add alternate CALL_SOURCE?

        # Load mapped regions
        map_tree_h1= collections.defaultdict(intervaltree.IntervalTree)
        map_tree_h2= collections.defaultdict(intervaltree.IntervalTree)

        df_map_h1 = pd.read_csv(input.bed_align1, sep='\t')
        df_map_h2 = pd.read_csv(input.bed_align2, sep='\t')

        for index, row in df_map_h1.iterrows():
            map_tree_h1[row['#CHROM']][row['POS']:row['END']] = True

        for index, row in df_map_h2.iterrows():
            map_tree_h2[row['#CHROM']][row['POS']:row['END']] = True

        # Add large SVs to map trees
        df_lg_del_h1 = pd.read_csv(input.bed_lg_del_h1, sep='\t', low_memory=False)

        for index, row in df_lg_del_h1.iterrows():
            map_tree_h1[row['#CHROM']][row['POS']:row['END']] = True

        df_lg_del_h2 = pd.read_csv(input.bed_lg_del_h2, sep='\t', low_memory=False)

        for index, row in df_lg_del_h2.iterrows():
            map_tree_h2[row['#CHROM']][row['POS']:row['END']] = True


        # Get genotypes setting no-call for non-mappable regions
        df['GT_H1'] = df.apply(asmlib.call.get_gt, hap='h1', map_tree=map_tree_h1, axis=1)
        df['GT_H2'] = df.apply(asmlib.call.get_gt, hap='h2', map_tree=map_tree_h2, axis=1)

        df['GT'] = df.apply(lambda row: '{}|{}'.format(row['GT_H1'], row['GT_H2']), axis=1)

        # Save BED
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
        bed_insdel=temp('temp/{asm_name}/bed/{hap}/svindel_insdel.bed.gz'),
        bed_inv=temp('temp/{asm_name}/bed/{hap}/sv_inv.bed.gz'),
        bed_snv=temp('temp/{asm_name}/bed/{hap}/snv_snv.bed.gz')
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


        # Read INV calls
        df_inv = pd.concat(
            [
                pd.read_csv(input.bed_inv, sep='\t', low_memory=False),
                pd.read_csv(input.bed_lg_inv, sep='\t', low_memory=False)
            ],
            axis=0
        ).sort_values(['#CHROM', 'POS'])

        if min_inv is not None:
            df_inv = df_inv.loc[df_inv['SVLEN'] >= min_inv]

        if max_inv is not None:
            df_inv = df_inv.loc[df_inv['SVLEN'] <= max_inv]

        # Initialize filter with inversions
        filter_tree = collections.defaultdict(intervaltree.IntervalTree)

        for index, row in df_inv.iterrows():
            filter_tree[row['#CHROM']][row['POS']:row['END']] = row['ID']

        # Read large variants and filter by inversions
        df_lg_ins = pd.read_csv(input.bed_lg_ins, sep='\t', low_memory=False)
        df_lg_del = pd.read_csv(input.bed_lg_del, sep='\t', low_memory=False)

        df_lg_ins = df_lg_ins.loc[
            df_lg_ins.apply(
                lambda row: len(filter_tree[row['#CHROM'][row['POS']:row['END']]]) == 0,
                axis=1
            )
        ]

        df_lg_del = df_lg_del.loc[
            df_lg_del.apply(
                lambda row: len(filter_tree[row['#CHROM'][row['POS']:row['END']]]) == 0,
                axis=1
            )
        ]

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
        df_cigar_insdel = df_cigar_insdel.loc[
            df_cigar_insdel.apply(
                lambda row: len(filter_tree[row['#CHROM'][row['POS']:row['END']]]) == 0,
                axis=1
            )
        ]

        df_snv = df_snv.loc[
            df_snv.apply(
                lambda row: len(filter_tree[row['#CHROM'][row['POS']:row['END']]]) == 0,
                axis=1
            )
        ]

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
        df_insdel.to_csv(output.bed_insdel, sep='\t', index=False, compression='gzip')
        df_inv.to_csv(output.bed_inv, sep='\t', index=False, compression='gzip')
        df_snv.to_csv(output.bed_snv, sep='\t', index=False, compression='gzip')



# rule call_integrate_sources:
#     input:
#         bed='temp/{asm_name}/bed/pre_merge/pre_inv_correction/{vartype}_{svtype}_{hap}.bed.gz',
#         bed_inv='temp/{asm_name}/bed/pre_merge/sv_inv_{hap}.bed.gz',
#         bed_lg_inv='results/{asm_name}/lg_sv/sv_inv_{hap}.bed.gz'
#     output:
#         bed=temp('temp/{asm_name}/bed/pre_merge/{vartype}_{svtype}_{hap}.bed.gz'),
#         bed_dropped='results/{asm_name}/bed/dropped/interinv_{vartype}_{svtype}_{hap}.bed.gz'
#     wildcard_constraints:
#         vartype='sv|indel|snv',
#         svtype='ins|del|snv'
#     run:
#
#         # Read
#         df = pd.read_csv(input.bed, sep='\t')
#
#         df_inv = pd.concat(
#             [
#                 pd.read_csv(input.bed_inv, sep='\t', usecols=['#CHROM', 'POS', 'END']),
#                 pd.read_csv(input.bed_lg_inv, sep='\t', usecols=['#CHROM', 'POS', 'END'])
#             ],
#             axis=0
#         )
#
#         # Build tree
#         invtree = collections.defaultdict(intervaltree.IntervalTree)
#
#         for index, row in pd.read_csv(input.bed_inv, sep='\t').iterrows():
#             invtree[row['#CHROM']][row['POS']:row['END']] = (row['ID'])
#
#         # Filter
#         filter_list = df.apply(
#             lambda row: len(invtree[row['#CHROM']][row['POS']:row['END']]) > 0, axis=1
#         )
#
#         # Write dropped
#         df.loc[filter_list].to_csv(output.bed_dropped, sep='\t', index=False, compression='gzip')
#
#         # Write
#         df.loc[~ filter_list].to_csv(output.bed, sep='\t', index=False, compression='gzip')

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
        df_align = pd.read_csv(input.bed, sep='\t')

        df_align = df_align.loc[df_align['INDEX'] % CALL_CIGAR_BATCH_COUNT == batch]

        # Call
        df_snv, df_insdel = asmlib.cigarcall.make_insdel_snv_calls(df_align, REF_FA, input.tig_fa_name, wildcards.hap)

        # Write
        df_insdel.to_csv(output.bed_insdel, sep='\t', index=False, compression='gzip')
        df_snv.to_csv(output.bed_snv, sep='\t', index=False, compression='gzip')
