"""
Call variants from aligned contigs.
"""

localrules: call_merge_batch_table

FILTER_REASON = {
    'TIG_FILTER': 'TIG_FILTER: Contig filter region (tig_filter_pattern)',
    'COMPOUND': 'COMPOUND: Inside larger variant',
    'COMPOUND_INV': 'COMPOUND_INV: Inside a larger inversion',
    'INV_MIN': 'INV_MIN: Less than min INV size ({})',
    'INV_MAX': 'INV_MAX: Exceeds max INV size ({})'
}

MERGE_BATCH_COUNT = 20


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
        df_ins = pd.read_csv(input.bed_ins, sep='\t', low_memory=False, keep_default_na=False)
        df_del = pd.read_csv(input.bed_del, sep='\t', low_memory=False, keep_default_na=False)
        df_inv = pd.read_csv(input.bed_inv, sep='\t', low_memory=False, keep_default_na=False)

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
                SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')

            del(df['SEQ'])

            # Arrange columns
            df = svpoplib.variant.order_variant_columns(
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
        df = pd.read_csv(input.bed_snv, sep='\t', low_memory=False, keep_default_na=False)
        df.to_csv(output.bed_snv_snv, sep='\t', index=False, compression='gzip')


#
# Merge haplotypes
#

# pg_variant_bed
#
# Concatenate variant BED files from batched merges.
rule call_merge_haplotypes:
    input:
        bed_batch=lambda wildcards: [
            'temp/{asm_name}/bed/bychrom/{vartype_svtype}/{batch}.bed.gz'.format(
                asm_name=wildcards.asm_name, vartype_svtype=wildcards.vartype_svtype, batch=batch
            ) for batch in range(MERGE_BATCH_COUNT)
        ]
    output:
        bed=temp('temp/{asm_name}/bed/merged/{vartype_svtype}.bed.gz')
    run:

        df = pd.concat(
            [pd.read_csv(file_name, sep='\t') for file_name in input.bed_batch],
            axis=0
        ).sort_values(
            ['#CHROM', 'POS', 'END', 'ID']
        ).to_csv(
            output.bed, sep='\t', index=False, compression='gzip'
        )


# call_merge_haplotypes_chrom
#
# Merge by batches.
rule call_merge_haplotypes_batch:
    input:
        tsv='data/ref/merge_batch.tsv.gz',
        bed_var_h1='results/{asm_name}/bed/pre_merge/h1/{vartype_svtype}.bed.gz',
        bed_var_h2='results/{asm_name}/bed/pre_merge/h2/{vartype_svtype}.bed.gz',
        callable_h1='results/{asm_name}/callable/callable_regions_h1_500.bed.gz',
        callable_h2='results/{asm_name}/callable/callable_regions_h2_500.bed.gz'
    output:
        bed='temp/{asm_name}/bed/bychrom/{vartype_svtype}/{batch}.bed.gz'
    params:
        merge_threads=lambda wildcards: int(get_config(wildcards, 'merge_threads', 12))
    run:

        # Read batch table
        df_batch = pd.read_csv(input.tsv, sep='\t')
        df_batch = df_batch.loc[df_batch['BATCH'] == int(wildcards.batch)]

        # Get variant type
        var_svtype_list = wildcards.vartype_svtype.split('_')

        if len(var_svtype_list) != 2:
            raise RuntimeError('Wildcard "vartype_svtype" must be two elements separated by an underscore: {}'.format(wildcards.var_svtype))

        # Get configured merge definition
        config_def = pavlib.call.get_merge_params(wildcards, get_config(wildcards))

        print('Merging with def: ' + config_def)
        sys.stdout.flush()

        # Merge
        df_list = list()

        for chrom in df_batch['CHROM']:
            df_list.append(
                pavlib.call.merge_haplotypes(
                    input.bed_var_h1, input.bed_var_h2,
                    input.callable_h1, input.callable_h2,
                    config_def,
                    threads=params.merge_threads,
                    chrom=chrom,
                    is_inv=var_svtype_list[1] == 'inv'
                )
            )

        # Concat and save
        if len(df_list) > 0:
            pd.concat(
                df_list, axis=0
            ).to_csv(
                output.bed, sep='\t', index=False, compression='gzip'
            )
        else:
            with open(output.bed, 'wt') as out_file:
                pass



# call_merge_batch_table
#
# Create a table of merge batch assignments
rule call_merge_batch_table:
    input:
        tsv='data/ref/contig_info.tsv.gz'
    output:
        tsv='data/ref/merge_batch.tsv.gz'
    run:

        # Read and sort
        df = pd.read_csv(
            'data/ref/contig_info.tsv.gz', sep='\t'
        ).sort_values(
            'LEN', ascending=False
        ).set_index(
            'CHROM'
        )[['LEN']]

        df['BATCH'] = -1

        # Get a list of assignments for each batch
        list_chr = collections.defaultdict(list)
        list_size = collections.Counter()

        def get_smallest():
            """
            Get the next smallest bin.
            """

            min_index = 0

            for i in range(MERGE_BATCH_COUNT):

                if list_size[i] == 0:
                    return i

                if list_size[i] < list_size[min_index]:
                    min_index = i

            return min_index

        for chrom in df.index:
            i = get_smallest()
            df.loc[chrom, 'BATCH'] = i
            list_size[i] += df.loc[chrom, 'LEN']

        # Check
        if np.any(df['BATCH'] < 0):
            raise RuntimeError('Failed to assign all reference contigs to batches (PROGRAM BUG)')

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')


#
# Integrate variant calls
#

# call_mappable_bed
#
# Make a table of mappable regions by merging aligned loci with loci covered by alignment-truncating events.
# "flank" parameter is an integer describing how far away records may be to merge (similar to the "bedtools merge"
# "slop" parameter). The flank is not added to the regions that are output.
rule call_mappable_bed:
    input:
        bed_align='results/{asm_name}/align/aligned_tig_{hap}.bed.gz',
        bed_lg_del='temp/{asm_name}/lg_sv/sv_del_{hap}.bed.gz',
        bed_lg_ins='temp/{asm_name}/lg_sv/sv_ins_{hap}.bed.gz',
        bed_lg_inv='temp/{asm_name}/lg_sv/sv_inv_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/callable/callable_regions_{hap}_{flank}.bed.gz'
    run:

        # Get flank param
        try:
            flank = np.int32(wildcards.flank)

        except ValueError:
            raise RuntimeError('Flank parameter is not an integer: {flank}'.format(**wildcards))

        # Merge
        df = pavlib.util.region_merge(
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
        bed_lg_ins='temp/{asm_name}/lg_sv/sv_ins_{hap}.bed.gz',
        bed_lg_del='temp/{asm_name}/lg_sv/sv_del_{hap}.bed.gz',
        bed_lg_inv='temp/{asm_name}/lg_sv/sv_inv_{hap}.bed.gz',
        bed_inv='temp/{asm_name}/inv_caller/sv_inv_{hap}.bed.gz'
    output:
        bed_ins='results/{asm_name}/bed/pre_merge/{hap}/svindel_ins.bed.gz',
        bed_del='results/{asm_name}/bed/pre_merge/{hap}/svindel_del.bed.gz',
        bed_inv='results/{asm_name}/bed/pre_merge/{hap}/sv_inv.bed.gz',
        bed_snv='results/{asm_name}/bed/pre_merge/{hap}/snv_snv.bed.gz',
        bed_inv_drp='results/{asm_name}/bed/dropped/sv_inv_{hap}.bed.gz',
        bed_insdel_drp='results/{asm_name}/bed/dropped/svindel_insdel_{hap}.bed.gz',
        bed_snv_drp='results/{asm_name}/bed/dropped/snv_snv_{hap}.bed.gz',
    params:
        inv_min=lambda wildcards: get_config(wildcards, 'inv_min', 0),
        inv_max=lambda wildcards: get_config(wildcards, 'inv_max', 1e10),
        inv_inner=lambda wildcards: pavlib.util.as_bool(get_config(wildcards, 'inv_inner', False)),
        redundant_callset=lambda wildcards: pavlib.util.as_bool(get_config(wildcards, 'redundant_callset', False))
    run:

        local_config = get_config(wildcards)

        # Init dropped variant lists
        inv_drp_list = list()
        insdel_drp_list = list()
        snv_drp_list = list()

        # Set parameters
        if params.inv_min is not None and params.inv_min != 'unlimited':
            inv_min = int(params.inv_min)
        else:
            inv_min = None

        if params.inv_max is not None and params.inv_max != 'unlimited':
            inv_max = int(params.inv_max)
        else:
            inv_max = None

        # Read tig filter (if present)
        tig_filter_tree = None

        if 'tig_filter_pattern' in local_config:

            tig_filter_file = local_config['tig_filter_pattern'].format(**wildcards)

            if os.path.isfile(tig_filter_file):
                tig_filter_tree = collections.defaultdict(intervaltree.IntervalTree)
                df_filter = pd.read_csv(tig_filter_file, sep='\t', header=None, comment='#', usecols=(0, 1, 2))
                df_filter.columns = ['#CHROM', 'POS', 'END']

                for index, row in df_filter.iterrows():
                    tig_filter_tree[row['#CHROM']][row['POS']:row['END']] = True


        ### INV ###

        # Read INV calls
        df_inv = pd.concat(
            [
                pd.read_csv(input.bed_inv, sep='\t', low_memory=False, keep_default_na=False),
                pd.read_csv(input.bed_lg_inv, sep='\t', low_memory=False, keep_default_na=False)
            ],
            axis=0
        ).sort_values(
            ['#CHROM', 'POS', 'END', 'ID']
        ).reset_index(drop=True)

        df_inv['ID'] = svpoplib.variant.version_id(df_inv['ID'])

        # INV: Filter by contig blacklisted regions
        df_inv, df_inv_drp = pavlib.call.filter_by_tig_tree(df_inv, tig_filter_tree)

        df_inv_drp['REASON'] = FILTER_REASON['TIG_FILTER']
        inv_drp_list.append(df_inv_drp)

        # INV: Filter by size
        if inv_min is not None:
            inv_drp_list.append(df_inv.loc[df_inv['SVLEN'] < inv_min].copy())
            inv_drp_list[-1]['REASON'] = FILTER_REASON['INV_MIN'].format(inv_min)
            df_inv = df_inv.loc[df_inv['SVLEN'] >= inv_min].copy()

        if inv_max is not None:
            inv_drp_list.append(df_inv.loc[df_inv['SVLEN'] > inv_max].copy())
            inv_drp_list[-1]['REASON'] = FILTER_REASON['INV_MAX'].format(inv_max)
            df_inv = df_inv.loc[df_inv['SVLEN'] <= inv_max].copy()

        # INV: Filter compound INVs
        if not params.inv_inner:

            # Build filter tree based on INV loci
            inv_tree = collections.defaultdict(intervaltree.IntervalTree)
            inv_index_set = set()
            inv_drp_index_set = set()

            for index, row in df_inv.sort_values(['SVLEN', 'POS']).iterrows():
                if len(inv_tree[row['#CHROM']][row['POS']:row['END']]) == 0:
                    inv_index_set.add(index)
                    inv_tree[row['#CHROM']][row['POS']:row['END']] = row['ID']
                else:
                    inv_drp_index_set.add(index)

            if inv_drp_index_set:
                inv_drp_list.append(df_inv.loc[inv_drp_index_set].sort_values(['#CHROM', 'POS']))
                inv_drp_list[-1]['REASON'] = FILTER_REASON['COMPOUND_INV']

            df_inv = df_inv.loc[inv_index_set].sort_values(['#CHROM', 'POS']).copy()


        ### Init filter tree ###

        # Create large variant filter tree (for filtering small variants inside larger ones)
        filter_tree = collections.defaultdict(intervaltree.IntervalTree)

        # Initialize large variant filter with INV locations
        if not params.inv_inner:
            for index, row in df_inv.iterrows():
                filter_tree[row['#CHROM']][row['POS']:row['END']] = row['TIG_REGION'].split(':', 1)[0]


        ### INS/DEL - Large ###

        # Large: Read
        df_lg_ins = pd.read_csv(input.bed_lg_ins, sep='\t', low_memory=False, keep_default_na=False)
        df_lg_del = pd.read_csv(input.bed_lg_del, sep='\t', low_memory=False, keep_default_na=False)

        # Large: Filter by tig regions
        df_lg_ins, df_ins_drp = pavlib.call.filter_by_tig_tree(df_lg_ins, tig_filter_tree)
        df_lg_del, df_del_drp = pavlib.call.filter_by_tig_tree(df_lg_del, tig_filter_tree)

        df_ins_drp['REASON'] = FILTER_REASON['TIG_FILTER']
        df_del_drp['REASON'] = FILTER_REASON['TIG_FILTER']

        insdel_drp_list.append(df_ins_drp)
        insdel_drp_list.append(df_del_drp)

        # Large: Compound filter
        df_lg_ins, df_ins_drp = pavlib.call.filter_by_ref_tree(df_lg_ins, filter_tree, match_tig=params.redundant_callset)
        df_lg_ins, df_del_drp = pavlib.call.filter_by_ref_tree(df_lg_ins, filter_tree, match_tig=params.redundant_callset)

        df_ins_drp['REASON'] = FILTER_REASON['COMPOUND']
        df_del_drp['REASON'] = FILTER_REASON['COMPOUND']

        insdel_drp_list.append(df_ins_drp)
        insdel_drp_list.append(df_del_drp)

        # Large: Add DEL to filter regions
        for index, row in df_lg_del.iterrows():
            filter_tree[row['#CHROM']][row['POS']:row['END']] = row['TIG_REGION'].split(':', 1)[0]


        ### CIGAR calls ###

        # CIGAR: Read
        df_insdel = pd.read_csv(input.bed_cigar_insdel, sep='\t', low_memory=False, keep_default_na=False)
        df_snv = pd.read_csv(input.bed_cigar_snv, sep='\t', low_memory=False, keep_default_na=False)

        # CIGAR: Check column conformance among INS/DEL callsets (required for merging)
        if list(df_insdel.columns) != list(df_lg_ins.columns):
            raise RuntimeError('Columns from CIGAR and large SV INS callsets do not match')

        if list(df_insdel.columns) != list(df_lg_del.columns):
            raise RuntimeError('Columns from CIGAR and large SV DEL callsets do not match')

        # CIGAR: Filter by tig regions
        df_insdel, df_insdel_drp = pavlib.call.filter_by_tig_tree(df_insdel, tig_filter_tree)
        df_snv, df_snv_drp = pavlib.call.filter_by_tig_tree(df_snv, tig_filter_tree)

        df_insdel_drp['REASON'] = FILTER_REASON['TIG_FILTER']
        df_snv_drp['REASON'] = FILTER_REASON['TIG_FILTER']

        df_insdel_drp.append(df_insdel_drp)
        df_snv_drp.append(df_snv_drp)

        # CIGAR: Compound filter
        df_insdel, df_insdel_drp = pavlib.call.filter_by_ref_tree(df_insdel, filter_tree, match_tig=params.redundant_callset)
        df_snv, df_snv_drp = pavlib.call.filter_by_ref_tree(df_snv, filter_tree, match_tig=params.redundant_callset)

        df_insdel_drp['REASON'] = FILTER_REASON['COMPOUND']
        df_snv_drp['REASON'] = FILTER_REASON['COMPOUND']

        insdel_drp_list.append(df_insdel_drp)
        snv_drp_list.append(df_snv_drp)


        ### Concat ###

        # Concat: Merge insertion/deletion variants
        df_insdel = pd.concat(
            [
                df_lg_ins,
                df_lg_del,
                df_insdel
            ],
            axis=0
        ).sort_values(['#CHROM', 'POS'])

        # Concat: De-duplicate IDs (can occur with redundant callsets)
        df_insdel['ID'] = svpoplib.variant.version_id(df_insdel['ID'])
        df_inv['ID'] = svpoplib.variant.version_id(df_inv['ID'])
        df_snv['ID'] = svpoplib.variant.version_id(df_snv['ID'])


        ### Write ###

        # Write: dropped variants
        pd.concat(inv_drp_list).to_csv(output.bed_inv_drp, sep='\t', index=False, compression='gzip')
        pd.concat(insdel_drp_list).to_csv(output.bed_insdel_drp, sep='\t', index=False, compression='gzip')
        pd.concat(snv_drp_list).to_csv(output.bed_snv_drp, sep='\t', index=False, compression='gzip')

        # Write: Passed variants
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


#
# Call from CIGAR
#

# call_cigar_merge
#
# Merge discovery sets from each batch.
rule call_cigar_merge:
    input:
        bed_insdel=expand('temp/{{asm_name}}/cigar/batched/insdel_{{hap}}_{batch}.bed.gz', batch=range(pavlib.cigarcall.CALL_CIGAR_BATCH_COUNT)),
        bed_snv=expand('temp/{{asm_name}}/cigar/batched/snv.bed_{{hap}}_{batch}.gz', batch=range(pavlib.cigarcall.CALL_CIGAR_BATCH_COUNT))
    output:
        bed_insdel=temp('temp/{asm_name}/cigar/pre_inv/svindel_insdel_{hap}.bed.gz'),
        bed_snv=temp('temp/{asm_name}/cigar/pre_inv/snv_snv_{hap}.bed.gz')
    run:

        # INS/DEL
        df_insdel = pd.concat(
            [pd.read_csv(file_name, sep='\t', keep_default_na=False) for file_name in input.bed_insdel],
            axis=0
        )

        df_insdel['ID'] = svpoplib.variant.version_id(df_insdel['ID'])

        df_insdel.sort_values(
            ['#CHROM', 'POS', 'END', 'ID']
        ).to_csv(
            output.bed_insdel, sep='\t', index=False, compression='gzip'
        )

        # SNV
        df_snv = pd.concat(
            [pd.read_csv(file_name, sep='\t', keep_default_na=False) for file_name in input.bed_snv],
            axis=0
        )

        df_snv['ID'] = svpoplib.variant.version_id(df_snv['ID'])

        df_snv.sort_values(
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
        tig_fa_name='temp/{asm_name}/align/contigs_{hap}.fa.gz'
    output:
        bed_insdel=temp('temp/{asm_name}/cigar/batched/insdel_{hap}_{batch}.bed.gz'),
        bed_snv=temp('temp/{asm_name}/cigar/batched/snv.bed_{hap}_{batch}.gz')
    run:

        batch = int(wildcards.batch)

        os.makedirs(os.path.dirname(output.bed_insdel), exist_ok=True)  # Random crashes with "FileNotFoundError", Snakemake not creating output directory?

        # Read
        df_align = pd.read_csv(input.bed, sep='\t', dtype={'#CHROM': str}, keep_default_na=False)

        df_align = df_align.loc[df_align['CALL_BATCH'] == batch]

        # Call
        df_snv, df_insdel = pavlib.cigarcall.make_insdel_snv_calls(df_align, REF_FA, input.tig_fa_name, wildcards.hap)

        # Write
        df_insdel.to_csv(output.bed_insdel, sep='\t', index=False, compression='gzip')
        df_snv.to_csv(output.bed_snv, sep='\t', index=False, compression='gzip')
