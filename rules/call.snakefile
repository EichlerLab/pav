"""
Call variants from aligned contigs.
"""


#
# Definitions
#


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
        df = pd.read_csv(input.bed_snv, sep='\t', low_memory=False)
        df.to_csv(output.bed_snv_snv, sep='\t', index=False, compression='gzip')


#
# Merge haplotypes
#

# pg_variant_bed
#
# Make variant BED and FASTA. Write all variant calls regardless of consensus loci.
#
# If merging is not done by chromosome (default), then this rule reads from each haplotype and merges in one step
# (wildcards.merge_level = merged). If merging is done per chromosome (config['merge_by_chrom'] is True), then this
# rule calls itself recursively first by chromosome (wildcards.merge_level = bychrom) then by concatenating the merged
# chromosomes (wildcards.merge_level = merged). The code will know which step its on based on the wildcards and config.
rule call_merge_haplotypes:
    input:
        bed_var_h1='temp/{asm_name}/bed/integrated/h1/{vartype_svtype}.bed.gz',
        bed_var_h2='temp/{asm_name}/bed/integrated/h2/{vartype_svtype}.bed.gz',
        callable_h1='results/{asm_name}/callable/callable_regions_h1_500.bed.gz',
        callable_h2='results/{asm_name}/callable/callable_regions_h2_500.bed.gz',
        bed_chrom=lambda wildcards: [
            'temp/{asm_name}/bed/bychrom/{vartype_svtype}/{chrom}.bed.gz'.format(
                asm_name=wildcards.asm_name, vartype_svtype=wildcards.vartype_svtype, chrom=chrom
            ) for chrom in sorted(svpoplib.ref.get_df_fai(config['reference'] + '.fai').index)
        ] if pavlib.util.as_bool(config.get('merge_by_chrom', True)) else []
    output:
        bed=temp('temp/{asm_name}/bed/merged/{vartype_svtype}.bed.gz')
    params:
        ro_min=float(config.get('ro_min', 0.5)),
        offset_max=int(config.get('offset_max', 200)),
        merge_threads=int(config.get('merge_threads', 12)),
        merge_align=config.get('merge_align', None)
    run:

        if pavlib.util.as_bool(config.get('merge_by_chrom', True)):
            # Merge in one step

            # Get configured merge definition
            if wildcards.vartype_svtype == 'snv_snv':
                config_def = 'nr:refalt'
            else:
                config_def = 'nr:szro={}:offset={}'.format(int(params.ro_min * 100), params.offset_max)

                # Get alignment merging parameters
                if params.merge_align is not None:

                    if params.merge_align.lower() == 'true':
                        align_params = ':match=0.8,2,-1,-1,-0.25'

                    elif params.merge_align.lower() == 'false':
                        align_params = ''

                    else:
                        align_params = f';match={params.merge_align}'

                    config_def += align_params

            print('Merging with def: ' + config_def)
            sys.stdout.flush()

            # Merge
            df = pavlib.call.merge_haplotypes(
                input.bed_var_h1, input.bed_var_h2,
                input.callable_h1,input.callable_h2,
                config_def,
                threads=params.merge_threads,
                chrom=None
            )

            # Save BED
            df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

        else:
            # Concatenate merged chromosomes

            write_header = True

            with gzip.open(output.bed, 'wt') as out_file:
                for in_file_name in input.bed_chrom:

                    df_iter = pd.read_csv(
                        in_file_name,
                        sep='\t', iterator=True, chunksize=20000
                    )

                    for df in df_iter:
                        df.to_csv(
                            out_file,sep='\t', index=False, header=write_header
                        )

                        write_header = False


# call_merge_haplotypes_chrom
#
# Merge by chromosome. This rule is used if "merge_by_chrom" is True.
rule call_merge_haplotypes_chrom:
    input:
        bed_var_h1='temp/{asm_name}/bed/integrated/h1/{vartype_svtype}.bed.gz',
        bed_var_h2='temp/{asm_name}/bed/integrated/h2/{vartype_svtype}.bed.gz',
        callable_h1='results/{asm_name}/callable/callable_regions_h1_500.bed.gz',
        callable_h2='results/{asm_name}/callable/callable_regions_h2_500.bed.gz'
    output:
        bed='temp/{asm_name}/bed/bychrom/{vartype_svtype}/{chrom}.bed.gz'
    params:
        ro_min=float(config.get('ro_min', 0.5)),
        offset_max=int(config.get('offset_max', 200)),
        merge_threads=int(config.get('merge_threads', 12)),
        merge_align=config.get('merge_align', 'true')
    run:

        var_svtype_list = wildcards.vartype_svtype.split('_')

        if len(var_svtype_list) != 2:
            raise RuntimeError('Wildcard "vartype_svtype" must be two elements separated by an underscore: {}'.format(wildcards.var_svtype))

        # Get configured merge definition
        if wildcards.vartype_svtype == 'snv_snv':
            config_def = 'nr:refalt'
        else:
            config_def = 'nr:szro={}:offset={}'.format(int(params.ro_min * 100), params.offset_max)

            # Get alignment merging parameters
            if params.merge_align is not None:

                if params.merge_align.lower() == 'true':
                    align_params = ':match=0.8,2,-1,-4,-0.25,500000,9'

                elif params.merge_align.lower() == 'false':
                    align_params = ''

                else:
                    align_params = f':match={params.merge_align}'

                config_def += align_params

        print('Merging with def: ' + config_def)
        sys.stdout.flush()

        # Merge
        df = pavlib.call.merge_haplotypes(
            input.bed_var_h1, input.bed_var_h2,
            input.callable_h1, input.callable_h2,
            config_def,
            threads=params.merge_threads,
            chrom=wildcards.chrom,
            is_inv=var_svtype_list[1] == 'inv'
        )

        # Save BED
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


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
        bed_ins=temp('temp/{asm_name}/bed/integrated/{hap}/svindel_ins.bed.gz'),
        bed_del=temp('temp/{asm_name}/bed/integrated/{hap}/svindel_del.bed.gz'),
        bed_inv=temp('temp/{asm_name}/bed/integrated/{hap}/sv_inv.bed.gz'),
        bed_snv=temp('temp/{asm_name}/bed/integrated/{hap}/snv_snv.bed.gz'),
        bed_inv_dropped='results/{asm_name}/bed/dropped/sv_inv_{hap}.bed.gz'
    params:
        inv_min=config.get('inv_min', 300),
        inv_max=config.get('inv_max', 2000000),
        redundant_callset=pavlib.util.as_bool(config.get('redundant_callset', False))
    run:

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
            ['#CHROM', 'POS', 'END', 'ID']
        ).reset_index(drop=True)

        df_inv['ID'] = svpoplib.variant.version_id(df_inv['ID'])

        # Filter INV by size
        inv_dropped_list = list()

        if inv_min is not None:
            inv_dropped_list.append(df_inv.loc[df_inv['SVLEN'] < inv_min])
            inv_dropped_list[-1]['REASON'] = 'Less than min size'
            df_inv = df_inv.loc[df_inv['SVLEN'] >= inv_min]

        if inv_max is not None:
            inv_dropped_list.append(df_inv.loc[df_inv['SVLEN'] > inv_max])
            inv_dropped_list[-1]['REASON'] = 'Greater than max size'
            df_inv = df_inv.loc[df_inv['SVLEN'] <= inv_max]

        # Apply contig filter to INV
        df_inv = pavlib.call.filter_by_tig_tree(df_inv, tig_filter_tree)

        # Filter overlapping inversion calls
        inv_tree = collections.defaultdict(intervaltree.IntervalTree)
        inv_index_set = set()
        inv_dropped_index_set = set()

        for index, row in df_inv.sort_values(['SVLEN', 'POS']).iterrows():
            if len(inv_tree[row['#CHROM']][row['POS']:row['END']]) == 0:
                inv_index_set.add(index)
                inv_tree[row['#CHROM']][row['POS']:row['END']] = row['ID']
            else:
                inv_dropped_index_set.add(index)

        if inv_dropped_index_set:
            inv_dropped_list.append(df_inv.loc[inv_dropped_index_set].sort_values(['#CHROM', 'POS']))
            inv_dropped_list[-1]['REASON'] = 'Overlaps another inversion'

        df_inv = df_inv.loc[inv_index_set].sort_values(['#CHROM', 'POS'])

        # Write dropped inversions
        pd.concat(inv_dropped_list).to_csv(output.bed_inv_dropped, sep='\t', index=False, compression='gzip')

        # Initialize filter with inversions
        filter_tree = collections.defaultdict(intervaltree.IntervalTree)

        for index, row in df_inv.iterrows():
            filter_tree[row['#CHROM']][row['POS']:row['END']] = row['TIG_REGION'].split(':', 1)[0]

        # Read large variants and filter by inversions
        df_lg_ins = pd.read_csv(input.bed_lg_ins, sep='\t', low_memory=False)
        df_lg_del = pd.read_csv(input.bed_lg_del, sep='\t', low_memory=False)

        if df_lg_ins.shape[0] > 0:
            df_lg_ins = pavlib.call.filter_by_ref_tree(df_lg_ins, filter_tree, match_tig=params.redundant_callset)

            # Apply contig filter to large INS
            df_lg_ins = pavlib.call.filter_by_tig_tree(df_lg_ins, tig_filter_tree)

        if df_lg_del.shape[0] > 0:
            df_lg_del = pavlib.call.filter_by_ref_tree(df_lg_del, filter_tree, match_tig=params.redundant_callset)

            # Apply contig filter to large DEL
            df_lg_del = pavlib.call.filter_by_tig_tree(df_lg_del, tig_filter_tree)

            # Add large deletions to filter
            for index, row in df_lg_del.iterrows():
                filter_tree[row['#CHROM']][row['POS']:row['END']] = row['TIG_REGION'].split(':', 1)[0]

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
            df_cigar_insdel = pavlib.call.filter_by_ref_tree(df_cigar_insdel, filter_tree, match_tig=params.redundant_callset)

        if df_snv.shape[0] > 0:
            df_snv = pavlib.call.filter_by_ref_tree(df_snv, filter_tree, match_tig=params.redundant_callset)

        # Apply contig filter to small variants
        df_cigar_insdel = pavlib.call.filter_by_tig_tree(df_cigar_insdel, tig_filter_tree)
        df_snv = pavlib.call.filter_by_tig_tree(df_snv, tig_filter_tree)

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
            [pd.read_csv(file_name, sep='\t') for file_name in input.bed_insdel],
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
            [pd.read_csv(file_name, sep='\t') for file_name in input.bed_snv],
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
        df_align = pd.read_csv(input.bed, sep='\t', dtype={'#CHROM': str})

        df_align = df_align.loc[df_align['CALL_BATCH'] == batch]

        # Call
        df_snv, df_insdel = pavlib.cigarcall.make_insdel_snv_calls(df_align, REF_FA, input.tig_fa_name, wildcards.hap)

        # Write
        df_insdel.to_csv(output.bed_insdel, sep='\t', index=False, compression='gzip')
        df_snv.to_csv(output.bed_snv, sep='\t', index=False, compression='gzip')
