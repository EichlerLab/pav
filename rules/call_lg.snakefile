"""
Call alignment-truncating events (large SVs).
"""

# call_merge_lg
#
# Merge variant calls from large SVs.
rule call_merge_lg:
    input:
        bed=expand(
            'temp/{{asm_name}}/lg_sv/sv_{{svtype}}_{{hap}}_{batch}.bed.gz',
            batch=range(config.get('lg_batch_count', 10))
        )
    output:
        bed='results/{asm_name}/lg_sv/sv_{svtype}_{hap}.bed.gz'
    run:

        pd.concat(
            [pd.read_csv(file_name, sep='\t') for file_name in input.bed],
            axis=1
        ).T.to_csv(
            output.bed, sep='\t', index=False, compression='gzip'
        )


# call_lg_discover
#
# Call alignment-truncating SVs.
rule call_lg_discover:
    input:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz',
        tsv_group='temp/{asm_name}/lg_sv/batch_{hap}.tsv.gz',
        fa='results/{asm_name}/align/contigs_{hap}.fa.gz',
        fai='results/{asm_name}/align/contigs_{hap}.fa.gz.fai',
        bed_n='data/ref/n_gap.bed.gz'
    output:
        bed_ins=temp('temp/{asm_name}/lg_sv/sv_ins_{hap}_{batch}.bed.gz'),
        bed_del=temp('temp/{asm_name}/lg_sv/sv_del_{hap}_{batch}.bed.gz'),
        bed_inv=temp('temp/{asm_name}/lg_sv/sv_inv_{hap}_{batch}.bed.gz')
    log:
        log='results/{asm_name}/lg_sv/log/sv_ins_{hap}_{batch}.log'
    params:
        k_size=config.get('inv_k_size', 31),
        inv_threads=config.get('inv_threads_lg', config.get('inv_threads', 12)),
        inv_mem=config.get('inv_mem', '4G')
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t')
        df_tig_fai = analib.ref.get_df_fai(input.fai)

        # Subset to alignment records in this batch
        df_group = pd.read_csv(input.tsv_group, sep='\t')
        df_group = df_group.loc[df_group['BATCH'] == int(wildcards.batch)]

        group_set = set(df_group[['CHROM', 'TIG']].apply(tuple, axis=1))

        df = df.loc[df.apply(lambda row: (row['#CHROM'], row['QUERY_ID']) in group_set, axis=1)]

        # Get trees of N bases
        n_tree = collections.defaultdict(intervaltree.IntervalTree)

        df_n = pd.read_csv(input.bed_n, sep='\t')

        for index, row in df_n.iterrows():
            n_tree[row['#CHROM']][row['POS']:row['END']] = True

        # Make density table output directory
        density_out_dir = 'results/{asm_name}/inv_caller/density_table'.format(**wildcards)
        os.makedirs(density_out_dir, exist_ok=True)

        # Get large events
        with open(log.log, 'w') as log_file:
            df_ins, df_del, df_inv = asmlib.lgsv.scan_for_events(
                df, df_tig_fai, wildcards.hap, REF_FA, input.fa,
                k_size=params.k_size,
                n_tree=n_tree,
                threads=params.inv_threads,
                log=log_file,
                density_out_dir=density_out_dir,
            )

        # Write
        df_ins.to_csv(output.bed_ins, sep='\t', index=False, compression='gzip')
        df_del.to_csv(output.bed_del, sep='\t', index=False, compression='gzip')
        df_inv.to_csv(output.bed_inv, sep='\t', index=False, compression='gzip')


# call_lg_split
#
# Split chromosome/tig records into batches for chromosome/tig pairs with multiple alignments.
rule call_lg_split:
    input:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz'
    output:
        tsv=temp('temp/{asm_name}/lg_sv/batch_{hap}.tsv.gz')
    params:
        batch_count=config.get('lg_batch_count', 10)
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t')

        # Get ref/tig pairs with multiple mappings
        tig_map_count = collections.Counter(df[['#CHROM', 'QUERY_ID']].apply(tuple, axis=1))

        df_group_list = list()

        index = 0

        for chrom, tig in [(chrom, tig_id) for (chrom, tig_id), count in tig_map_count.items() if count > 1]:
            df_group_list.append(pd.Series(
                [chrom, tig, index % params.batch_count],
                index=['CHROM', 'TIG', 'BATCH']
            ))

            index += 1

        # Merge (CHROM, TIG, BATCH)
        df_group = pd.concat(df_group_list, axis=1).T

        # Write
        df_group.to_csv(output.tsv, sep='\t', index=False, compression='gzip')
