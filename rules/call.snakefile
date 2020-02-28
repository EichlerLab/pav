"""
Call variants from aligned contigs.
"""

#
# INS/DEL/SNV
#

### PrintGaps Post-processing ###

# pg_variant_bed
#
# Make variant BED and FASTA. Write all variant calls regardless of consensus loci.
rule pg_variant_bed:
    input:
       bed1='results/{asm_name}/bed/pre_merge/{vartype}_{svtype}_h1.bed.gz',
       bed2='results/{asm_name}/bed/pre_merge/{vartype}_{svtype}_h2.bed.gz',
       con_bed='results/{asm_name}/align/depth_1/regions_h12.bed'
    output:
        bed=protected('results/{asm_name}/bed/{vartype}_{svtype}_h12.bed.gz'),
        fa=protected('results/{asm_name}/bed/fa/{vartype}_{svtype}_h12.fa.gz')
    run:

        # Get configured merge definition
        if wildcards.vartype == 'snv':
            config_def = 'nrid'
        else:
            config_def = 'nr:szro={}:offset={}:roor'.format(int(RO_MIN * 100), OFFSET_MAX)

        print('Merging with def: ' + config_def)
        sys.stdout.flush()

        # Merge
        df = analib.svmerge.merge_variants(
            bed_list=[input.bed1, input.bed2],
            sample_names=['h1', 'h2'],
            strategy=config_def,
            threads=6
        )

        # Restructure columns
        del(df['DISC_CLASS'])

        df.columns = [re.sub('^MERGE_', 'HAP_', val) for val in df.columns]

        # Load consensus regions
        consensus_tree = collections.defaultdict(intervaltree.IntervalTree)

        df_con = pd.read_csv(input.con_bed, sep='\t', header=None, names=('#CHROM', 'POS', 'END'))

        for index, row in df.iterrows():
            consensus_tree[row['#CHROM']][row['POS']:row['END']] = True

        # Define a function to annotate consensus regions
        def con_match(row):
            for interval in consensus_tree[row['#CHROM']].overlap(row['POS'], row['END']):
                if (interval.begin <= row['POS']) and (interval.end >= row['END']):
                    return True

            return False

        # Annotate consensus regions
        df['HAP_CONSENSUS'] = df.apply(con_match, axis=1)

        # Save SEQ as a FASTA
        if 'SEQ' in df.columns:
            with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
                SeqIO.write(analib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')

            del(df['SEQ'])
        else:
            with open(output.fa, 'w') as out_file:
                pass

        # Save BED
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# pg_cluster_merge
#
# Merge variants from each chromosome.
rule pg_cluster_merge:
    input:
        bed=expand('temp/{{asm_name}}/pg/clustered/by_chrom/{{vartype}}_{{svtype}}_{{hap}}/{chrom}.bed.gz', chrom=analib.ref.get_df_fai(REF_FAI).index)
    output:
        bed=temp('results/{asm_name}/bed/pre_merge/{vartype}_{svtype}_{hap}.bed.gz')
    wildcard_constraints:
        vartype='sv|indel|snv',
        svtype='ins|del|snv',
        hap='h1|h2'
    run:

        # Get file list, exclude 0-byte files
        file_list = [file_name for file_name in input.bed if os.stat(file_name).st_size > 0]

        pd.concat(
            [
                pd.read_csv(file_name, sep='\t') for file_name in file_list
            ],
            axis=0
        ).sort_values(
            ['#CHROM', 'POS']
        ).to_csv(
            output.bed, sep='\t', index=False, compression='gzip'
        )

# pg_cluster
#
# Cluster variants. Overlapping alignments may generate the same call from more than one contig. This step merges them
# to one call, but all contigs that support it are annotated.
rule pg_cluster:
    input:
        bed='temp/{asm_name}/pg/raw/{vartype}_{hap}.bed',
        bed_tile='results/{asm_name}/align/central_tiling_path_{hap}.bed'
    output:
        bed=temp('temp/{asm_name}/pg/clustered/by_chrom/{vartype}_{svtype}_{hap}/{chrom}.bed.gz'),
        bed_dropped='results/{asm_name}/bed/dropped/aux_{vartype}_{svtype}_{hap}_{chrom}.bed.gz'
    wildcard_constraints:
        vartype='sv|indel|snv',
        svtype='ins|del|snv',
        hap='h1|h2'
    params:
        sge_opts='-l mfree=8G -l h_rt=72:00:00'
    run:

        # Read variants
        print('Reading...')
        sys.stdout.flush()

        df = pd.read_csv(input.bed, sep='\t')

        df = df.loc[df['#CHROM'] == wildcards.chrom]

        df = df.loc[df['SVTYPE'] == wildcards.svtype.upper()].copy()

        if df.shape[0] == 0:
            print('Empty variant set')

            with open(output.bed, 'w') as out_file:
                pass

            df.to_csv(output.bed_dropped, sep='\t', index=False)

            return


        if wildcards.vartype == 'snv':

            df['REF'] = df['REF'].apply(lambda val: val.upper())
            df['ALT'] = df['ALT'].apply(lambda val: val.upper())

            df['QUERY_END'] = df['QUERY_POS'] + 1

            df['SVLEN'] = 1

        # Transform columns (by ref), set ID, sort
        df['END'] = df.apply(lambda row: (row['POS'] + 1) if row['SVTYPE'] == 'INS' else row['END'], axis=1)
        df['ID'] = analib.variant.get_variant_id(df)
        df.sort_values(['#CHROM', 'POS'], inplace=True)

        df['HAP'] = wildcards.hap

        df = analib.variant.order_variant_columns(df, tail_cols=['HAP', 'QUERY_ID', 'QUERY_POS', 'QUERY_END', 'QUERY_STRAND'])

        # Remove duplicate ID/contig pairs (can occur if the same contig is mapped multiple times to a reference location,
        # e.g. a large duplication relative the reference. Happens rarely, but will crash the pipeline.

        df.drop_duplicates(['ID', 'QUERY_ID'], inplace=True)

        df.reset_index(drop=True, inplace=True)

        # Read tiling BED and make tree
        print('Tiling...')
        sys.stdout.flush()

        df_tile = pd.read_csv(input.bed_tile, sep='\t')

        tiling_tree = collections.defaultdict(intervaltree.IntervalTree)

        for index, row in df_tile.iterrows():
            if row['END'] - row['POS'] > 0:
                tiling_tree[row['#CHROM']][row['POS']:row['END']] = row['ID']

        def get_central_id(row):
            central_set = tiling_tree[row['#CHROM']][row['POS']:row['END']]

            if len(central_set) < 1:
                return False

            return row['QUERY_ID'] == list(central_set)[0].data

        df['IS_CENTRAL'] = df.apply(get_central_id, axis=1)

        # Setup clustering data structures
        print('Setup clustering...')

        cluster_key = set()  # Primary variants (from centrally located variants)
        dropped_key = set()  # Aux variants that do not intersect a central variant
        cluster_support = collections.defaultdict(set)

        n_processed = 0

        cluster_tree = collections.defaultdict(intervaltree.IntervalTree)

        # Separate variants into central and auxiliary
        # * Central: Located in the central tiling path: when contigs overlap, choose a variant
        #   from the contig where it is further from an alignment end.
        # * Auxiliary: Variants in a region with multiple contig alignments, but closer to the end of
        #   a contig than another variant. If these intersect a central variant, annotate the central variant.
        #   For auxiliary variants that do not intersect a central variant, discard.
        df_central = df.loc[df['IS_CENTRAL']]
        df_aux = df.loc[~ df['IS_CENTRAL']]

        df_central.drop_duplicates(['ID'], inplace=True)  # Contigs with multiple alignments generate the same variant more than once

        # Check for no central variants
        if df_central.shape[0] == 0:
            print('Empty variant set')

            with open(output.bed, 'w') as out_file:
                pass

            df.to_csv(output.bed_dropped, sep='\t', index=False)

            return

        # Add primary calls from most centrally-located contig alignments
        for index, row in df_central.iterrows():
            cluster_tree[row['#CHROM']][row['POS']:row['END']] = index
            cluster_key.add(index)

        # Cluster aux calls
        print('Clustering AUX...')

        for index, row in df_aux.iterrows():
            region_start = row['POS'] - OFFSET_MAX
            region_end = row['END'] + OFFSET_MAX

            # Report progress
            if n_processed % 1000 == 0:
                print('\t\t* {} of {}'.format(n_processed + 1, df_aux.shape[0]))

            n_processed += 1

            # Get cluster
            cluster_set = cluster_tree[row['#CHROM']][region_start:region_end]

            # Check for no cluster intersection
            if not cluster_set:
                dropped_key.add(index)
                continue

            cluster_id_list = [interval.data for interval in cluster_set]

            # Check for exact ID match
            for cluster_id in cluster_id_list:
                if row['ID'] == df.loc[cluster_id, 'ID']:
                    cluster_support[cluster_id].add(index)
                    index = -1  # Flag to stop processing this variant
                    break

            # Drop SNVs unless they match by ID
            if wildcards.vartype == 'snv' and index != -1:
                dropped_key.add(index)
                continue

            # Matched variant by ID, stop processing
            if index < 0:
                continue

            # Intersect this variant with the cluster
            df_source = df.loc[[index]]
            df_target = df.loc[cluster_id_list]

            df_target.drop_duplicates('ID', inplace=True)

            df_intersect = analib.variant.nearest_by_svlen_overlap(
                df_source, df_target,
                szro_min=RO_MIN,
                offset_max=OFFSET_MAX
            )

            if not df_intersect.shape[0] > 0:
                # No match, drop aux variant
                dropped_key.add(index)

            else:
                # Add support for existing cluster.
                cluster_support[
                    df_source.loc[df_source['ID'] == df_intersect.iloc[0]['ID']].index[0]
                ].add(index)

        # Process clusters into a callset
        print('Post-cluster merging...')
        sys.stdout.flush()

        df_merge = df.loc[cluster_key]

        df_merge_support = pd.concat(
            [
                pd.Series(
                    [
                        index,
                        1 + len(cluster_support[index]),
                        ','.join(
                            [
                                df_merge.loc[index]['QUERY_ID']
                            ] + (list(
                                df.loc[cluster_support[index], 'QUERY_ID']
                            ) if cluster_support[index] else [])
                        ),
                        ','.join(
                            [
                                '{QUERY_POS}:{QUERY_END}'.format(**df_merge.loc[index])
                            ] + (list(
                                df.loc[cluster_support[index]].apply(lambda row: '{QUERY_POS}:{QUERY_END}'.format(**row), axis=1)
                            ) if cluster_support[index] else [])
                        ),
                        ','.join(
                            [
                                df_merge.loc[index]['QUERY_STRAND']
                            ] + (list(
                                df.loc[cluster_support[index], 'QUERY_STRAND']
                            ) if cluster_support[index] else [])
                        )

                    ],
                    index=['INDEX', 'TIG_N', 'TIG_SUPPORT', 'TIG_COORD', 'TIG_STRAND']
                ) for index in df_merge.index
            ],
            axis=1
        ).T

        df_merge_support.set_index('INDEX', inplace=True)

        df_merge = pd.concat([df_merge, df_merge_support], axis=1).sort_values(['#CHROM', 'POS'])

        # Move SEQ to end
        if 'SEQ' in df_merge.columns:
            df_merge = analib.variant.order_variant_columns(df_merge, tail_cols=['SEQ', ])

        df_dropped = df.loc[sorted(dropped_key)]

        # Save BED
        print('Writing...')
        sys.stdout.flush()

        df_merge.to_csv(output.bed, sep='\t', index=False, compression='gzip')

        # Save dropped bed
        df_dropped.to_csv(output.bed_dropped, sep='\t', index=False, compression='gzip')


### PrintGaps ###

# pg_sv
#
# Get SVs for one alignment.
rule pg_sv:
    input:
        aln='results/{asm_name}/align/aligned_tig_{hap}.cram'
    output:
        bed=temp('temp/{asm_name}/pg/raw/sv_{hap}.bed')
    shell:
        """samtools view -h {input.aln} | """
        """python3 {PRINT_GAPS} """
            """{REF_FA} /dev/stdin """
            """--condense 20 """
            """--minq 0 """
            """--outFile {output.bed}"""

# pg_indel
#
# Get SVs for one sample.
rule pg_indel:
    input:
        aln='results/{asm_name}/align/aligned_tig_{hap}.cram'
    output:
        bed=temp('temp/{asm_name}/pg/raw/indel_{hap}.bed')
    shell:
        """samtools view -h {input.aln} | """
        """python3 {PRINT_GAPS} """
            """{REF_FA} /dev/stdin """
            """--minLength 0 --maxLength 50 """
            """--removeAdjacentIndels """
            """--onTarget """
            """--condense 0 """
            """--minq 0 """
            """--outFile {output.bed}"""

# pg_snv
#
# Get SNVs for one sample.
rule pg_snv:
    input:
        aln='results/{asm_name}/align/aligned_tig_{hap}.cram'
    output:
        bed=temp('temp/{asm_name}/pg/raw/snv_{hap}.bed')
    shell:
        """samtools view -h {input.aln} | """
        """python3 {PRINT_GAPS} """
            """{REF_FA} /dev/stdin """
            """--minLength 0 --maxLength 0 """
            """--minq 0 """
            """--condense 0 """
            """--snv {output.bed} """
            """> /dev/null"""
