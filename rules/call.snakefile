"""
Call variants from aligned contigs.
"""

#
# INS/DEL/SNV
#


#
# Finalize variant calls
#

# pg_variant_bed
#
# Make variant BED and FASTA. Write all variant calls regardless of consensus loci.
rule call_merge_haplotypes:
    input:
       bed1='temp/{asm_name}/bed/pre_merge/{vartype}_{svtype}_h1.bed.gz',
       bed2='temp/{asm_name}/bed/pre_merge/{vartype}_{svtype}_h2.bed.gz',
       con_bed='results/{asm_name}/align/depth_1/regions_h12.bed'
    output:
        bed='results/{asm_name}/bed/{vartype}_{svtype}_h12.bed.gz',
        fa='results/{asm_name}/fa/{vartype}_{svtype}_h12.fa.gz'
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

        del(df['HAP_SRC'])
        del(df['HAP_SRC_ID'])
        del(df['HAP_AC'])
        del(df['HAP_AF'])

        df.columns = ['HAP_SRC' if val == 'HAP_SAMPLES' else val for val in df.columns]

        # Load consensus regions
        consensus_tree = collections.defaultdict(intervaltree.IntervalTree)

        df_con = pd.read_csv(input.con_bed, sep='\t', header=None, names=('#CHROM', 'POS', 'END'))

        for index, row in df_con.iterrows():
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


#
# Call alignment-truncating events
#

rule call_align_trunc:
    input:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz'
    output:
        bed_ins=temp('temp/{asm_name}/bed/align_trunc/sv_ins_{hap}.bed.gz'),
        bed_del=temp('temp/{asm_name}/bed/align_trunc/sv_del_{hap}.bed.gz'),
        bed_inv=temp('temp/{asm_name}/bed/align_trunc/sv_inv_{hap}.bed.gz')
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t')





#
# Call inversions and filter inter-INV variants
#

# call_correct_inter_inv
#
# Filter variants from inside inversions
rule call_correct_inter_inv:
    input:
        bed='temp/{asm_name}/bed/pre_merge/pre_inv_correction/{vartype}_{svtype}_{hap}.bed.gz',
        bed_inv='temp/{asm_name}/bed/pre_merge/sv_inv_{hap}.bed.gz'
    output:
        bed=temp('temp/{asm_name}/bed/pre_merge/{vartype}_{svtype}_{hap}.bed.gz'),
        bed_dropped='results/{asm_name}/bed/dropped/interinv_{vartype}_{svtype}_{hap}.bed.gz'
    wildcard_constraints:
        vartype='sv|indel|snv',
        svtype='ins|del|snv'
    run:

        # Read
        df = pd.read_csv(input.bed, sep='\t')
        df_inv = pd.read_csv(input.bed_inv, sep='\t')

        # Build tree
        invtree = collections.defaultdict(intervaltree.IntervalTree)

        for index, row in pd.read_csv(input.bed_inv, sep='\t').iterrows():
            invtree[row['#CHROM']][row['POS']:row['END']] = (row['ID'])

        # Filter
        filter_list = df.apply(
            lambda row: len(invtree[row['#CHROM']][row['POS']:row['END']]) > 0, axis=1
        )

        # Write dropped
        df.loc[filter_list].to_csv(output.bed_dropped, sep='\t', index=False, compression='gzip')

        # Write
        df.loc[~ filter_list].to_csv(output.bed, sep='\t', index=False, compression='gzip')

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
        min_svlen=config.get('inv_min_svlen', 500)
    run:

        # Read inversions
        df = pd.read_csv(input.bed, sep='\t')

        # Filter
        df_drop = df.loc[df['SVLEN'] < params.min_svlen]
        df_drop.to_csv(output.bed_dropped, sep='\t', index=False, compression='gzip')

        df = df.loc[[index for index in df.index if index not in df_drop.index]]

        df.drop_duplicates('ID', inplace=True)

        # Set common columns
        tig_region_match = re.match('^([^:]+):(\d+)-(\d+)$', df.iloc[0]['RGN_TIG_OUTER'])

        df['HAP'] = wildcards.hap
        df['TIG_N'] = 1
        df['QUERY_ID'] = tig_region_match[1]
        df['QUERY_POS'] = tig_region_match[2]
        df['QUERY_END'] = tig_region_match[3]

        # Add SEQ column
        df['SEQ'] = df.apply(
            lambda row: asmlib.seq.fa_region(
                asmlib.seq.Region(row['#CHROM'], row['POS'], row['END']),
                REF_FA
            ).upper(),
            axis=1
        )

        # Write BED
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

#
# Merge SV (ins, del), indel, SNV
#

# call_variant_inter_align
#
# Filter variants by trimmed alignments. PrintGaps was run on the full alignment, including redundantly-mapped
# contigs and reference loci.
rule call_variant_inter_align:
    input:
        bed='temp/{asm_name}/pg/raw/{vartype}_{hap}.bed',
        bed_align='results/{asm_name}/align/aligned_tig_{hap}.bed.gz'
    output:
         bed=temp('temp/{asm_name}/bed/pre_merge/pre_inv_correction/{vartype}_{svtype}_{hap}.bed.gz'),
         bed_dropped='results/{asm_name}/bed/dropped/align-trim_{vartype}_{svtype}_{hap}.bed.gz'
    wildcard_constraints:
        vartype='sv|indel|snv',
        svtype='ins|del|snv',
        hap='h1|h2'
    run:

        # Read variants
        df = pd.read_csv(input.bed, sep='\t')

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

        # Add source
        df['CALL_SOURCE'] = 'CIGAR'

        # Order columns
        tail_cols = ['HAP', 'QUERY_ID', 'QUERY_POS', 'QUERY_END', 'QUERY_STRAND', 'CALL_SOURCE']

        if 'SEQ' in df.columns:
            tail_cols += ['SEQ',]

        df = analib.variant.order_variant_columns(df, tail_cols=tail_cols)

        # Filter to accepted alignments (after trimming multiply mapped tig and reference regions)
        df_align = pd.read_csv(input.bed_align, sep='\t')

        tiling_tree = collections.defaultdict(intervaltree.IntervalTree)

        for index, row in df_align.iterrows():
            tiling_tree[row['#CHROM']][row['POS']:row['END']] = (row['QUERY_ID'], row['QUERY_POS'], row['QUERY_END'])

        def is_pass_region(row):
            for pass_record in tiling_tree[row['#CHROM']][row['POS']:row['END']]:
                pass_record_data = pass_record.data

                if (
                    row['QUERY_ID'] == pass_record_data[0] and
                    row['QUERY_POS'] >= pass_record_data[1] and
                    row['QUERY_END'] <= pass_record_data[2]
                ):
                    return True

            return False


        df['PASS_REGION'] = df.apply(is_pass_region, axis=1)

        # Split variants
        df_dropped = df.loc[~ df['PASS_REGION']].copy()
        df = df.loc[df['PASS_REGION']].copy()

        del(df['PASS_REGION'])
        del(df_dropped['PASS_REGION'])

        # Check for duplicate variant IDs (should not occur)
        # dup_id = {key for key,val in collections.Counter(df['ID']).items() if val > 1}
        #
        # if dup_id:
        #     raise RuntimeError('Found {} duplicate variant IDs after filtering: {}{}'.format(
        #         len(dup_id),
        #         ', '.join(sorted(dup_id)[:3]),
        #         '...' if len(dup_id) > 3 else ''
        #     ))
        df.drop_duplicates('ID', inplace=True)  # DBGTMP: Fix process so duplicates are never generated (filtering is limited)

        # Save BED
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

        # Save dropped bed
        df_dropped.to_csv(output.bed_dropped, sep='\t', index=False, compression='gzip')


#
# PrintGaps
#

# call_printgaps_sv
#
# Get SVs for one alignment.
rule call_printgaps_sv:
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

# call_printgaps_indel
#
# Get SVs for one sample.
rule call_printgaps_indel:
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

# call_printgaps_snv
#
# Get SNVs for one sample.
rule call_printgaps_snv:
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
