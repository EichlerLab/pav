"""
Process alignments and alignment tiling paths.
"""

# align_genomecov
#
# Get genome coverage.
rule align_genomecov:
    input:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/align/genomecov_{hap}.bed.gz'
    shell:
        """{{ \n"""
        """    echo -e "#CHROM\tPOS\tEND\tDEPTH"; \n"""
        """    bedtools genomecov -bga -i {input.bed} -g {REF_FAI}; \n"""
        """}} | """
        """gzip > {output.bed}"""

# align_merge_h12_read_bed
#
# Alignment table for all reads.
rule align_merge_h12_read_bed:
    input:
        bed1='results/{asm_name}/align/aligned_tig_h1.bed.gz',
        bed2='results/{asm_name}/align/aligned_tig_h2.bed.gz'
    output:
        bed='results/{asm_name}/align/aligned_tig_h12.bed.gz'
    run:

        # Read
        df1 = pd.read_csv(input.bed1, sep='\t', low_memory=False)
        df1['HAP'] = 'h1'

        df2 = pd.read_csv(input.bed2, sep='\t', low_memory=False)
        df1['HAP'] = 'h2'

        # Merge & sort
        df = pd.concat([df1, df2], axis=0)

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write
        df.to_csv(output.bed, sep='\t', index=False)


# align_sort_cram
#
# Sort alignments.
rule align_sort_cram:
    input:
        sam='temp/{asm_name}/align/aligned_tig_{hap}.sam'
    output:
        aln='results/{asm_name}/align/aligned_tig_{hap}.cram',
        alni='results/{asm_name}/align/aligned_tig_{hap}.cram.crai'
    run:

        # Make temp
        temp_dir = tempfile.mkdtemp(prefix='pg_align_map_')

        try:
            shell(
                """samtools sort -@ 4 -T {temp_dir}/sort_ {input.sam} | """
                """samtools view -T {REF_FA} -O CRAM -o {output.aln}; """
                """samtools index {output.aln}; """
                """touch -r {output.aln} {output.alni}"""
            )

        finally:

            # Remove temp
            shutil.rmtree(temp_dir)


# align_cut_tig_overlap
#
# Cut contig alignment overlaps
rule align_cut_tig_overlap:
    input:
        bed='results/{asm_name}/align/bed/aligned_tig_uncut_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz'
    params:
        chrom_cluster=asmlib.util.as_bool(config.get('chrom_cluster', True))  # Assembly was clustered by chromosome and first part of chromosome name before "_" is the cluster name.
    run:

        # Read uncut alignments
        df = pd.read_csv(input.bed, sep='\t')

        # Add fields for the number of bases that are removed from each end
        df['CUT_REF_L'] = 0
        df['CUT_REF_R'] = 0
        df['CUT_TIG_L'] = 0
        df['CUT_TIG_R'] = 0

        # Sort by contig alignment length
        df['QUERY_LEN'] = df['QUERY_END'] - df['QUERY_POS']
        df['SUB_LEN'] = df['END'] - df['POS']

        df.sort_values(['QUERY_ID', 'QUERY_LEN'], ascending=(True, False), inplace=True)

        df.reset_index(inplace=True, drop=True)

        # Find max cluster match for each chromosome
        if params.chrom_cluster:
            df['CLUSTER'] = df['QUERY_ID'].apply(lambda val: val.split('_')[0])
            max_cluster = {chrom: asmlib.align.get_max_cluster(df, chrom) for chrom in set(df['#CHROM'])}

            df['CLUSTER_MATCH'] = df.apply(lambda row: row['CLUSTER'] == max_cluster[row['#CHROM']], axis=1)
            df['CLUSTER_MATCH'] = df.apply(lambda row: row['CLUSTER_MATCH'] if max_cluster[row['#CHROM']] is not None else np.nan, axis=1)

        else:
            df['CLUSTER_MATCH'] = np.nan

        # Resolve overlapping contig regions aligned (one contig region aligned more than once)
        iter_index_l = 0
        index_max = df.shape[0]

        while iter_index_l < index_max:
            iter_index_r = iter_index_l + 1

            while iter_index_r < index_max and df.loc[iter_index_l, 'QUERY_ID'] == df.loc[iter_index_r, 'QUERY_ID']:

                # Skip if one record was already removed
                if df.loc[iter_index_l, 'INDEX'] < 0 or df.loc[iter_index_r, 'INDEX'] < 0:
                    iter_index_r += 1
                    continue

                # Get indices ordered by contig placement
                if df.loc[iter_index_l, 'QUERY_POS'] <= df.loc[iter_index_r, 'QUERY_POS']:
                    index_l = iter_index_l
                    index_r = iter_index_r
                else:
                    index_l = iter_index_r
                    index_r = iter_index_l

                # Check for overlaps
                if df.loc[index_r, 'QUERY_POS'] < df.loc[index_l, 'QUERY_END']:
                    # Found overlapping records
                    # print('Tig Overlap: {}-{} ({}:{}-{} vs {}:{}-{}) [iter {}, {}]'.format(
                    #     df.loc[index_l, 'INDEX'], df.loc[index_r, 'INDEX'],
                    #     df.loc[index_l, 'QUERY_ID'], df.loc[index_l, 'QUERY_POS'], df.loc[index_l, 'QUERY_END'],
                    #     df.loc[index_r, 'QUERY_ID'], df.loc[index_r, 'QUERY_POS'], df.loc[index_r, 'QUERY_END'],
                    #     iter_index_l, iter_index_r
                    # ))

                    # Check for record fully contained within another
                    if df.loc[index_r, 'QUERY_END'] <= df.loc[index_l, 'QUERY_END']:
                        # print('\t* Fully contained')

                        df.loc[index_r, 'INDEX'] = -1

                    else:

                        record_l, record_r = asmlib.align.trim_alignments(df.loc[index_l], df.loc[index_r], 'query')

                        if record_l is not None and record_r is not None:
                            df.loc[index_l] = record_l
                            df.loc[index_r] = record_r

                        # print('\t* Trimmed')

                # Next r record
                iter_index_r += 1

            # Next l record
            iter_index_l += 1


        # # Remove discarded records and re-sort
        #
        # df = df.loc[df['INDEX'] >= 0]
        #
        # df['QUERY_LEN'] = df['QUERY_END'] - df['QUERY_POS']
        #
        # df.sort_values(['#CHROM', 'QUERY_LEN'], ascending=(True, False), inplace=True)
        #
        # df.reset_index(inplace=True, drop=True)
        #
        # # Resolve overlapping same-contig alignments (more than one contig region mapped to the same reference region)
        # iter_index_l = 0
        # index_max = df.shape[0]
        #
        # while iter_index_l < index_max:
        #     iter_index_r = iter_index_l + 1
        #
        #     while (
        #             iter_index_r < index_max and
        #             df.loc[iter_index_l, 'QUERY_ID'] == df.loc[iter_index_r, 'QUERY_ID'] and
        #             df.loc[iter_index_l, '#CHROM'] == df.loc[iter_index_r, '#CHROM']
        #     ):
        #
        #         # Skip if one record was already removed
        #         if df.loc[iter_index_l, 'INDEX'] < 0 or df.loc[iter_index_r, 'INDEX'] < 0:
        #             iter_index_r += 1
        #             continue
        #
        #         # Get indices ordered by contig placement
        #         if df.loc[iter_index_l, 'POS'] <= df.loc[iter_index_r, 'POS']:
        #             index_l = iter_index_l
        #             index_r = iter_index_r
        #         else:
        #             index_l = iter_index_r
        #             index_r = iter_index_l
        #
        #         # Check for overlaps
        #         if df.loc[index_r, 'POS'] < df.loc[index_l, 'END']:
        #             # Found overlapping records
        #             # print('Ref-Tig-Same Overlap: {}-{} ({}:{}-{} vs {}:{}-{}) ({}:{}-{}, {}-{}) [iter {}, {}]'.format(
        #             #     df.loc[index_l, 'INDEX'], df.loc[index_r, 'INDEX'],
        #             #     df.loc[index_l, 'QUERY_ID'], df.loc[index_l, 'QUERY_POS'], df.loc[index_l, 'QUERY_END'],
        #             #     df.loc[index_r, 'QUERY_ID'], df.loc[index_r, 'QUERY_POS'], df.loc[index_r, 'QUERY_END'],
        #             #     df.loc[index_l, '#CHROM'], df.loc[index_l, 'POS'], df.loc[index_l, 'END'], df.loc[index_r, 'POS'], df.loc[index_r, 'END'],
        #             #     iter_index_l, iter_index_r
        #             # ))
        #
        #             # Check for record fully contained within another
        #             if df.loc[index_r, 'END'] <= df.loc[index_l, 'END']:
        #                 # print('\t* Fully contained')
        #
        #                 df.loc[index_r, 'INDEX'] = -1
        #
        #             else:
        #
        #                 record_l, record_r = asmlib.align.trim_alignments(df.loc[index_l], df.loc[index_r], 'subject')
        #
        #                 if record_l is not None and record_r is not None:
        #                     df.loc[index_l] = record_l
        #                     df.loc[index_r] = record_r
        #
        #                 # print('\t* Trimmed')
        #
        #         # Next r record
        #         iter_index_r += 1
        #
        #     # Next l record
        #     iter_index_l += 1


        # Remove discarded records and re-sort

        df = df.loc[df['INDEX'] >= 0]

        df['QUERY_LEN'] = df['QUERY_END'] - df['QUERY_POS']

        df.sort_values(['#CHROM', 'QUERY_LEN'], ascending=(True, False), inplace=True)

        df.reset_index(inplace=True, drop=True)

        # Resolve overlapping contig alignments relative to the reference
        iter_index_l = 0
        index_max = df.shape[0]

        while iter_index_l < index_max:
            iter_index_r = iter_index_l + 1

            while (
                    iter_index_r < index_max and
                    df.loc[iter_index_l, '#CHROM'] == df.loc[iter_index_r, '#CHROM']
            ):

                # Skip if one record was already removed
                if df.loc[iter_index_l, 'INDEX'] < 0 or df.loc[iter_index_r, 'INDEX'] < 0:
                    iter_index_r += 1
                    continue

                # Get indices ordered by contig placement
                if df.loc[iter_index_l, 'POS'] <= df.loc[iter_index_r, 'POS']:
                    index_l = iter_index_l
                    index_r = iter_index_r
                else:
                    index_l = iter_index_r
                    index_r = iter_index_l

                # Check for overlaps
                if df.loc[index_r, 'POS'] < df.loc[index_l, 'END']:
                    # Found overlapping records
                    # print('Ref-Tig-Other Overlap: {}-{} ({}:{}-{} vs {}:{}-{}) ({}:{}-{}, {}-{}) [iter {}, {}]'.format(
                    #     df.loc[index_l, 'INDEX'], df.loc[index_r, 'INDEX'],
                    #     df.loc[index_l, 'QUERY_ID'], df.loc[index_l, 'QUERY_POS'], df.loc[index_l, 'QUERY_END'],
                    #     df.loc[index_r, 'QUERY_ID'], df.loc[index_r, 'QUERY_POS'], df.loc[index_r, 'QUERY_END'],
                    #     df.loc[index_l, '#CHROM'], df.loc[index_l, 'POS'], df.loc[index_l, 'END'], df.loc[index_r, 'POS'], df.loc[index_r, 'END'],
                    #     iter_index_l, iter_index_r
                    # ))

                    # Check for record fully contained within another
                    if df.loc[index_r, 'END'] <= df.loc[index_l, 'END']:
                        # print('\t* Fully contained')

                        df.loc[index_r, 'INDEX'] = -1

                    else:

                        record_l, record_r = asmlib.align.trim_alignments(df.loc[index_l], df.loc[index_r], 'subject')

                        if record_l is not None and record_r is not None:
                            df.loc[index_l] = record_l
                            df.loc[index_r] = record_r

                        # print('\t* Trimmed')

                # Next r record
                iter_index_r += 1

            # Next l record
            iter_index_l += 1


        # Clean and re-sort
        df = df.loc[df['INDEX'] >= 0]

        df = df.loc[(df['END'] - df['POS']) > 0]  # Should never occur, but don't allow 0-length records
        df = df.loc[(df['QUERY_END'] - df['QUERY_POS']) > 0]

        df.sort_values('INDEX', inplace=True)

        del(df['QUERY_LEN'])
        del(df['SUB_LEN'])

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')



# align_get_read_bed
#
# Get alignment BED for one part (one aligned cell or split BAM) in one assembly.
rule align_get_read_bed:
    input:
        sam='temp/{asm_name}/align/aligned_tig_{hap}.sam'
    output:
        bed='results/{asm_name}/align/bed/aligned_tig_uncut_{hap}.bed.gz'
    wildcard_constraints:
        hap='h(0|1|2)'
    run:

        # Get records
        clip_l = 0
        clip_r = 0

        record_list = list()

        with pysam.AlignmentFile(input.sam, 'rb') as in_file:
            for record in in_file:

                # Skipped unmapped reads
                if record.is_unmapped:
                    continue

                # Read tags
                tags = dict(record.get_tags())

                # Get clipping
                cigar_tuples = record.cigartuples

                l_index = 0 if cigar_tuples[0][0] != 5 else 1
                r_index = -1 if cigar_tuples[-1][0] != 5 else -2

                clip_l = cigar_tuples[l_index][1] if cigar_tuples[l_index][0] == 4 else 0
                clip_r = cigar_tuples[r_index][1] if cigar_tuples[r_index][0] == 4 else 0

                # Disallow alignment match (M) in CIGAR (requires =X for base match/mismatch)
                if 'M' in record.cigarstring:
                    raise RuntimeError((
                        'Found alignment match CIGAR operation (M) for record {} (Start = {}:{}): '
                        'Alignment requires CIGAR base-level match/mismatch (=X)'
                    ).format(record.query_name, record.reference_name, record.reference_start))

                # Save record
                record_list.append(pd.Series(
                    [
                        record.reference_name,
                        record.reference_start,
                        record.reference_end,

                        record.query_name,
                        record.query_alignment_start,
                        record.query_alignment_end,

                        tags['RG'] if 'RG' in tags else 'NA',

                        record.mapping_quality,
                        clip_l,
                        clip_r,

                        str(record.is_reverse),
                        '0x{:04x}'.format(record.flag),

                        wildcards.hap,
                        record.cigarstring
                    ],
                    index=[
                        '#CHROM', 'POS', 'END',
                        'QUERY_ID', 'QUERY_POS', 'QUERY_END',
                        'RG',
                        'MAPQ', 'CLIP_L', 'CLIP_R',
                        'REV', 'FLAGS', 'HAP',
                        'CIGAR'
                    ]
                ))

        # Merge records
        df = pd.concat(record_list, axis=1).T

        df.sort_values(['#CHROM', 'POS', 'END', 'QUERY_ID'], ascending=[True, True, False, True], inplace=True)

        df.insert(3, 'INDEX', list(range(df.shape[0])))

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# align_map
#
# Map contigs as SAM. Pull read information from the SAM before sorting and writing CRAM since tool tend to change
# "=X" to "M" in the CIGAR.
rule align_map:
    input:
        fa='align/{asm_name}/contigs_{hap}.fa.gz'
    output:
        sam=temp('temp/{asm_name}/align/aligned_tig_{hap}.sam')
    shell:
        """minimap2 """
            """-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 """
            """--secondary=no -a -t 20 --eqx -Y """
            """-O 5,56 -E 4,1 -B 5 """
            """{REF_FA} {input.fa} """
            """> {output.sam}"""

# align_index
#
# Get FASTA files.
rule align_index:
    input:
        fa=align_input_fasta
    output:
        fa='align/{asm_name}/contigs_{hap}.fa.gz'
    run:

        # Determine if file is BGZF compressed
        is_bgzf = False

        try:
            with Bio.bgzf.open(input.fa, 'r') as in_file_test:
                is_bgzf = True

        except ValueError:
            pass

        # Copy or compress
        if is_bgzf:

            # Copy file if already compressed
            shutil.copyfile(input.fa, output.fa)

        else:
            # Compress to BGZF

            is_gz = False

            try:
                with gzip.open(input.fa, 'r') as in_file_test:

                    line = next(in_file_test)

                    is_gz = True

            except OSError:
                pass

            if is_gz:
                # Re-compress to BGZF

                with gzip.open(input.fa, 'rb') as in_file:
                    with Bio.bgzf.open(output.fa, 'wb') as out_file:
                        for line in in_file:
                            out_file.write(line)

            else:
                # Compress plain text

                with open(input.fa, 'r') as in_file:
                    with Bio.bgzf.open(output.fa, 'wb') as out_file:
                        for line in in_file:
                            out_file.write(line)

        # Index
        shell("""samtools faidx {output.fa}""")
