"""
Process alignments and alignment tiling paths.
"""

# align_tiling_bed
#
# Make a tiling path for the most central aligned locations in each contig. This defines the regions
# for each contig where variants should be called.
rule align_tiling_bed:
    input:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/align/central_tiling_path_{hap}.bed'
    wildcard_constraints:
        hap='h(1|2)'
    run:

        # Read contig alignments and prioritize
        df_tig = pd.read_csv(input.bed, sep='\t')

        df_tig.sort_values(['#CHROM', 'POS'], inplace=True)

        df_tig.reset_index(drop=True, inplace=True)

        this_chrom = None

        last_interval = None  # [0] POS, [1] END, [2] ID

        interval_list = list()

        # Traverse sorted alignment records
        for index, row in df_tig.iterrows():

            # Handle chrom change
            if row['#CHROM'] != this_chrom:
                if last_interval is not None:
                    interval_list.append(pd.Series(
                        [this_chrom, last_interval[0], last_interval[1], last_interval[2]],
                        index=['#CHROM', 'POS', 'END', 'ID']
                    ))

                this_chrom = row['#CHROM']
                last_interval = None

            # Compare to last interval
            if last_interval is not None:

                if row['END'] <= last_interval[1]:
                    continue  # Alignment is shorter than the last interval, cannot factor into tiling path

                if last_interval[1] > row['POS']:
                    # Contig alignments overlap

                    # Find endpoint
                    endpoint = (last_interval[1] + row['POS']) // 2

                    interval_list.append(pd.Series(
                        [this_chrom, last_interval[0], endpoint, last_interval[2]],
                        index=['#CHROM', 'POS', 'END', 'ID']
                    ))

                    # Set last interval
                    last_interval = (endpoint, row['END'], row['ID'])

                else:
                    # Contig alignments do not overlap
                    interval_list.append(pd.Series(
                        [this_chrom, last_interval[0], last_interval[1], last_interval[2]],
                        index=['#CHROM', 'POS', 'END', 'ID']
                    ))

                    # Set last interval
                    last_interval = (row['POS'], row['END'], row['ID'])

            else:
                # Set last interval
                last_interval = (row['POS'], row['END'], row['ID'])

        # Add final interval
        if last_interval is not None:
            interval_list.append(pd.Series(
                [this_chrom, last_interval[0], last_interval[1], last_interval[2]],
                index=['#CHROM', 'POS', 'END', 'ID']
            ))

        # Merge and write BED
        df_tile = pd.concat(interval_list, axis=1).T

        df_tile.to_csv(output.bed, sep='\t', index=False)

# align_single_hap_win_merge
#
# Intersect haplotype regions.
rule align_single_hap_win_merge:
    input:
        bed1='results/{asm_name}/align/depth_1/regions_h1.bed',
        bed2='results/{asm_name}/align/depth_1/regions_h2.bed'
    output:
        bed='results/{asm_name}/align/depth_1/regions_h12.bed'
    shell:
        """bedtools intersect -header -a {input.bed1} -b {input.bed2} """
        """> {output.bed}"""

# align_single_hap_win
#
# Get locations represented by a single contig.
rule align_single_hap_win:
    input:
        bed='results/{asm_name}/align/genomecov_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/align/depth_1/regions_{hap}.bed'
    wildcard_constraints:
        hap='h(1|2)'
    shell:
        """zcat {input.bed} | """
        """awk -vOFS="\\t" '($4 == 1) {{print $1, $2, $3}}' | """
        """bedtools merge -d {CONTIG_ALIGN_MERGE_DIST} """
        """> {output.bed}"""

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
        df2 = pd.read_csv(input.bed2, sep='\t', low_memory=False)

        # Merge & sort
        df = pd.concat([df1, df2], axis=0)

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Write
        df.to_csv(output.bed, sep='\t', index=False)

# align_get_read_bed
#
# Get alignment BED for one part (one aligned cell or split BAM) in one assembly.
rule align_get_read_bed:
    input:
        aln='results/{asm_name}/align/aligned_tig_{hap}.cram'
    output:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz'
    wildcard_constraints:
        hap='h(0|1|2)'
    run:

        # Get records
        clip_l = 0
        clip_r = 0

        record_list = list()

        with pysam.AlignmentFile(input.aln, 'rb') as in_file:
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

                # Save record
                record_list.append(
                pd.Series(
                    [
                        record.reference_name,
                        record.reference_start,
                        record.reference_end,
                        record.query_name,

                        tags['RG'] if 'RG' in tags else 'NA',

                        record.mapping_quality,
                        clip_l,
                        clip_r,

                        str(record.is_reverse),
                        '0x{:04x}'.format(record.flag),

                        wildcards.hap
                    ],
                    index=[
                        '#CHROM', 'POS', 'END', 'ID',
                        'RG',
                        'MAPQ', 'CLIP_L', 'CLIP_R',
                        'REV', 'FLAGS', 'HAP'
                    ]
                ))

        # Merge records
        df = pd.concat(record_list, axis=1).T

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# align_map
#
# Map contigs
rule align_map:
    input:
        asm=align_input_fasta
    output:
        aln='results/{asm_name}/align/aligned_tig_{hap}.cram',
        alni='results/{asm_name}/align/aligned_tig_{hap}.cram.crai'
    run:

        # Make temp
        temp_dir = tempfile.mkdtemp(prefix='pg_align_map_')

        try:
            shell(
                """minimap2 """
                    """--secondary=no -a -t 20 --eqx -Y """
                    """-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 """
                    """-O 5,56 -E 4,1 -B 5 """
                    """{REF_FA} {input.asm} | """
                """samtools sort -@ 4 -T {temp_dir}/sort_ | """
                """samtools view -T {REF_FA} -O CRAM -o {output.aln}; """
                """samtools index {output.aln}; """
                """touch -r {output.aln} {output.alni}"""
            )

        finally:

            # Remove temp
            shutil.rmtree(temp_dir)
