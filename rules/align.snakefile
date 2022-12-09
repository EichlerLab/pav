"""
Process alignments and alignment tiling paths.
"""

#
# Definitions
#

def _align_map_cpu(wildcards, config):

    if 'map_threads' in config:
        try:
            return int(config['map_threads'])
        except ValueError as ex:
            raise ValueError('Config parameter "map_threads" is not an integer: {map_threads}'.format(**config))

    return 12


#
# Alignment generation and processing
#

# align_cut_tig_overlap
#
# Cut contig alignment overlaps
rule align_cut_tig_overlap:
    input:
        bed='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.bed.gz',
        tig_fai='temp/{asm_name}/align/contigs_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz'
    params:
        min_trim_tig_len=lambda wildcards: np.int32(get_config(wildcards, 'min_trim_tig_len', 1000)),  # Minimum aligned tig length
        redundant_callset=lambda wildcards: pavlib.util.as_bool(get_config(wildcards, 'redundant_callset', False))
    run:

        # Trim alignments
        df = pavlib.align.trim_alignments(
            pd.read_csv(input.bed, sep='\t'),  # Untrimmed alignment BED
            params.min_trim_tig_len,  # Minimum contig length
            input.tig_fai,  # Path to alignment FASTA FAI
            match_tig=params.redundant_callset  # Redundant callset, trim reference space only for records with matching IDs
        )

        # Add batch ID for CIGAR calling (calls in batches)
        df['CALL_BATCH'] = df['INDEX'].apply(lambda val: val % pavlib.cigarcall.CALL_CIGAR_BATCH_COUNT)

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# align_get_read_bed
#
# Get alignment BED for one part (one aligned cell or split BAM) in one assembly.
rule align_get_read_bed:
    input:
        sam='temp/{asm_name}/align/pre-cut/aligned_tig_{hap}.sam.gz',
        tig_fai='temp/{asm_name}/align/contigs_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.bed.gz',
        align_head='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.headers.gz'
    params:
        chrom_cluster=lambda wildcards: pavlib.util.as_bool(get_config(wildcards, 'chrom_cluster', False))  # Assembly was clustered by chromosome and first part of chromosome name before "_" is the cluster name.
    wildcard_constraints:
        hap='h(0|1|2)'
    run:

        # Write an empty file if SAM is emtpy
        if os.stat(input.sam).st_size == 0:

            pd.DataFrame(
                [],
                columns=[
                    '#CHROM', 'POS', 'END',
                    'INDEX',
                    'QUERY_ID', 'QUERY_POS', 'QUERY_END',
                    'QUERY_TIG_POS', 'QUERY_TIG_END',
                    'RG', 'AO',
                    'MAPQ',
                    'REV', 'FLAGS', 'HAP',
                    'CIGAR'
                ]
            ).to_csv(
                output.bed, sep='\t', index=False, compression='gzip'
            )

            with open(output.align_head, 'w') as out_file:
                pass

            return

        # Read FAI
        df_tig_fai = svpoplib.ref.get_df_fai(input.tig_fai)
        df_tig_fai.index = df_tig_fai.index.astype(str)

        # Read alignments as a BED file.
        df = pavlib.align.get_align_bed(input.sam, df_tig_fai, wildcards.hap, params.chrom_cluster)

        # Write SAM headers
        with gzip.open(input.sam, 'rt') as in_file:
            with gzip.open(output.align_head, 'wt') as out_file:

                line = next(in_file)

                while True:

                    if not line.strip():
                        continue

                    if not line.startswith('@'):
                        break

                    out_file.write(line)

                    try:
                        line = next(in_file)
                    except StopIteration:
                        break

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# align_map
#
# Map contigs as SAM. Pull read information from the SAM before sorting and writing CRAM since tool tend to change
# "=X" to "M" in the CIGAR.
rule align_map:
    input:
        ref_fa='data/ref/ref.fa.gz',
        fa=lambda wildcards: 'temp/{asm_name}/align/contigs_{hap}.fa.gz' if get_config(wildcards, 'aligner', 'minimap2') != 'lra' else 'temp/{asm_name}/align/contigs_{hap}.fa',
        gli=lambda wildcards: 'data/ref/ref.fa.gz.gli' if get_config(wildcards, 'aligner', 'minimap2') == 'lra' else [],
        mmi=lambda wildcards: 'data/ref/ref.fa.gz.mmi' if get_config(wildcards, 'aligner', 'minimap2') == 'lra' else []
    output:
        sam=temp('temp/{asm_name}/align/pre-cut/aligned_tig_{hap}.sam.gz')
    params:
        cpu=lambda wildcards: _align_map_cpu(wildcards, get_config(wildcards)),
        aligner=lambda wildcards: get_config(wildcards, 'aligner', 'minimap2'),
        minimap2_params=lambda wildcards: get_config(wildcards, 'minimap2_params', '-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5')
    threads: lambda wildcards: _align_map_cpu(wildcards, get_config(wildcards))
    run:

        # Get aligner
        if params.aligner not in {'minimap2', 'lra'}:
            raise RuntimeError('Unknown "aligner" parameter in config: {}'.format(params.aligner))

        aligner = params.aligner

        # Write an empty file if input is empty
        if os.stat(input.fa).st_size == 0:
            with open(output.sam, 'w') as out_file:
                pass

        else:

            # Align
            if aligner == 'minimap2':
                shell(
                    """minimap2 """
                        """{params.minimap2_params} """
                        """--secondary=no -a -t {threads} --eqx -Y """
                        """{input.ref_fa} {input.fa} | """
                        """awk -vOFS="\\t" '($1 !~ /^@/) {{$10 = "*"; $11 = "*"}} {{print}}' | """
                        """gzip > {output.sam}"""
                )

            if aligner == 'lra':
                shell(
                    """lra align {input.ref_fa} {input.fa} -CONTIG -p s -t {threads} | """
                    """awk -vOFS="\\t" '($1 !~ /^@/) {{$10 = "*"; $11 = "*"}} {{print}}' | """
                    """gzip > {output.sam}"""
                )

# align_uncompress_tig
#
# Uncompress contig for aligners that cannot read gzipped FASTAs.
rule align_uncompress_tig:
    input:
        fa='temp/{asm_name}/align/contigs_{hap}.fa.gz'
    output:
        fa='temp/{asm_name}/align/contigs_{hap}.fa'
    run:

        if os.stat(input.fa).st_size > 0:
            shell(
                """zcat {input.fa} > {output.fa}"""
            )
        else:
            with open(output.fa, 'w') as out_file:
                pass

# align_get_tig_fa
#
# Get FASTA files.
rule align_get_tig_fa:
    input:
        fa=lambda wildcards: pavlib.pipeline.get_rule_input_list(wildcards.asm_name, wildcards.hap, ASM_TABLE, get_config(wildcards))
    output:
        fa=temp('temp/{asm_name}/align/contigs_{hap}.fa.gz'),
        fai=temp('temp/{asm_name}/align/contigs_{hap}.fa.gz.fai')
    run:

        # Get input files
        input_list = input.fa if 'fa' in input.keys() else []

        input_tuples, fofn_list = pavlib.pipeline.expand_input(
            pavlib.pipeline.get_asm_input_list(wildcards.asm_name, wildcards.hap, ASM_TABLE, config)
        )

        # Report input sources
        if input_tuples is not None:
            for file_name, file_format in input_tuples:
                print(f'Input: {wildcards.asm_name} {wildcards.hap}: {file_name} ({file_format})')
        else:
            print(f'No input sources: {wildcards.asm_name} {wildcards.hap}')

        # Merge/write FASTA (or empty files if there is no input)
        pavlib.pipeline.input_tuples_to_fasta(input_tuples, output.fa)

        # Index
        if os.stat(output.fa).st_size > 0:
            shell("""samtools faidx {output.fa}""")

        else:
            with open(output.fai, 'w') as out_file:
                pass


#
# Utilities
#
# Not needed for PAV, but useful rules for troubleshooting.

# align_get_cram_postcut
#
# Reconstruct CRAM from alignment BED files after trimming redundantly mapped bases (post-cut).
rule align_get_cram_postcut:
    input:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz',
        fa='temp/{asm_name}/align/contigs_{hap}.fa.gz',
        align_head='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.headers.gz',
        ref_fa='data/ref/ref.fa.gz'

    output:
        cram='results/{asm_name}/align/aligned_tig_{hap}.cram'
    shell:
        """python3 {PIPELINE_DIR}/scripts/reconstruct_sam.py """
            """--bed {input.bed} --fasta {input.fa} --headers {input.align_head} | """
        """samtools view -T {input.ref_fa} -O CRAM -o {output.cram} && """
        """samtools index {output.cram}"""

# align_get_cram_precut
#
# Reconstruct CRAM from alignment BED files before trimming redundantly mapped bases (post-cut).
rule align_get_cram_precut:
    input:
        bed='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.bed.gz',
        fa='temp/{asm_name}/align/contigs_{hap}.fa.gz',
        align_head='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.headers.gz',
        ref_fa='data/ref/ref.fa.gz'
    output:
        cram='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.cram'
    shell:
        """python3 {PIPELINE_DIR}/scripts/reconstruct_sam.py """
            """--bed {input.bed} --fasta {input.fa} --headers {input.align_head} | """
        """samtools view -T {input.ref_fa} -O CRAM -o {output.cram} && """
        """samtools index {output.cram}"""
