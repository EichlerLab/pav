"""
Process alignments and alignment tiling paths.
"""


import gzip
import os
import pandas as pd

import pavlib
import svpoplib


global config
global expand
global shell
global temp
global get_config
global ASM_TABLE
global REF_FA
global REF_FAI


#
# Rules
#

# Run all alignments
localrules: align_all

rule align_all:
    input:
        bed_depth=lambda wildcards: pavlib.pipeline.expand_pattern(
            'results/{asm_name}/align/trim-{trim}/depth_tig_{hap}.bed.gz', ASM_TABLE,
            trim=('none', 'tig', 'tigref')
        )

# Create a depth BED file for alignments.
rule align_depth_bed:
    input:
        bed='results/{asm_name}/align/trim-{trim}/aligned_tig_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/align/trim-{trim}/depth_tig_{hap}.bed.gz'
    run:

        pavlib.align.align_bed_to_depth_bed(
            pd.read_csv(input.bed, sep='\t', dtype={'#CHROM': str}),
            svpoplib.ref.get_df_fai(REF_FAI)
        ).to_csv(
            output.bed, sep='\t', index=False, compression='gzip'
        )

# Cut contig alignment overlaps in reference coordinates
rule align_trim_tigref:
    input:
        bed='results/{asm_name}/align/trim-tig/aligned_tig_{hap}.bed.gz',
        tig_fai='temp/{asm_name}/align/contigs_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/trim-tigref/aligned_tig_{hap}.bed.gz'
    params:
        min_trim_tig_len=lambda wildcards: int(get_config(wildcards, 'min_trim_tig_len', 1000)),  # Minimum aligned tig length
        redundant_callset=lambda wildcards: pavlib.util.as_bool(get_config(wildcards, 'redundant_callset', False))
    run:

        # Trim alignments
        df = pavlib.align.trim_alignments(
            pd.read_csv(input.bed, sep='\t', dtype={'#CHROM': str}),  # Untrimmed alignment BED
            params.min_trim_tig_len,  # Minimum contig length
            input.tig_fai,  # Path to alignment FASTA FAI
            match_tig=params.redundant_callset,  # Redundant callset, trim reference space only for records with matching IDs
            mode='ref'
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# Cut contig alignment overlaps in contig coordinates
rule align_trim_tig:
    input:
        bed='results/{asm_name}/align/trim-none/aligned_tig_{hap}.bed.gz',
        tig_fai='temp/{asm_name}/align/contigs_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/trim-tig/aligned_tig_{hap}.bed.gz'
    params:
        min_trim_tig_len=lambda wildcards: int(get_config(wildcards, 'min_trim_tig_len', 1000))  # Minimum aligned tig length
    run:

        # Trim alignments
        df = pavlib.align.trim_alignments(
            pd.read_csv(input.bed, sep='\t', dtype={'#CHROM': str}),  # Untrimmed alignment BED
            params.min_trim_tig_len,  # Minimum contig length
            input.tig_fai,  # Path to alignment FASTA FAI
            mode='tig'
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# Get alignment BED for one part (one aligned cell or split BAM) in one assembly.
rule align_get_read_bed:
    input:
        sam='temp/{asm_name}/align/trim-none/aligned_tig_{hap}.sam.gz',
        tig_fai='temp/{asm_name}/align/contigs_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/trim-none/aligned_tig_{hap}.bed.gz',
        align_head='results/{asm_name}/align/trim-none/aligned_tig_{hap}.headers.gz'
    run:

        # Write an empty file if SAM is emtpy
        if os.stat(input.sam).st_size == 0:

            pd.DataFrame(
                [],
                columns=[
                    '#CHROM', 'POS', 'END',
                    'INDEX',
                    'QRY_ID', 'QRY_POS', 'QRY_END',
                    'QRY_MAPPED',
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
        df = pavlib.align.get_align_bed(input.sam, df_tig_fai, wildcards.hap)

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

        # Add batch ID for CIGAR calling (calls in batches)
        df['CALL_BATCH'] = df['INDEX'].apply(lambda val: val % pavlib.cigarcall.CALL_CIGAR_BATCH_COUNT)

        # Add trimming fields
        df['TRIM_REF_L'] = 0
        df['TRIM_REF_R'] = 0
        df['TRIM_QRY_L'] = 0
        df['TRIM_QRY_R'] = 0

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# Map contigs as SAM. Pull read information from the SAM before sorting and writing CRAM since tool tend to change
# "=X" to "M" in the CIGAR.
rule align_map:
    input:
        ref_fa='data/ref/ref.fa.gz',
        fa=lambda wildcards: 'temp/{asm_name}/align/contigs_{hap}.fa.gz' if get_config(wildcards, 'aligner', 'minimap2') != 'lra' else 'temp/{asm_name}/align/contigs_{hap}.fa',
        fai=lambda wildcards: 'temp/{asm_name}/align/contigs_{hap}.fa.gz.fai' if get_config(wildcards, 'aligner', 'minimap2') != 'lra' else [],
        gzi=lambda wildcards: 'temp/{asm_name}/align/contigs_{hap}.fa.gz.gzi' if get_config(wildcards, 'aligner', 'minimap2') != 'lra' else [],
        gli=lambda wildcards: 'data/ref/ref.fa.gz.gli' if get_config(wildcards, 'aligner', 'minimap2') == 'lra' else [],
        mmi=lambda wildcards: 'data/ref/ref.fa.gz.mms' if get_config(wildcards, 'aligner', 'minimap2') == 'lra' else []
    output:
        sam=temp('temp/{asm_name}/align/trim-none/aligned_tig_{hap}.sam.gz')
    params:
        aligner=lambda wildcards: get_config(wildcards, 'aligner', 'minimap2'),
        minimap2_params=lambda wildcards: get_config(wildcards, 'minimap2_params', '-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5')
    threads: 24
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

# Get FASTA files.
rule align_get_tig_fa:
    input:
        fa=lambda wildcards: pavlib.pipeline.get_rule_input_list(wildcards.asm_name, wildcards.hap, ASM_TABLE, get_config(wildcards))
    output:
        fa=temp('temp/{asm_name}/align/contigs_{hap}.fa.gz'),
        fai=temp('temp/{asm_name}/align/contigs_{hap}.fa.gz.fai'),
        gzi=temp('temp/{asm_name}/align/contigs_{hap}.fa.gz.gzi')
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

        # Link or generate FASTA
        is_link = False

        if len(input_tuples) == 1 and input_tuples[0][1] == 'fasta' and input_tuples[0][0].lower().endswith('.gz'):
            os.symlink(input_tuples[0][0], output.fa)
            is_link = True

        else:
            # Merge/write FASTA (or empty files if there is no input)
            pavlib.pipeline.input_tuples_to_fasta(input_tuples, output.fa)
            is_link = False

        # Index
        if is_link and os.path.isfile(input_tuples[0][1] + '.fai') and os.path.isfile(input_tuples[0][1] + '.gzi'):
            os.symlink(os.path.abspath(input_tuples[0][1] + '.fai'), output.fai)
            os.symlink(os.path.abspath(input_tuples[0][1] + '.gzi'), output.gzi)

        else:
            shell("""samtools faidx {output.fa}""")


#
# Utilities
# Not needed for PAV, but useful rules for troubleshooting.
#

def _align_cram_all(wildcards):

    if 'trim' in config:
        trim = config['trim'].strip().split(',')
    else:
        trim = ('tigref',)

    return pavlib.pipeline.expand_pattern(
        'results/{asm_name}/align/cram/pav_tig_trim-{trim}_{hap}.cram', ASM_TABLE, trim=trim
    )

# Get CRAM files
localrules: align_cram_all

rule align_cram_all:
    input:
        cram=_align_cram_all


# Reconstruct CRAM from alignment BED files after trimming redundantly mapped bases (post-cut).
rule align_get_cram:
    input:
        bed='results/{asm_name}/align/trim-{trim}/aligned_tig_{hap}.bed.gz',
        fa='temp/{asm_name}/align/contigs_{hap}.fa.gz',
        align_head='results/{asm_name}/align/trim-none/aligned_tig_{hap}.headers.gz',
        ref_fa='data/ref/ref.fa.gz'
    output:
        cram='results/{asm_name}/align/cram/pav_tig_trim-{trim}_{hap}.cram',
        crai='results/{asm_name}/align/cram/pav_tig_trim-{trim}_{hap}.cram.crai'
    params:
        sam_tag=lambda wildcards: fr'@PG\tID:PAV-{wildcards.trim}\tPN:PAV\tVN:{pavlib.constants.get_version_string()}\tDS:PAV Alignment trimming {pavlib.align.TRIM_DESC[wildcards.trim]}'
    shell:
        """python3 {PIPELINE_DIR}/scripts/reconstruct_sam.py """
            """--bed {input.bed} --fasta {input.fa} --headers {input.align_head} --tag "{params.sam_tag}" | """
        """samtools view -T {input.ref_fa} -O CRAM -o {output.cram} && """
        """samtools index {output.cram}"""
