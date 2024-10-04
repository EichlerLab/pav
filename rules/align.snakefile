"""
Process alignments and alignment tiling paths.
"""


import Bio.bgzf
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
global get_override_config
global ASM_TABLE
global PIPELINE_DIR
global REF_FA
global REF_FAI


#
# Rules
#

# Run all alignments
localrules: align_all

rule align_all:
    input:
        bed_align=lambda wildcards: pavlib.pipeline.expand_pattern(
            'results/{asm_name}/align/trim-{trim}/align_qry_{hap}.bed.gz', ASM_TABLE, config,
            trim=('none', 'qry', 'qryref')
        ),
        bed_depth=lambda wildcards: pavlib.pipeline.expand_pattern(
            'results/{asm_name}/align/trim-{trim}/depth_ref_{hap}.bed.gz', ASM_TABLE, config,
            trim=('none', 'qry', 'qryref')
        )

# Create a depth BED file for alignments.
rule align_depth_bed:
    input:
        bed='results/{asm_name}/align/trim-{trim}/align_qry_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/align/trim-{trim}/depth_ref_{hap}.bed.gz'
    run:

        pavlib.align.util.align_bed_to_depth_bed(
            pd.read_csv(input.bed, sep='\t', dtype={'#CHROM': str, 'QRY_ID': str}),
            svpoplib.ref.get_df_fai(REF_FAI)
        ).to_csv(
            output.bed, sep='\t', index=False, compression='gzip'
        )

# Cut alignment overlaps in reference coordinates
rule align_trim_qryref:
    input:
        bed='results/{asm_name}/align/trim-qry/align_qry_{hap}.bed.gz',
        qry_fai='temp/{asm_name}/align/query/query_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/trim-qryref/align_qry_{hap}.bed.gz'
    params:
        align_score=lambda wildcards: get_config('align_score_model', wildcards),
        redundant_callset=lambda wildcards: get_config('redundant_callset', wildcards)
    run:

        # Trim alignments
        df = pavlib.align.trim.trim_alignments(
            pd.read_csv(input.bed, sep='\t', dtype={'#CHROM': str, 'QRY_ID': str}),  # Untrimmed alignment BED
            input.qry_fai,  # Path to query FASTA FAI
            match_qry=params.redundant_callset,  # Redundant callset, trim reference space only for records with matching IDs
            mode='ref',
            score_model=params.align_score
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')

# Cut alignment overlaps in query coordinates
rule align_trim_qry:
    input:
        bed='results/{asm_name}/align/trim-none/align_qry_{hap}.bed.gz',
        qry_fai='temp/{asm_name}/align/query/query_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/trim-qry/align_qry_{hap}.bed.gz'
    params:
        align_score=lambda wildcards: get_config('align_score_model', wildcards)
    run:

        # Trim alignments
        df = pavlib.align.trim.trim_alignments(
            pd.read_csv(input.bed, sep='\t', dtype={'#CHROM': str, 'QRY_ID': str}),  # Untrimmed alignment BED
            input.qry_fai,  # Path to query FASTA FAI
            mode='qry',
            score_model=params.align_score
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# Get alignment BED for one part (one aligned cell or split BAM) in one assembly.
rule align_get_bed:
    input:
        sam='temp/{asm_name}/align/trim-none/align_qry_{hap}.sam.gz',
        qry_fai='temp/{asm_name}/align/query/query_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/trim-none/align_qry_{hap}.bed.gz',
        align_head='results/{asm_name}/align/trim-none/align_qry_{hap}.headers.gz'
    params:
        align_score=lambda wildcards: get_config('align_score_model', wildcards)
    run:

        # Write an empty file if SAM is emtpy
        if os.stat(input.sam).st_size == 0:

            pd.DataFrame(
                [],
                columns=[
                    '#CHROM', 'POS', 'END',
                    'TIER', 'INDEX',
                    'QRY_ID', 'QRY_POS', 'QRY_END',
                    'RG', 'AO',
                    'MAPQ',
                    'REV', 'FLAGS', 'HAP',
                    'CIGAR', 'SCORE',
                    'TRIM_REF_L', 'TRIM_REF_R', 'TRIM_QRY_L', 'TRIM_QRY_R'
                ]
            ).to_csv(
                output.bed, sep='\t', index=False, compression='gzip'
            )

            with open(output.align_head, 'w') as out_file:
                pass

            return

        # Read FAI
        df_qry_fai = svpoplib.ref.get_df_fai(input.qry_fai)
        df_qry_fai.index = df_qry_fai.index.astype(str)

        # Read alignments as a BED file.
        df = pavlib.align.util.get_align_bed(input.sam, df_qry_fai, wildcards.hap, score_model=params.align_score)

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

        # Add trimming fields
        df['TRIM_REF_L'] = 0
        df['TRIM_REF_R'] = 0
        df['TRIM_QRY_L'] = 0
        df['TRIM_QRY_R'] = 0

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# Map query as SAM. Pull read information from the SAM before sorting and writing CRAM since tool tend to change
# "=X" to "M" in the CIGAR.
rule align_map:
    input:
        ref_fa='data/ref/ref.fa.gz',
        seq_files=lambda wildcards: pavlib.config.get_aligner_input(wildcards, config, ASM_TABLE)
    output:
        sam=temp('temp/{asm_name}/align/trim-none/align_qry_{hap}.sam.gz')
    params:
        aligner=lambda wildcards: pavlib.config.get_aligner(wildcards, config, ASM_TABLE),
        aligner_params=lambda wildcards: pavlib.config.get_aligner_params(wildcards, config, ASM_TABLE)
    threads: 4
    run:

        # Get alignment command
        if params.aligner == 'minimap2':
            align_cmd = (
                f"""minimap2 """
                    f"""{params.aligner_params} """
                    f"""--secondary=no -a -t {threads} --eqx -Y """
                    f"""{input.ref_fa} {input.seq_files[0]}"""
            )

        elif params.aligner == 'lra':
            align_cmd = (
                f"""lra align {input.ref_fa} {input.seq_files[0]} -CONTIG -p s -t {threads}"""
                f"""{params.aligner_params}"""
            )

        else:
            raise RuntimeError(f'Unknown alignment program (aligner parameter): {params.aligner}')

        # Run alignment
        if os.stat(input.seq_files[0]).st_size > 0:

            # Run alignment
            print(f'Aligning {wildcards.asm_name}-{wildcards.hap}: {align_cmd}', flush=True)

            with Bio.bgzf.BgzfWriter(output.sam, 'wt') as out_file:
                for line in shell(align_cmd, iterable=True):
                    if not line.startswith('@'):
                        line = line.split('\t')
                        line[9] = '*'
                        line[10] = '*'
                        line = '\t'.join(line)

                    out_file.write(line)
                    out_file.write('\n')

        else:
            # Write an empty file if input is empty
            with open(output.sam, 'w') as out_file:
                pass

# Uncompress query sequences for aligners that cannot read gzipped FASTAs.
rule align_uncompress_qry:
    input:
        fa='temp/{asm_name}/align/query/query_{hap}.fa.gz'
    output:
        fa='temp/{asm_name}/align/query/query_{hap}.fa'
    run:

        if os.stat(input.fa).st_size > 0:
            shell(
                """zcat {input.fa} > {output.fa}"""
            )
        else:
            with open(output.fa, 'w') as out_file:
                pass

# Get FASTA files.
rule align_get_qry_fa:
    input:
        fa=lambda wildcards: pavlib.pipeline.get_rule_input_list(wildcards.asm_name, wildcards.hap, ASM_TABLE, get_override_config(wildcards.asm_name))
    output:
        fa=temp('temp/{asm_name}/align/query/query_{hap}.fa.gz'),
        fai=temp('temp/{asm_name}/align/query/query_{hap}.fa.gz.fai'),
        gzi=temp('temp/{asm_name}/align/query/query_{hap}.fa.gz.gzi')
    run:

        # Get input files
        # noinspection PyUnresolvedReferences
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
            os.symlink(os.path.abspath(input_tuples[0][0]), output.fa)
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
# Export alignments (optional feature)
#

def _align_export_all(wildcards):

    if 'trim' in config:
        trim_set = set(config['trim'].strip().split(','))
    else:
        trim_set = {'qryref'}

    if 'export_fmt' in config:
        ext_set = set(config['export_fmt'].strip().split(','))
    else:
        ext_set = {'cram'}

    ext_set = set([ext if ext != 'sam' else 'sam.gz' for ext in ext_set])

    if 'asm_name' in config:
        asm_set = set(config['asm_name'].strip().split(','))
    else:
        asm_set = None

    if 'hap' in config:
        hap_set = set(config['hap'].strip().split(','))
    else:
        hap_set = None

    return pavlib.pipeline.expand_pattern(
        'results/{asm_name}/align/export/pav_align_trim-{trim}_{hap}.{ext}',
        ASM_TABLE, config,
        asm_name=asm_set, hap=hap_set, trim=trim_set, ext=ext_set
    )

# Get CRAM files
localrules: align_export_all

rule align_export_all:
    input:
        cram=_align_export_all


# Reconstruct CRAM from alignment BED files after trimming redundantly mapped bases (post-cut).
rule align_export:
    input:
        bed='results/{asm_name}/align/trim-{trim}/align_qry_{hap}.bed.gz',
        fa='temp/{asm_name}/align/query/query_{hap}.fa.gz',
        align_head='results/{asm_name}/align/trim-none/align_qry_{hap}.headers.gz',
        ref_fa='data/ref/ref.fa.gz'
    output:
        align='results/{asm_name}/align/export/pav_align_trim-{trim}_{hap}.{ext}'
    params:
        sam_tag=lambda wildcards: fr'@PG\tID:PAV-{wildcards.trim}\tPN:PAV\tVN:{pavlib.const.get_version_string()}\tDS:PAV Alignment trimming {pavlib.align.util.TRIM_DESC[wildcards.trim]}'
    run:

        if wildcards.ext == 'cram':
            out_fmt = 'CRAM'
            do_bgzip = False
            do_index = True
            do_tabix = False

        elif wildcards.ext == 'bam':
            out_fmt = 'BAM'
            do_bgzip = False
            do_index = True
            do_tabix = False

        elif wildcards.ext == 'sam.gz':
            out_fmt = 'SAM'
            do_bgzip = True
            do_index = False
            do_tabix = True

        else:
            raise RuntimeError(f'Unknown output format extension: {wildcards.ext}: (Allowed: "cram", "bam", "sam.gz")')

        # Export

        if not do_bgzip:
            shell(
                """python3 {PIPELINE_DIR}/scripts/reconstruct_sam.py """
                    """--bed {input.bed} --fasta {input.fa} --headers {input.align_head} --tag "{params.sam_tag}" | """
                """samtools view -T {input.ref_fa} -O {out_fmt} -o {output.align}"""
            )
        else:
            shell(
                """python3 {PIPELINE_DIR}/scripts/reconstruct_sam.py """
                    """--bed {input.bed} --fasta {input.fa} --headers {input.align_head} --tag "{params.sam_tag}" | """
                """samtools view -T {input.ref_fa} -O {out_fmt} | """
                """bgzip > {output.align}"""
            )

        # Index
        if do_index:
            shell(
                """samtools index {output.align}"""
            )

        if do_tabix:
            shell(
                """tabix {output.align}"""
            )
