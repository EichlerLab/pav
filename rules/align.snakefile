"""
Process alignments and alignment tiling paths.
"""


import gzip
import numpy as np
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
        bed_depth=lambda wildcards: pavlib.pipeline.expand_pattern(
            'results/{asm_name}/align/t{tier}_{trim}/align_qry_{hap}.bed.gz', ASM_TABLE,
            trim=('none', 'qry', 'qryref'), tier=('1', '2')
        ),
        bed_depth_t0=lambda wildcards: pavlib.pipeline.expand_pattern(
            'results/{asm_name}/align/t0_{trim}/align_qry_{hap}.bed.gz', ASM_TABLE,
            trim=('qry', 'qryref')
        )
        # bed_depth=lambda wildcards: pavlib.pipeline.expand_pattern(
        #     'results/{asm_name}/align/t{tier}_{trim}/depth_ref_{hap}.bed.gz', ASM_TABLE,
        #     trim=('none', 'qry', 'qryref'), tier=('1', '2')
        # )

# Create a depth BED file for alignments.
rule align_depth_bed:
    input:
        bed='results/{asm_name}/align/t{tier}_{trim}/align_qry_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/align/t{tier}_{trim}/depth_ref_{hap}.bed.gz'
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
        bed='results/{asm_name}/align/t{tier}_qry/align_qry_{hap}.bed.gz',
        qry_fai='temp/{asm_name}/align/query/query_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/t{tier}_qryref/align_qry_{hap}.bed.gz'
    params:
        align_score=lambda wildcards: get_config(wildcards,'align_score_model', pavlib.align.score.DEFAULT_ALIGN_SCORE_MODEL),
        redundant_callset=lambda wildcards: pavlib.util.as_bool(get_config(wildcards, 'redundant_callset', False))
        # min_trim_qry_len=lambda wildcards: int(get_config(wildcards, 'min_trim_qry_len', 1000)),  # Minimum aligned length
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


# Create tier 1 alignments: Truncate tier 2 alignments to eliminate tig-space overlaps with tier 1
rule align_join_tiers:
    input:
        bed_t1='results/{asm_name}/align/t1_qry/align_qry_{hap}.bed.gz',
        bed_t2='results/{asm_name}/align/t2_qry/align_qry_{hap}.bed.gz',
        fai='temp/{asm_name}/align/query/query_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/t0_qry/align_qry_{hap}.bed.gz'
    params:
        align_score=lambda wildcards: get_config(wildcards,'align_score_model', pavlib.align.score.DEFAULT_ALIGN_SCORE_MODEL)
    run:

        score_model = pavlib.align.score.get_score_model(params.align_score)

        # Read
        df1 = pd.read_csv(input.bed_t1, sep='\t', dtype={'#CHROM': str, 'QRY_ID': str})
        df1.sort_values(['QRY_ID', 'QRY_POS', 'QRY_END'], inplace=True)

        df2 = pd.read_csv(input.bed_t2, sep='\t', dtype={'#CHROM': str, 'QRY_ID': str})
        df2.sort_values(['QRY_ID', 'QRY_POS', 'QRY_END'], inplace=True)

        df_fai = svpoplib.ref.get_df_fai(input.fai)

        # Adjust indexes for merge
        if 'INDEX_ORG' not in df1.columns:
            df1['INDEX_ORG'] = df1['INDEX']

        if 'INDEX_ORG' not in df2.columns:
            df2['INDEX_ORG'] = df2['INDEX']

        df2['INDEX'] += int(10**np.ceil(np.log10(max(df1['INDEX']))))

        df2['RETAIN'] = True  # If False after trimming, record is discarded

        # Set index
        if len(set(df1['INDEX'])) != len(df1):
            raise RuntimeError(f'Duplicate values in INDEX in tier 1 table: {input.bed_t1}')

        if len(set(df2['INDEX'])) != len(df2):
            raise RuntimeError(f'Duplicate values in INDEX in tier 2 table: {input.bed_t2}')

        df1.set_index('INDEX', inplace=True, drop=False)
        df1.index.name = 'INDEX_T1'

        df2.set_index('INDEX', inplace=True, drop=False)
        df2.index.name = 'INDEX_T2'

        # Get records from tier 2 that cover sequence not aligned in tier 1
        for qry_id in set(df2['QRY_ID']):

            # Step through each alignment record
            index1_list = list(df1.loc[df1['QRY_ID'] == qry_id].index)[::-1]
            index2_list = list(df2.loc[df2['QRY_ID'] == qry_id].index)[::-1]

            index1 = index1_list.pop() if index1_list else None
            index2 = index2_list.pop() if index2_list else None

            while index1 is not None and index2 is not None:

                if not df2.loc[index2, 'RETAIN']:
                    raise RuntimeError('Program Bug: Encountered dropped record in traversal')

                # Skip tier2 records inside a tier1 record
                if df2.loc[index2, 'QRY_POS'] >= df1.loc[index1, 'QRY_POS'] and \
                        df2.loc[index2, 'QRY_END'] <= df1.loc[index1, 'QRY_END']:
                    # |-------- 1 --------|
                    #     |---- 2 ----|

                    df2.loc[index2, 'RETAIN'] = False

                    index2 = index2_list.pop() if index2_list else None
                    continue

                # Skip tier2 records spanning a tier1 record
                # Assume tier1 is a better representation than tier2, even if shorter. Tier1 alignments should be
                # better in general, but may be room for improvement in the future.
                if df2.loc[index2, 'QRY_POS'] <= df1.loc[index1, 'QRY_POS'] and \
                        df2.loc[index2, 'QRY_END'] >= df1.loc[index1, 'QRY_END']:
                    #     |---- 1 ----|
                    # |-------- 2 --------|

                    df2.loc[index2, 'RETAIN'] = False

                    index2 = index2_list.pop() if index2_list else None
                    continue

                # Skip non-overlapping records
                if df2.loc[index2, 'QRY_END'] <= df1.loc[index1, 'QRY_POS']:
                    #                 |---- 1 ----|
                    # |---- 2 ----|

                    index2 = index2_list.pop() if index2_list else None
                    continue

                if df2.loc[index2, 'QRY_POS'] >= df1.loc[index1, 'QRY_END']:
                    # |---- 1 ----|
                    #                 |---- 2 ----|

                    index1 = index1_list.pop() if index1_list else None
                    continue

                # Process overlap
                record = None
                overlap_bp = None
                trunc_side = None

                if df2.loc[index2, 'QRY_POS'] < df1.loc[index1, 'QRY_POS']:
                    #           |---- 1 ----|
                    # |---- 2 ----|
                    overlap_bp = df2.loc[index2, 'QRY_END'] - df1.loc[index1, 'QRY_POS']
                    trunc_side = 'r'

                    record = df2.loc[index2].copy()

                    try:
                        record = pavlib.align.trim.truncate_alignment_record(
                            record, overlap_bp, trunc_side, score_model=score_model
                        )
                    except RuntimeError as e:
                        raise RuntimeError(f'Alignment overlap truncation error in {qry_id}: tier 1 index {index1}, tier 2 index {index2}: overlap_bp={overlap_bp}, trunc_side={trunc_side}: {e}') from e

                elif df2.loc[index2, 'QRY_END'] > df1.loc[index1, 'QRY_END']:
                    # |---- 1 ----|
                    #           |---- 2 ----|
                    overlap_bp = df1.loc[index1, 'QRY_END'] - df2.loc[index2, 'QRY_POS']
                    trunc_side = 'l'

                    record = df2.loc[index2].copy()

                    try:
                        record = pavlib.align.trim.truncate_alignment_record(
                            record, overlap_bp, trunc_side, score_model=score_model
                        )
                    except RuntimeError as e:
                        raise RuntimeError(f'Alignment overlap truncation error in {qry_id}: tier 1 index {index1}, tier 2 index {index2}: overlap_bp={overlap_bp}, trunc_side={trunc_side}: {e}') from e
                else:
                    raise RuntimeError(f'Overlap logic error in {qry_id}: tier 1 index {index1}, tier 2 index {index2}')

                # Check truncation
                if record is None:
                    raise RuntimeError(f'Alignment overlap truncation error in {qry_id}: tier 1 index {index1}, tier 2 index {index2}: Full record removed by truncation')

                try:
                    pavlib.align.util.check_record(record, df_fai)
                except RuntimeError as e:
                    raise RuntimeError(f'Alignment overlap truncation error in {qry_id}: tier 1 index {index1}, tier 2 index {index2}: overlap_bp={overlap_bp}, trunc_side={trunc_side}: {e}') from e

                # Modify record
                df2.loc[index2] = record

                # Next record
                if trunc_side == 'l':
                    index1 = index1_list.pop() if index1_list else None
                elif trunc_side == 'r':
                    index2 = index2_list.pop() if index2_list else None
                else:
                    raise RuntimeError(f'Program Bug: Alignment overlap truncation error in {qry_id}: tier 1 index {index1}, tier 2 index {index2}: Unknown truncation side: {trunc_side}')

        # Merge tables and write
        df2 = df2.loc[df2['RETAIN']]
        del df2['RETAIN']

        df = pd.concat([df1, df2], axis=0).sort_values(['#CHROM', 'POS', 'END', 'QRY_ID', 'QRY_POS', 'QRY_END'])

        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


# Cut alignment overlaps in query coordinates
rule align_trim_qry:
    input:
        bed='results/{asm_name}/align/t{tier}_none/align_qry_{hap}.bed.gz',
        qry_fai='temp/{asm_name}/align/query/query_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/t{tier}_qry/align_qry_{hap}.bed.gz'
    wildcard_constraints:
        tier=r'[1-9]+[0-9]*'  # Tier 0 not allowed
    params:
        align_score=lambda wildcards: get_config(wildcards, 'align_score_model', pavlib.align.score.DEFAULT_ALIGN_SCORE_MODEL)
    # params:
    #     min_trim_qry_len=lambda wildcards: int(get_config(wildcards, 'min_trim_qry_len', 1000))  # Minimum aligned length
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
        sam='temp/{asm_name}/align/t{tier}_none/align_qry_{hap}.sam.gz',
        qry_fai='temp/{asm_name}/align/query/query_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/t{tier}_none/align_qry_{hap}.bed.gz',
        align_head='results/{asm_name}/align/t{tier}_none/align_qry_{hap}.headers.gz'
    params:
        align_score=lambda wildcards: get_config(wildcards, 'align_score_model', pavlib.align.score.DEFAULT_ALIGN_SCORE_MODEL)
    wildcard_constraints:
        tier=r'[12]'
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
                    'CALL_BATCH',
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
        df = pavlib.align.util.get_align_bed(input.sam, df_qry_fai, wildcards.hap, tier=int(wildcards.tier), score_model=params.align_score)

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


# Map query as SAM. Pull read information from the SAM before sorting and writing CRAM since tool tend to change
# "=X" to "M" in the CIGAR.
rule align_map:
    input:
        ref_fa='data/ref/ref.fa.gz',
        fa=lambda wildcards: 'temp/{asm_name}/align/query/query_{hap}.fa.gz' if pavlib.align.params.get_aligner(wildcards.tier, config) != 'lra' else 'temp/{asm_name}/align/query/query_{hap}.fa',
        fai=lambda wildcards: 'temp/{asm_name}/align/query/query_{hap}.fa.gz.fai' if pavlib.align.params.get_aligner(wildcards.tier, config) != 'lra' else [],
        gzi=lambda wildcards: 'temp/{asm_name}/align/query/query_{hap}.fa.gz.gzi' if pavlib.align.params.get_aligner(wildcards.tier, config) != 'lra' else [],
        gli=lambda wildcards: 'data/ref/ref.fa.gz.gli' if pavlib.align.params.get_aligner(wildcards.tier, config) == 'lra' else [],
        mmi=lambda wildcards: 'data/ref/ref.fa.gz.mms' if pavlib.align.params.get_aligner(wildcards.tier, config) == 'lra' else []
    output:
        sam=temp('temp/{asm_name}/align/t{tier}_none/align_qry_{hap}.sam.gz')
    wildcard_constraints:
        tier=r'[12]'
    threads: 4
    run:

        # Get alignment command
        aligner = pavlib.align.params.get_aligner(wildcards.tier, config)
        aligner_params = pavlib.align.params.get_aligner_params(wildcards.tier, config, aligner)

        if aligner == 'minimap2':
            align_cmd = (
                f"""minimap2 """
                    f"""{aligner_params} """
                    f"""--secondary=no -a -t {threads} --eqx -Y """
                    f"""{input.ref_fa} {input.fa}"""
            )

        elif aligner == 'lra':
            align_cmd = (
                f"""lra align {input.ref_fa} {input.fa} -CONTIG -p s -t {threads}"""
                f"""{aligner_params}"""
            )

        else:
            raise RuntimeError(f'Unknown alignment program (aligner parameter): {aligner}')

        # Run alignment
        # noinspection PyTypeChecker
        if os.stat(input.fa).st_size > 0:

            # Run alignment
            print(f'Aligning {wildcards.asm_name}-{wildcards.hap}: {align_cmd}', flush=True)

            shell(
                align_cmd + (
                    """ | awk -vOFS="\\t" '($1 !~ /^@/) {{$10 = "*"; $11 = "*"}} {{print}}' | """
                    """gzip > {output.sam}"""
                )
            )

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
        fa=lambda wildcards: pavlib.pipeline.get_rule_input_list(wildcards.asm_name, wildcards.hap, ASM_TABLE, get_config(wildcards))
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
# Export alignments (optional feature)
#

def _align_export_all(wildcards):

    if 'trim' in config:
        trim_set = set(config['trim'].strip().split(','))
    else:
        trim_set = {'qryref'}

    if 'tier' in config:
        tier_set = set(config['tier'].strip().split(','))
    else:
        tier_set = {'0', '1', '2'}

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
        'results/{asm_name}/align/export/pav_align_tier-{tier}_trim-{trim}_{hap}.{ext}',
        ASM_TABLE,
        asm_name=asm_set, hap=hap_set, trim=trim_set, tier=tier_set, ext=ext_set
    )

# Get CRAM files
localrules: align_export_all

rule align_export_all:
    input:
        cram=_align_export_all


# Reconstruct CRAM from alignment BED files after trimming redundantly mapped bases (post-cut).
rule align_export:
    input:
        bed='results/{asm_name}/align/t{tier}_{trim}/align_qry_{hap}.bed.gz',
        fa='temp/{asm_name}/align/query/query_{hap}.fa.gz',
        align_head='results/{asm_name}/align/t{tier}_none/align_qry_{hap}.headers.gz',
        ref_fa='data/ref/ref.fa.gz'
    output:
        align='results/{asm_name}/align/export/pav_align_tier-{tier}_trim-{trim}_{hap}.{ext}'
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
