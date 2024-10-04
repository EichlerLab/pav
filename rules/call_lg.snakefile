"""
Call alignment-truncating events (large SVs).
"""

import os
import pandas as pd

import pavlib

global REF_FA
global get_config

global expand
global shell
global temp
global ASM_TABLE

# Call all large SVs
localrules: call_lg_all

rule call_lg_all:
    input:
        bed=lambda wildcards: pavlib.pipeline.expand_pattern(
            'temp/{asm_name}/lg_sv/t{tier}/svindel_ins_{hap}.bed.gz', ASM_TABLE, config,
            tier=(0, 1, 2)
        )


# Call alignment-truncating SVs.
rule call_lg_discover:
    input:
        bed_qry='results/{asm_name}/align/t{tier}_qry/align_qry_{hap}.bed.gz',
        bed_qryref='results/{asm_name}/align/t{tier}_qryref/align_qry_{hap}.bed.gz',
        fa_qry='temp/{asm_name}/align/query/query_{hap}.fa.gz',
        fai_qry='temp/{asm_name}/align/query/query_{hap}.fa.gz.fai',
        fa_ref='data/ref/ref.fa.gz',
        fai_ref='data/ref/ref.fa.gz.fai'
    output:
        bed_ins=temp('temp/{asm_name}/lg_sv/t{tier}/svindel_ins_{hap}.bed.gz'),
        bed_del=temp('temp/{asm_name}/lg_sv/t{tier}/svindel_del_{hap}.bed.gz'),
        bed_inv=temp('temp/{asm_name}/lg_sv/t{tier}/sv_inv_{hap}.bed.gz'),
        bed_cpx=temp('temp/{asm_name}/lg_sv/t{tier}/sv_cpx_{hap}.bed.gz'),
        bed_cpx_seg=temp('temp/{asm_name}/lg_sv/t{tier}/segment_cpx_{hap}.bed.gz'),
        bed_cpx_ref=temp('temp/{asm_name}/lg_sv/t{tier}/reftrace_cpx_{hap}.bed.gz')
    params:
        align_score=lambda wildcards: get_config(wildcards, 'align_score_model', pavlib.align.score.DEFAULT_ALIGN_SCORE_MODEL),
        min_anchor_score=lambda wildcards: get_config(wildcards, 'min_anchor_score', "1000bp"),
        lg_dot_graph=lambda wildcards: get_config(wildcards, 'lg_dot_graph', 'False')
    #threads: 12
    run:

        # Set graph file output
        if pavlib.util.as_bool(params.lg_dot_graph):
            dot_basename = f'temp/{wildcards.asm_name}/lg_sv/t{wildcards.tier}/graph/lgsv_{wildcards.hap}_'
            dot_dirname = os.path.dirname(dot_basename)
            os.makedirs(dot_dirname, exist_ok=True)
        else:
            dot_basename = None

        # Get score model
        score_model = pavlib.align.score.get_score_model(params.align_score)

        # Get minimum anchor score
        if isinstance(params.min_anchor_score, str):
            min_anchor_score_str = params.min_anchor_score.strip()

            if len(min_anchor_score_str) == 0:
                raise RuntimeError('Parameter "min_anchor_score" is empty')

            if min_anchor_score_str.lower().endswith('bp'):
                min_anchor_score_bp = abs(float(min_anchor_score_str[:-2]))

                if min_anchor_score_bp == 0:
                    raise RuntimeError(f'Parameter "min_anchor_score" is zero: {params.min_anchor_score}')

                min_anchor_score = score_model.match(min_anchor_score_bp)

            else:
                try:
                    min_anchor_score = abs(float(params.min_anchor_score))
                except ValueError:
                    raise RuntimeError(f'Parameter "min_anchor_score" is a string that does not represent a numeric value: type={min_anchor_score}')

        else:
            try:
                # noinspection PyTypeChecker
                min_anchor_score = float(params.min_anchor_score)
            except ValueError:
                raise RuntimeError(f'Parameter "min_anchor_score" is not a string or numeric: type={type(min_anchor_score)}')

            if min_anchor_score < 0:
                raise RuntimeError(f'Parameter "min_anchor_score" is negative: {min_anchor_score}')

        # Read alignments - Trim QRY
        df_align_qry = pd.read_csv(
            input.bed_qry,
            sep='\t',
            dtype={'#CHROM': str, 'QRY_ID': str}
        )

        df_align_qry.sort_values(['QRY_ID', 'QRY_POS', 'QRY_END'], inplace=True)
        df_align_qry.reset_index(inplace=True, drop=True)

        # Read alignments - Trim QRY/REF
        df_align_qryref = pd.read_csv(
            input.bed_qryref,
            sep='\t',
            index_col='INDEX',
            dtype={'#CHROM': str, 'QRY_ID': str}
        )

        # Set caller resources
        caller_resources = pavlib.lgsv.util.CallerResources(
            df_align_qry, df_align_qryref,
            input.fa_qry, input.fa_ref,
            wildcards.hap, score_model
        )

        # Call
        lgsv_list = pavlib.lgsv.call.call_from_align(caller_resources, min_anchor_score=min_anchor_score, dot_basename=dot_basename)

        # Create tables
        df_list = {
            'INS': list(),
            'DEL': list(),
            'INV': list(),
            'CPX': list()
        }

        for var in lgsv_list:
            row = var.row()

            if row['SVTYPE'] not in df_list.keys():
                raise RuntimeError(f'Unexpected SVTYPE: "{row["SVTYPE"]}"')

            df_list[row['SVTYPE']].append(row)

        if len(df_list['INS']) > 0:
            df_ins = pd.concat(df_list['INS'], axis=1).T
        else:
            df_ins = pd.DataFrame([], columns=pavlib.lgsv.variant.InsertionVariant(None, None).row().index)

        if len(df_list['DEL']) > 0:
            df_del = pd.concat(df_list['DEL'], axis=1).T
        else:
            df_del = pd.DataFrame([], columns=pavlib.lgsv.variant.DeletionVariant(None, None).row().index)

        if len(df_list['INV']) > 0:
            df_inv = pd.concat(df_list['INV'], axis=1).T
        else:
            df_inv = pd.DataFrame([], columns=pavlib.lgsv.variant.InversionVariant(None, None).row().index)

        if len(df_list['CPX']) > 0:
            df_cpx = pd.concat(df_list['CPX'], axis=1).T
        else:
            df_cpx = pd.DataFrame([], columns=pavlib.lgsv.variant.ComplexVariant(None, None).row().index)

        df_ins.sort_values(['#CHROM', 'POS', 'END', 'ID', 'QRY_REGION'], inplace=True)
        df_del.sort_values(['#CHROM', 'POS', 'END', 'ID', 'QRY_REGION'], inplace=True)
        df_inv.sort_values(['#CHROM', 'POS', 'END', 'ID', 'QRY_REGION'], inplace=True)
        df_cpx.sort_values(['#CHROM', 'POS', 'END', 'ID', 'QRY_REGION'], inplace=True)

        # Write variant tables
        df_ins.to_csv(output.bed_ins, sep='\t', index=False, compression='gzip')
        df_del.to_csv(output.bed_del, sep='\t', index=False, compression='gzip')
        df_inv.to_csv(output.bed_inv, sep='\t', index=False, compression='gzip')
        df_cpx.to_csv(output.bed_cpx, sep='\t', index=False, compression='gzip')

        # Write segment and reference trace tables
        df_segment_list = list()
        df_reftrace_list = list()

        for var in lgsv_list:
            if var.svtype != 'CPX':
                continue

            df_segment = var.interval.df_segment.copy()
            df_segment.insert(3, 'ID', var.variant_id)

            df_reftrace = var.df_ref_trace.copy()
            df_reftrace.insert(3, 'ID', var.variant_id)

            df_segment_list.append(df_segment)
            df_reftrace_list.append(df_reftrace)

        if len(df_segment_list) > 0:
            df_segment = pd.concat(df_segment_list, axis=0)
        else:
            df_segment = pd.DataFrame([], columns=(
                pavlib.lgsv.interval.SEGMENT_TABLE_FIELDS[:3] + ['ID'] + pavlib.lgsv.interval.SEGMENT_TABLE_FIELDS[3:]
            ))

        if len(df_reftrace_list) > 0:
            df_reftrace = pd.concat(df_reftrace_list, axis=0)
        else:
            df_reftrace = pd.DataFrame([], columns=(
                pavlib.lgsv.variant.REF_TRACE_COLUMNS[:3] + ['ID'] + pavlib.lgsv.variant.REF_TRACE_COLUMNS[3:]
            ))

        df_segment.to_csv(output.bed_cpx_seg, sep='\t', index=False, compression='gzip')
        df_reftrace.to_csv(output.bed_cpx_ref, sep='\t', index=False, compression='gzip')

        # Compress graph dot files
        if dot_basename is not None:
            dot_tar_filename = f'results/{wildcards.asm_name}/lg_sv/lgsv_graph_{wildcards.hap}_t{wildcards.tier}.tar'

            os.makedirs(os.path.dirname(dot_tar_filename), exist_ok=True)

            shell(
                """tar -cf {dot_tar_filename} {dot_dirname}/*"""
            )
