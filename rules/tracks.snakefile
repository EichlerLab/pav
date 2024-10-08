"""
Prepare UCSC tracks for data.
"""

import matplotlib as mpl
import numpy as np
import os
import pandas as pd

import pavlib
import svpoplib

global ASM_TABLE
global PIPELINE_DIR
global REF_FAI
global get_config
global temp


#
# Definitions
#

VARTYPE_TO_SVTYPE_TUPLE = {
    'snv': ('snv',),
    'sv': ('ins', 'del'),
    'indel': ('ins', 'del')
}

ALIGN_COLORMAP = 'viridis'

VAR_BED_PREMERGE_PATTERN = 'results/{{asm_name}}/bed_hap/{filter}/{{hap}}/{vartype}_{svtype}.bed.gz'

CPX_COLOR_DICT = {
    'INS': '64,64,255',
    'DEL': '255,64,64',
    'INV': '96,255,96',

    'DUP': '255,64,255',
    'TRP': '217,54,217',
    'QUAD': '179,45,179',
    'HDUP': '96,32,96',

    'INVDUP': '82,217,94',
    'INVTRP': '70,179,70',
    'INVQUAD': '55,140,55',
    'INVHDUP': '40,102,40',

    'MIXDUP': '217, 133, 37',
    'MIXTRP': '179, 109, 30',
    'MIXQUAD': '139, 85, 24',
    'MIXHDUP': '102, 63, 18',

    'NML': '51,51,51',

    'UNMAPPED': '115, 115, 115'
}


def _track_get_input_bed(wildcards):
    """
    Get one or more input files for tracks. If "svtype" is "all", collect all relevant input files.

    :param wildcards: Wildcards.

    :return: List of input file(s).
    """

    # Get input file variant type
    if wildcards.vartype in {'sv', 'indel', 'svindel'}:


        if wildcards.svtype in {'ins', 'del'}:
            input_vartype = 'svindel'
            input_svtype = [wildcards.vartype]
        elif wildcards.svtype == 'insdel':
            input_vartype = 'svindel'
            input_svtype = ['ins', 'del']
        elif wildcards.svtype == 'inv':
            input_vartype = 'sv'
            input_svtype = ['inv']
        else:
            raise RuntimeError(f'Unknown svtype {wildcards.svtype} for variant type {wildcards.vartype} (expected "ins", "del", or "insdel" - "inv" allowed for vartype "sv")')

        if 'inv' in input_svtype and input_vartype not in {'sv', 'svindel'}:
            raise RuntimeError(f'Bad svtype {wildcards.svtype} for variant type {wildcards.vartype}: vartype must include SVs to output inversions')

    elif wildcards.vartype == 'snv':
        if wildcards.svtype != 'snv':
            raise RuntimeError(f'Unknown svtype {wildcards.svtype} for variant type {wildcards.vartype} (expected "snv")')

        input_vartype = 'snv'
        input_svtype = ['snv']

    else:
        raise RuntimeError(f'Unrecognized variant type: {wildcards.vartype}')

    # Get filter list
    if wildcards.filter == 'pass':
        input_filter = ['pass']
    elif wildcards.filter == 'fail':
        input_filter = ['fail']
    elif wildcards.filter == 'all':
        input_filter = ['pass', 'fail']
    else:
        raise RuntimeError(f'Unknown input filter (wildcards.filter): "{wildcards.filter}"')

    # Return list of input files
    return [
        VAR_BED_PREMERGE_PATTERN.format(vartype=input_vartype, svtype=svtype, filter=filter)
            for svtype in input_svtype for filter in input_filter
    ]


#
# All
#

# All tracks
rule tracks_all:
    input:
        bb_call=lambda wildcards: pavlib.pipeline.expand_pattern(
            'tracks/{asm_name}/variant/pre_merge/pass/{varsvtype}_{hap}.bb', ASM_TABLE, config,
            varsvtype=['sv_insdel', 'sv_inv', 'indel_insdel', 'snv_snv']
        ),
        bb_align=lambda wildcards: pavlib.pipeline.expand_pattern(
            'tracks/{asm_name}/align/align_qry_trim-{trim}.bb', ASM_TABLE, config,
            trim=('none', 'qry', 'qryref')
        ),
        bb_invflag=lambda wildcards: pavlib.pipeline.expand_pattern(
            'tracks/{asm_name}/inv_flag/inv_flag.bb', ASM_TABLE, config
        )

rule tracks_bb:
    input:
        bed='temp/{asm_name}/tracks/{subdir}/{filename}.bed',
        asfile='temp/{asm_name}/tracks/{subdir}/{filename}.as',
        fai=REF_FAI
    output:
        bb='tracks/{asm_name}/{subdir}/{filename}.bb'
    shell:
        """bedToBigBed -tab -as={input.asfile} -type=bed9+ {input.bed} {input.fai} {output.bb}"""

#
# Variant calls
#

rule tracks_hap_call_all:
    input:
        bb=lambda wildcards: pavlib.pipeline.expand_pattern(
            'tracks/{asm_name}/variant/pre_merge/pass/{varsvtype}_{hap}.bb', ASM_TABLE, config,
            varsvtype=['sv_insdel', 'sv_inv', 'indel_insdel', 'snv_snv']
        )

# # BigBed for one variant set.
# rule tracks_hap_call_bb:
#     input:
#         bed='temp/{asm_name}/tracks/bed_pre_merge/{filter}/{vartype}_{svtype}_{hap}.bed',
#         asfile='temp/{asm_name}/tracks/bed_pre_merge/{filter}/{vartype}_{svtype}_{hap}.as',
#         fai=REF_FAI
#     output:
#         bb='tracks/{asm_name}/variant/pre_merge/{filter}/{vartype}_{svtype}_{hap}.bb'
#     shell:
#         """bedToBigBed -tab -as={input.asfile} -type=bed9+ {input.bed} {input.fai} {output.bb}"""

# Tracks for one variant set.
rule tracks_hap_call:
    input:
        bed=_track_get_input_bed,
        fai=REF_FAI
    output:
        bed=temp('temp/{asm_name}/tracks/bed_pre_merge/{filter}/{vartype}_{svtype}_{hap}.bed'),
        asfile=temp('temp/{asm_name}/tracks/bed_pre_merge/{filter}/{vartype}_{svtype}_{hap}.as')
    wildcard_constraints:
        filter='pass|fail|all'
    run:

        if wildcards.filter != 'pass':
            raise NotImplementedError(f'Tracks containing non-PASS variants is not yet supported: {wildcards.filter}')

        field_table_file_name = os.path.join(PIPELINE_DIR, 'files/tracks/variant_track_fields.tsv')

        # Read variants
        df = pd.concat(
            [pd.read_csv(file_name, sep='\t', dtype={'#CHROM': str}, low_memory=False) for file_name in input.bed],
            axis=0
        ).reset_index(drop=True)

        if wildcards.vartype == 'sv':
            df = df.loc[df['SVLEN'] >= 50].copy()
        elif wildcards.vartype == 'indel':
            df = df.loc[df['SVLEN'] < 50].copy()

        df.sort_values(['#CHROM', 'POS', 'END'], inplace=True)

        # Select table columns
        for col in ('QRY_ID', 'QRY_POS', 'QRY_END', 'QRY_MAPPED', 'QRY_STRAND', 'SEQ'):
            if col in df.columns:
                del(df[col])

        # Read FAI and table columns
        df_fai = svpoplib.ref.get_df_fai(input.fai)

        # Filter columns that have track annotations
        field_set = set(
            pd.read_csv(
                field_table_file_name,
                sep='\t', header=0
            )['FIELD']
        )

        df = df.loc[:, [col for col in df.columns if col in field_set]]

        # Make BigBed
        track_name = 'VariantTable'
        track_description = '{asm_name} - {vartype}-{svtype} - {hap}'.format(**wildcards)

        svpoplib.tracks.variant.make_bb_track(df, df_fai, output.bed, output.asfile, track_name, track_description, field_table_file_name)


#
# Alignments
#

# Generate all alignment tracks
rule tracks_align_all:
    input:
        bb=lambda wildcards: pavlib.pipeline.expand_pattern(
            'tracks/{asm_name}/align/align_qry_trim-{trim}.bb', ASM_TABLE, config,
            trim=('none', 'qry', 'qryref')
        )

# # Alignment track BED to BigBed.
# rule tracks_align_bb:
#     input:
#         bed='temp/{asm_name}/tracks/align/align_qry_trim-{trim}.bed',
#         asfile='temp/{asm_name}/tracks/align/align_qry_trim-{trim}.as'
#     output:
#         bb='tracks/{asm_name}/align/align_qry_trim-{trim}.bb'
#     shell:
#         """bedToBigBed -tab -as={input.asfile} -type=bed9+ {input.bed} {REF_FAI} {output.bb}"""

# Alignment tracks.
rule tracks_align:
    input:
        bed=lambda wildcards: [
            f'results/{wildcards.asm_name}/align/trim-{wildcards.trim}/align_qry_{hap}.bed.gz'
                for hap in pavlib.pipeline.get_hap_list(wildcards.asm_name, ASM_TABLE)
        ]
    output:
        bed=temp('temp/{asm_name}/tracks/align/align_qry_trim-{trim}.bed'),
        asfile=temp('temp/{asm_name}/tracks/align/align_qry_trim-{trim}.as')
    run:

        # Get track description
        if wildcards.trim == 'none':
            track_desc_short = f'PavAlignNone'
            track_description = f'PAV Align (Trim NONE)'

        elif wildcards.trim == 'qry':
            track_desc_short = f'PavAlignQry'
            track_description = f'PAV Align (Trim QRY)'

        elif wildcards.trim == 'qryref':
            track_desc_short = f'PavAlignQryref'
            track_description = f'PAV Align (Trim QRY/REF)'

        else:
            raise RuntimeError('Unknown trim wildcard: '.format(wildcards.trim))

        # Read field table
        df_as = pd.read_csv(
            os.path.join(PIPELINE_DIR, 'files/tracks/alignment_track_fields.tsv'),
            sep='\t'
        ).set_index('FIELD')

        # Read alignments
        df = pd.concat(
            [pd.read_csv(file_name, sep='\t', dtype={'#CHROM': str, 'QRY_ID': str}) for file_name in input.bed],
            axis=0
        )

        del df['CIGAR']

        # Rename SCORE to PAV_SCORE
        df.columns = [{'SCORE': 'PAV_SCORE'}.get(col, col) for col in df.columns]

        # Sort
        df.sort_values(['#CHROM', 'POS', 'END', 'QRY_ID'], inplace=True)

        # Add BED fields
        df['POS_THICK'] = df['POS']
        df['END_THICK'] = df['END']
        df['ID'] = df.apply(lambda row: '{QRY_ID} - {INDEX} ({HAP}-{QRY_ORDER})'.format(**row), axis=1)
        df['SCORE'] = 1000
        df['STRAND'] = df['REV'].apply(lambda val: '-' if val else '+')

        # Set Color
        hap_list = pavlib.pipeline.get_hap_list(wildcards.asm_name, ASM_TABLE)
        colormap_index = np.linspace(0, 0.9999, len(hap_list))
        colormap = mpl.colormaps[ALIGN_COLORMAP]

        hap_color = {  # Color to RGB string (e.g. "(0.267004, 0.004874, 0.329415, 1.0)" from colormap to "68,1,84")
            hap_list[i]: ','.join([str(int(col * 255)) for col in mpl.colors.to_rgb(colormap(colormap_index[i]))])
                for i in range(len(hap_list))
        }

        # hap_color = {
        #     hap_list[i]: colormap(colormap_index[i]) for i in range(len(hap_list))
        # }

        df['COL'] = df['HAP'].apply(lambda val: hap_color[val])

        # Sort columns
        head_cols = ['#CHROM', 'POS', 'END', 'ID', 'SCORE', 'STRAND', 'POS_THICK', 'END_THICK', 'COL']
        tail_cols = [col for col in df.columns if col not in head_cols]

        df = df[head_cols + tail_cols]

        # Check AS fields
        missing_fields = [col for col in df.columns if col not in df_as.index]

        if missing_fields:
            raise RuntimeError('Missing {} fields in AS definition: {}{}'.format(
                len(missing_fields), ', '.join(missing_fields[:3]), '...' if len(missing_fields) else ''
            ))

        # Write AS file
        with open(output.asfile, 'w') as out_file:

            # Heading
            out_file.write('table Align{}\n"{}"\n(\n'.format(track_desc_short, track_description))

            # Column definitions
            for col in df.columns:
                out_file.write('{TYPE} {NAME}; "{DESC}"\n'.format(**df_as.loc[col]))

            # Closing
            out_file.write(')\n')

        # Write BED
        df.to_csv(output.bed, sep='\t', index=False)

#
# Inversion flagged loci
#

rule tracks_invflag_all:
    input:
        bb=lambda wildcards: pavlib.pipeline.expand_pattern(
            'tracks/{asm_name}/inv_flag/inv_flag_{hap}.bb', ASM_TABLE, config
        )

# BigBed for one variant set.
rule tracks_invflag_bb:
    input:
        bed='temp/{asm_name}/tracks/inv_flag/inv_flag_{hap}.bed',
        asfile='temp/{asm_name}/tracks/inv_flag/inv_flag_{hap}.as',
        fai=REF_FAI
    output:
        bb='tracks/{asm_name}/inv_flag/inv_flag_{hap}.bb'
    shell:
        """bedToBigBed -tab -as={input.asfile} -type=bed9+ {input.bed} {input.fai} {output.bb}"""

# Tracks for one variant set.
rule tracks_invflag_bed:
    input:
        bed=lambda wildcards: 'results/{asm_name}/inv_caller/flagged_regions_{hap}_parts-{part_count}.bed.gz'.format(
            part_count=get_config('inv_sig_part_count', wildcards), **wildcards
        ),
        fai=REF_FAI
    output:
        bed=temp('temp/{asm_name}/tracks/inv_flag/inv_flag_{hap}.bed'),
        asfile=temp('temp/{asm_name}/tracks/inv_flag/inv_flag_{hap}.as')
    run:

        field_table_file_name = os.path.join(PIPELINE_DIR, 'files/tracks/inv_flag_fields.tsv')

        color = {
            True: '0,0,0',
            False: '120,120,120'
        }

        # Read variants
        df = pd.read_csv(input.bed, sep='\t')

        df.sort_values(['#CHROM', 'POS', 'END'], inplace=True)

        # Set color
        df['COL'] = df['PARTITION'].apply(lambda val: color[val >= 0])

        # Read FAI and table columns
        df_fai = svpoplib.ref.get_df_fai(input.fai)

        # Filter columns that have track annotations
        field_set = set(
            pd.read_csv(
                field_table_file_name,
                sep='\t', header=0
            )['FIELD']
        )

        df = df.loc[:, [col for col in df.columns if col in field_set]]

        # Make BigBed
        track_name = 'InvFlagTable'
        track_description = '{asm_name} - {hap}'.format(**wildcards)

        svpoplib.tracks.variant.make_bb_track(df, df_fai, output.bed, output.asfile, track_name, track_description, field_table_file_name)

#
# LG-SV - CPX
#


rule tracks_lgsv_cpx_all:
    input:
        bb=lambda wildcards: pavlib.pipeline.expand_pattern(
            'tracks/{asm_name}/lgsv/lgsv_cpx_{hap}.bb', ASM_TABLE, config
        )

rule tracks_lgsv_cpx_bed:
    input:
        bed_cpx='temp/{asm_name}/lg_sv/sv_cpx_{hap}.bed.gz',
        bed_seg='temp/{asm_name}/lg_sv/segment_cpx_{hap}.bed.gz',
        bed_ref='temp/{asm_name}/lg_sv/reftrace_cpx_{hap}.bed.gz',
        fai=REF_FAI
    output:
        bed=temp('temp/{asm_name}/tracks/lgsv/lgsv_cpx_{hap}.bed'),
        asfile=temp('temp/{asm_name}/tracks/lgsv/lgsv_cpx_{hap}.as')
    run:

        field_table_file_name = os.path.join(PIPELINE_DIR, 'files/tracks/lgsv_track_fields.tsv')

        # Read variants
        df_cpx = pd.read_csv(input.bed_cpx, sep='\t', dtype={'#CHROM': str}, low_memory=False)
        df_cpx.sort_values(['#CHROM', 'POS', 'END'], inplace=True)

        # Read structure tables
        df_rt_all = pd.read_csv(
            input.bed_ref, sep='\t', low_memory=False,
            usecols=['#CHROM', 'POS', 'END', 'ID', 'TYPE', 'DEPTH', 'INDEX', 'FWD_COUNT', 'REV_COUNT']
        )[['#CHROM', 'POS', 'END', 'ID', 'TYPE', 'DEPTH', 'INDEX', 'FWD_COUNT', 'REV_COUNT']]

        df_seg_all = pd.read_csv(input.bed_seg, sep='\t', low_memory=False)

        # Read FAI and table columns
        df_fai = svpoplib.ref.get_df_fai(input.fai)

        # Build track table (mix of record types for each complex variant - DEL, DUP, etc)
        df_track_list = list()

        for index, row in df_cpx.iterrows():

            df_rt = df_rt_all.loc[df_rt_all['ID'] == row['ID']].copy()
            df_seg = df_seg_all.loc[df_seg_all['ID'] == row['ID']]

            # if df_rt.shape[0] == 0:
            #     raise RuntimeError(f'Variant {row["ID"]} not found in the CPX reference trace table')

            if df_seg.shape[0] == 0:
                raise RuntimeError(f'Variant {row["ID"]} not found in the CPX segment table')

            # Add unmapped segments to reference trace
            pos = None
            trace_list = list()

            for i in range(df_seg.shape[0]):
                if not df_seg.iloc[i]['IS_ALIGNED'] and not df_seg.iloc[i]['IS_ANCHOR']:
                    if pos is None:
                        raise RuntimeError(f'Segment at position {i} is not aligned and is not preceded by aligned segments')

                    row_seg = df_seg.iloc[i]

                    trace_list.append(
                        pd.Series(
                            [row['#CHROM'], pos, pos + row_seg['LEN_QRY'], row['ID'], 'UNMAPPED', 0, '.', 0, 0],
                            index=['#CHROM', 'POS', 'END', 'ID', 'TYPE', 'DEPTH', 'INDEX', 'FWD_COUNT', 'REV_COUNT']
                        )
                    )

                else:
                    if df_seg.iloc[i]['#CHROM'] == row['#CHROM']:
                        pos = df_seg.iloc[i]['END']

            if len(trace_list) > 0:
                df_rt = pd.concat([df_rt, pd.concat(trace_list, axis=1).T], axis=0)

            for col in ('QRY_REGION', 'QRY_STRAND', 'SEG_N', 'STRUCT_REF', 'STRUCT_QRY', 'VAR_SCORE', 'ANCHOR_SCORE_MIN', 'ANCHOR_SCORE_MAX'):
                df_rt[col] = row[col]

            df_rt.reset_index(drop=True, inplace=True)

            df_rt['ID'] = df_rt.apply(lambda row_rt: row_rt['ID'] + f' ({row_rt.name + 1} / {df_rt.shape[0]})', axis=1)

            df_track_list.append(df_rt)

        if len(df_track_list) > 0:
            df = pd.concat(df_track_list, axis=0).reset_index(drop=True)
        else:
            df = pd.DataFrame([], columns=['#CHROM', 'POS', 'END', 'ID', 'TYPE', 'DEPTH', 'INDEX', 'FWD_COUNT', 'REV_COUNT'])

        # Truncate records that extend off the ends of chromosomes
        if df.shape[0] > 0:
            df_err = df.loc[df.apply(lambda row: row['POS'] < 0, axis=1)]

            if np.any(df_err):
                raise RuntimeError(f'Found {df_err.shape[0]} records with negative POS')

            df['END'] = df.apply(lambda row: np.min([row['END'], df_fai.loc[row['#CHROM']]]), axis=1)

        # Add bed track fields
        df['SCORE'] = 1000
        df['STRAND'] = '.'
        df['POS_THICK'] = df['POS']
        df['END_THICK'] = df['END']

        if df.shape[0] > 0:
            df['COL'] = df.apply(lambda row: CPX_COLOR_DICT[row['TYPE']], axis=1)
        else:
            df['COL'] = np.nan

        # Arrange columns
        head_cols = ['#CHROM', 'POS', 'END', 'ID', 'SCORE', 'STRAND', 'POS_THICK', 'END_THICK', 'COL']
        tail_cols = [col for col in df.columns if col not in head_cols]

        df = df.loc[:, head_cols + tail_cols]

        df.sort_values(['#CHROM', 'POS', 'ID'], inplace=True)

        df['INDEX'] = df['INDEX'].fillna('.')


        ### Define AS columns (AutoSQL, needed to make a BigBed) ###

        # noinspection PyTypeChecker
        df_as = pd.read_csv(
            field_table_file_name,
            sep='\t', header=0,
            dtype={'DEFAULT': object},
            na_values=[''], keep_default_na=False
        )

        df_as.set_index('FIELD', inplace=True, drop=False)

        missing_list = [col for col in tail_cols if col not in set(df_as['FIELD'])]

        if missing_list:
            raise RuntimeError('Missing AS definitions for columns: {}'.format(', '.join(missing_list)))

        # Reformat columns
        for col in df.columns:
            if col == '#CHROM':
                if np.any(pd.isnull(df[col])):
                    raise RuntimeError(f'Error formatting {col}: Found null values in this column (not allowed)')

                continue

            if 'DEFAULT' in df_as.columns and not pd.isnull(df_as.loc[col, 'DEFAULT']):
                default_val = df_as.loc[col, 'DEFAULT']
            else:
                default_val = '.'

            format_type = svpoplib.tracks.variant.TYPE_DICT.get(df_as.loc[col, 'TYPE'], str)

            try:
                df[col] = svpoplib.tracks.variant.format_column(df[col], svpoplib.tracks.variant.TYPE_DICT.get(df_as.loc[col, 'TYPE'], str), default_val=default_val)
            except Exception as ex:
                raise RuntimeError('Error formatting {} as {}: {}'.format(col, df_as.loc[col, 'TYPE'], ex))

        # Write
        track_name = 'LGSVVariantTableCPX'
        track_description = '{asm_name} - LG-SV CPX - {hap}'.format(**wildcards)

        with open(output.asfile, 'w') as out_file:
            # Heading
            out_file.write('table {}\n"{}"\n(\n'.format(track_name, track_description))

            # Column definitions
            for col in head_cols + tail_cols:
                out_file.write('{TYPE} {NAME}; "{DESC}"\n'.format(**df_as.loc[col]))

            # Closing
            out_file.write(')\n')

        df.to_csv(output.bed, sep='\t', na_rep='.', index=False)


#
# LG-SV - INS/DEL/INV
#


rule tracks_lgsv_insdelinv_all:
    input:
        bb=lambda wildcards: pavlib.pipeline.expand_pattern(
            'tracks/{asm_name}/lgsv/lgsv_insdelinv_{hap}.bb', ASM_TABLE, config
        )

rule tracks_lgsv_insdelinv_bed:
    input:
        bed_ins='temp/{asm_name}/lg_sv/svindel_ins_{hap}.bed.gz',
        bed_del='temp/{asm_name}/lg_sv/svindel_del_{hap}.bed.gz',
        bed_inv='temp/{asm_name}/lg_sv/sv_inv_{hap}.bed.gz',
        fai=REF_FAI
    output:
        bed=temp('temp/{asm_name}/tracks/lgsv/lgsv_insdelinv_{hap}.bed'),
        asfile=temp('temp/{asm_name}/tracks/lgsv/lgsv_insdelinv_{hap}.as')
    run:

        field_table_file_name = os.path.join(PIPELINE_DIR, 'files/tracks/lgsv_track_fields.tsv')

        # Read variants
        df = pd.concat(
            [
                pd.read_csv(filename, sep='\t', dtype={'#CHROM': str}, low_memory=False)
                    for filename in [input.bed_ins, input.bed_del, input.bed_inv]
            ]
        ).sort_values(['#CHROM', 'POS', 'END'])

        df = df.loc[df['END'] > df['POS']]  # TEMP - Caller needs to be fixed

        # Read FAI and table columns
        df_fai = svpoplib.ref.get_df_fai(input.fai)

        # Filter columns that have track annotations
        field_set = set(
            pd.read_csv(
                field_table_file_name,
                sep='\t', header=0
            )['FIELD']
        )

        df = df.loc[:, [col for col in df.columns if col in field_set]]

        # Make BigBed
        track_name = 'LGSVVariantTableSV'
        track_description = '{asm_name} - LG-SV INS/DEL/INV - {hap}'.format(**wildcards)

        svpoplib.tracks.variant.make_bb_track(df, df_fai, output.bed, output.asfile, track_name, track_description, field_table_file_name)
