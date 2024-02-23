"""
Prepare UCSC tracks for data.
"""

import matplotlib as mpl
import numpy as np
import os
import pandas as pd

import pavlib
import svpoplib

global PIPELINE_DIR
global REF_FAI

global temp
global ASM_TABLE



#
# Definitions
#

VARTYPE_TO_SVTYPE_TUPLE = {
    'snv': ('snv'),
    'sv': ('ins', 'del'),
    'indel': ('ins', 'del')
}

ALIGN_COLORMAP = 'viridis'

# ALIGN_COLOR = {
#     'h1': '160,00,144', # Pink
#     'h2': '64,00,160'   # Purple
# }

VAR_BED_PREMERGE_PATTERN = 'results/{{asm_name}}/bed/pre_merge/{{hap}}/{vartype}_{svtype}.bed.gz'

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

    # Return list of input files
    return [
        VAR_BED_PREMERGE_PATTERN.format(vartype=input_vartype, svtype=svtype)
            for svtype in input_svtype
    ]


#
# Variant calls
#

# BigBed for one variant set.
rule tracks_hap_call_bb:
    input:
        bed='temp/{asm_name}/tracks/bed_pre_merge/{vartype}_{svtype}_{hap}.bed',
        asfile='temp/{asm_name}/tracks/bed_pre_merge/{vartype}_{svtype}_{hap}.as',
        fai=REF_FAI
    output:
        bb='tracks/{asm_name}/variant/pre_merge/{vartype}_{svtype}_{hap}.bb'
    shell:
        """bedToBigBed -tab -as={input.asfile} -type=bed9+ {input.bed} {input.fai} {output.bb}"""

# Tracks for one variant set.
rule tracks_hap_call:
    input:
        bed=_track_get_input_bed,
        fai=REF_FAI
    output:
        bed=temp('temp/{asm_name}/tracks/bed_pre_merge/{vartype}_{svtype}_{hap}.bed'),
        asfile=temp('temp/{asm_name}/tracks/bed_pre_merge/{vartype}_{svtype}_{hap}.as')
    run:

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

# Alignment track BED to BigBed.
rule tracks_align_bb:
    input:
        bed='temp/{asm_name}/tracks/align/tig_align_trim-{trim}.bed',
        asfile='temp/{asm_name}/tracks/align/tig_align_trim-{trim}.as'
    output:
        bb='tracks/{asm_name}/align/tig_align_trim-{trim}.bb'
    shell:
        """bedToBigBed -tab -as={input.asfile} -type=bed9+ {input.bed} {REF_FAI} {output.bb}"""

# Alignment tracks.
rule tracks_align:
    input:
        bed_h1='results/{asm_name}/align/trim-{trim}/aligned_tig_h1.bed.gz',
        bed_h2='results/{asm_name}/align/trim-{trim}/aligned_tig_h2.bed.gz'
    output:
        bed=temp('temp/{asm_name}/tracks/align/tig_align_trim-{trim}.bed'),
        asfile=temp('temp/{asm_name}/tracks/align/tig_align_trim-{trim}.as')
    wildcard_constraints:
        trim='none|tig|tigref'
    run:

        # Get track description
        if wildcards.trim == 'none':
            track_desc_short = 'PreTrim'
            track_description = 'Pre-trimmed alignments'

        elif wildcards.trim == 'tig':
            track_desc_short = 'TrimTig'
            track_description = 'Alignments (contig trimmed)'

        elif wildcards.trim == 'tigref':
            track_desc_short = 'TrimTigRef'
            track_description = 'Alignments (contig & ref trimmed)'

        else:
            raise RuntimeError('Unknown trim wildcard: '.format(wildcards.trim))

        # Read field table
        df_as = pd.read_csv(
            os.path.join(PIPELINE_DIR, 'files/tracks/alignment_track_fields.tsv'),
            sep='\t'
        ).set_index('FIELD')

        # Read alignments
        df = pd.concat(
            [pd.read_csv(file_name, sep='\t', dtype={'#CHROM': str, 'QRY_ID': str}) for file_name in [input.bed_h1, input.bed_h2]],
            axis=0
        )

        del df['CIGAR']

        # Set contig order
        df.sort_values(['QRY_ID', 'QRY_POS', 'QRY_END'], inplace=True)
        df['QRY_ORDER'] = -1

        last_qry = None
        qry_order = 0

        for index, row in df.iterrows():
            if row['QRY_ID'] != last_qry:
                last_qry = row['QRY_ID']
                qry_order = 0

            df.loc[index, 'QRY_ORDER'] = qry_order
            qry_order += 1

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
