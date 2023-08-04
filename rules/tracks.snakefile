"""
Prepare UCSC tracks for data.
"""

#
# Definitions
#

VARTYPE_TO_SVTYPE_TUPLE = {
    'snv': ('snv'),
    'sv': ('ins', 'del'),
    'indel': ('ins', 'del')
}

ALIGN_COLOR = {
    'h1': '160,00,144', # Pink
    'h2': '64,00,160'   # Purple
}

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

def _get_align_bed(wildcards):

    if wildcards.align_stage == 'pre-cut':
        return [
            'results/{asm_name}/align/pre-cut/aligned_tig_h1.bed.gz'.format(asm_name=wildcards.asm_name),
            'results/{asm_name}/align/pre-cut/aligned_tig_h2.bed.gz'.format(asm_name=wildcards.asm_name)
        ]

    if wildcards.align_stage == 'post-cut':
        return [
            'results/{asm_name}/align/aligned_tig_h1.bed.gz'.format(asm_name=wildcards.asm_name),
            'results/{asm_name}/align/aligned_tig_h2.bed.gz'.format(asm_name=wildcards.asm_name)
        ]

    raise RuntimeError('Unknown align_stage wildcard: '.format(wildcards.align_stage))

#
# Variant calls
#

# tracks_hap_call_bb
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

# tracks_hap_call
#
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
            [pd.read_csv(file_name, sep='\t') for file_name in input.bed],
            axis=0
        ).reset_index(drop=True)

        if wildcards.vartype == 'sv':
            df = df.loc[df['SVLEN'] >= 50].copy()
        elif wildcards.vartype == 'indel':
            df = df.loc[df['SVLEN'] < 50].copy()

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Select table columns
        if 'QUERY_ID' in df.columns:
            del(df['QUERY_ID'])

        if 'QUERY_POS' in df.columns:
            del(df['QUERY_POS'])

        if 'QUERY_END' in df.columns:
            del(df['QUERY_END'])

        if 'QUERY_STRAND' in df.columns:
            del(df['QUERY_STRAND'])

        if 'SEQ' in df.columns:
            del(df['SEQ'])

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

# tracks_align_bb
#
# Alignment track BED to BigBed.
rule tracks_align_bb:
    input:
        bed='temp/{asm_name}/tracks/align/{align_stage}.bed',
        asfile='temp/{asm_name}/tracks/align/{align_stage}.as',
        fai=REF_FAI
    output:
        bb='tracks/{asm_name}/align/{align_stage}.bb'
    shell:
        """bedToBigBed -tab -as={input.asfile} -type=bed9+ {input.bed} {input.fai} {output.bb}"""

# tracks_align
#
# Alignment tracks.
rule tracks_align:
    input:
        bed=_get_align_bed
    output:
        bed=temp('temp/{asm_name}/tracks/align/{align_stage}.bed'),
        asfile=temp('temp/{asm_name}/tracks/align/{align_stage}.as')
    wildcard_constraints:
        align_stage='(pre|post)-cut'
    run:

        # Get track description
        if wildcards.align_stage == 'pre-cut':
            track_desc_short = 'PreTrim'
            track_description = 'Pre-trimmed alignments'

        elif wildcards.align_stage == 'post-cut':
            track_desc_short = 'PostTrim'
            track_description = 'Post-trimmed alignments'

        else:
            raise RuntimeError('Unknown align_stage wildcard: '.format(wildcards.align_stage))

        # Read field table
        df_as = pd.read_csv(
            os.path.join(PIPELINE_DIR, 'files/tracks/alignment_track_fields.tsv'),
            sep='\t'
        ).set_index('FIELD')

        # Read variants
        df = pd.concat(
            [pd.read_csv(file_name, sep='\t') for file_name in input.bed],
            axis=0
        ).sort_values(['#CHROM', 'POS', 'END', 'QUERY_ID'])

        # Add BED fields
        df['POS_THICK'] = df['POS']
        df['END_THICK'] = df['END']
        df['ID'] = df['QUERY_ID']
        df['SCORE'] = 1000
        df['COL'] = df['HAP'].apply(lambda val: ALIGN_COLOR.get(val, '0,0,0'))
        df['STRAND'] = df['REV'].apply(lambda val: '-' if val else '+')

        del(df['CIGAR'])

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
