"""
Prepare UCSC tracks for data.
"""

VARTYPE_TO_SVTYPE_TUPLE = {
    'snv': ('snv'),
    'sv': ('ins', 'del'),
    'indel': ('ins', 'del')
}

def _track_get_input_bed(wildcards):
    """
    Get one or more input files for tracks. If "svtype" is "all", collect all relevant input files.

    :param wildcards: Wildcards.

    :return: List of input file(s).
    """

    if wildcards.svtype == 'all':

        if wildcards.vartype not in {'sv', 'indel'}:
            raise RuntimeError('Cannot collect "all" for variant type: {} (must be a valid variant type with more than one svtype)'.format(wildcards.vartype))

        return [
            'results/{{asm_name}}/bed/pre_merge/{{vartype}}_{svtype}_{{hap}}.bed.gz'.format(svtype=svtype)
            for svtype in VARTYPE_TO_SVTYPE_TUPLE[wildcards.vartype]
        ]

    return ['results/{asm_name}/bed/pre_merge/{vartype}_{svtype}_{hap}.bed.gz']

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

        field_table_file_name = os.path.join(PIPELINE_DIR, 'files/tracks/ucsc_track_fields.tsv')

        # Read variants
        df = pd.concat(
            [pd.read_csv(file_name, sep='\t') for file_name in input.bed],
            axis=0
        )

        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # Select table columns
        del(df['QUERY_ID'])
        del(df['QUERY_POS'])
        del(df['QUERY_END'])
        del(df['QUERY_STRAND'])

        if wildcards.vartype not in {'snv', 'indel'}:
            del(df['SEQ'])

        # Read FAI and table columns
        df_fai = analib.ref.get_df_fai(input.fai)

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

        analib.tracks.variant.make_bb_track(df, df_fai, output.bed, output.asfile, track_name, track_description, field_table_file_name)
