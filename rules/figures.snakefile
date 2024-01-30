"""
Generate figures from variant discovery.

The rules in this snakefile may not generate variant calls. The pipeline should be run first, then variants can be
plotted with rules in this file. These may be useful for reporting results or troubleshooting.
"""

#
# Inversions
#


#######################
### Input functions ###
#######################

def fig_input_inv_tuples(asm_name, filter, hap_list):

    tuple_list = list()

    for hap in hap_list:
        for index, row in pd.read_csv(
            f'results/{asm_name}/bed/{filter}/{hap}/sv_inv.bed.gz', sep='\t', usecols=('ID', 'CALL_SOURCE')
        ).iterrows():
            call_source, den_flag = row['CALL_SOURCE'].split('-')

            tuple_list.append((
                row['ID'], hap, call_source, den_flag
            ))

    return tuple_list

def fig_input_inv_all_dot(wildcards):
    """
    Get an iterator for all dotplot output file names for an assembly.
    """

    for inv_id, hap, call_source, den_flag in fig_input_inv_tuples(
        wildcards.asm_name, wildcards.filter,
        pavlib.pipeline.get_hap_list(wildcards.asm_name, ASM_TABLE)
    ):
        yield f'figures/inv/{wildcards.asm_name}/{wildcards.filter}/{inv_id}_{hap}_dot.{wildcards.ext}'

def fig_input_inv_all_den(wildcards):
    """
    Get an iterator for all dotplot output file names for an assembly.
    """

    for inv_id, hap, call_source, den_flag in fig_input_inv_tuples(
        wildcards.asm_name, wildcards.filter,
        pavlib.pipeline.get_hap_list(wildcards.asm_name, ASM_TABLE)
    ):
        if den_flag != 'NODEN':
            yield f'figures/inv/{wildcards.asm_name}/{wildcards.filter}/{inv_id}_{hap}_den.{wildcards.ext}'



###################
### All figures ###
###################

# Make density and dotplot figures for all inversions.
localrules: figures_inv_all

rule figures_inv_all:
    input:
        fig_dot=fig_input_inv_all_dot,
        fig_den=fig_input_inv_all_den
    output:
        flag=touch('temp/flag/fig/{asm_name}_{filter}_{ext}.flag')


###############
### Density ###
###############

# Make dot and density plot for an inversion call. May plot with inverted repeat whiskers (whiskers=whisk) or without
# (whiskers=nowhisk). These whiskers show k-mers unique to each inverted repeat flanking the inversion and whether they
# are in reference orientation (pointing up) or inverted orientation (pointing down).
rule figures_inv_den:  # Replacing with one rule per dot or density
    input:
        tig_fa='temp/{asm_name}/align/contigs_{hap}.fa.gz',
        bed_inv='results/{asm_name}/bed/{filter}/{hap}/sv_inv.bed.gz'
    output:
        fig_den='figures/inv/{asm_name}/{filter}/{inv_id}_{hap}_den.{ext}'
    run:

        # Read inversion calls and get record
        df_inv = pd.read_csv(input.bed_inv, sep='\t', index_col='ID')

        if wildcards.inv_id not in df_inv.index:
            raise RuntimeError(
                f'Inversion not found in {wildcards.inv_id} callset '
                f'(assembly={wildcards.asm_name}, hap={wildcards.hap})'
            )

        inv_row = df_inv.loc[wildcards.inv_id]

        # Get density BED path
        call_source, den_flag = inv_row['CALL_SOURCE'].split('-')

        if den_flag != 'DEN':
            raise RuntimeError(
                'Cannot generate a density plot for an inversion that was not called with density: '
                f'assembly={wildcards.asm_name}, hap={wildcards.hap}, call source={inv_row["CALL_SOURCE"]}'
            )

        if call_source == 'FLAG':
            density_table_file_name = f'results/{wildcards.asm_name}/inv_caller/density_table/density_{wildcards.inv_id}_{wildcards.hap}.tsv.gz'

        elif call_source == 'ALNTRUNC':
            density_table_file_name = f'results/{wildcards.asm_name}/inv_caller/density_table_lg/density_{wildcards.inv_id}_{wildcards.hap}.tsv.gz'

        else:
            raise RuntimeError(
                'Cannot retrieve density table: Unrecognized call source: '
                f'assembly={wildcards.asm_name}, hap={wildcards.hap}, call source={inv_row["CALL_SOURCE"]}'
            )

        if not os.path.isfile(density_table_file_name):
            raise RuntimeError(
                'Missing density table: '
                f'assembly={wildcards.asm_name}, hap={wildcards.hap}, call source={inv_row["CALL_SOURCE"]}: '
                + density_table_file_name
            )

        df_density = pd.read_csv(density_table_file_name, sep='\t')

        # Get inversion call
        inv_call = pavlib.inv.get_inv_from_record(inv_row, df_density)

        # Get contig sequence
        with pysam.FastaFile(input.tig_fa) as fa_file:
            seq_tig = fa_file.fetch(
                inv_call.region_tig_discovery.chrom,
                inv_call.region_tig_discovery.pos,
                inv_call.region_tig_discovery.end
            )

        # Make plots
        fig_density = pavlib.plot.kmer_density_plot(
            inv_call, hap=wildcards.hap, flank_whiskers=False
        )

        # Write plots
        fig_density.savefig(output.fig_den, bbox_inches='tight')


# Inversion dotplot.
rule figures_inv_dot:
    input:
        bed_inv='results/{asm_name}/bed/{filter}/{hap}/sv_inv.bed.gz',
        tig_fa='temp/{asm_name}/align/contigs_{hap}.fa.gz',
        tig_fai='temp/{asm_name}/align/contigs_{hap}.fa.gz.fai'
    output:
        fig_dot='figures/inv/{asm_name}/{filter}/{inv_id}_{hap}_dot.{ext}'
    run:

        pad_region = 0.20  # Pad each end by this proportion when extracting regions.

        # Read
        df_inv = pd.read_csv(input.bed_inv, sep='\t', index_col='ID')

        if wildcards.inv_id not in df_inv.index:
            raise RuntimeError(
                f'Inversion not found: {wildcards.inv_id} '
                f'(assembly={wildcards.asm_name}, hap={wildcards.hap})'
            )

        inv_row = df_inv.loc[wildcards.inv_id]

        # Get FAIs
        ref_fai = svpoplib.ref.get_df_fai(REF_FAI)
        tig_fai = svpoplib.ref.get_df_fai(input.tig_fai)

        # Get inversion call object
        inv_call = pavlib.inv.get_inv_from_record(inv_row, None)

        # Get contig and reference region
        tig_outer_len = inv_call.region_tig_outer.end - inv_call.region_tig_outer.pos
        ref_outer_len = inv_call.region_ref_outer.end - inv_call.region_ref_outer.pos

        tig_pad = int(tig_outer_len * pad_region)
        ref_pad = int(ref_outer_len * pad_region)

        # Expand contig region
        tig_pos = inv_call.region_tig_outer.pos
        tig_end = inv_call.region_tig_outer.end

        if inv_call.region_tig_discovery is not None:
            if inv_call.region_tig_discovery.pos < tig_pos:
                tig_pos = inv_call.region_tig_discovery.pos

            if inv_call.region_tig_discovery.end > tig_end:
                tig_end = inv_call.region_tig_discovery.end

        tig_pos = np.max([
            0,
            tig_pos - tig_pad,
        ])

        tig_end = np.min([
            tig_end + tig_pad,
            tig_fai[inv_call.region_tig_outer.chrom]
        ])

        # Expand reference region
        ref_pos = inv_call.region_ref_outer.pos
        ref_end = inv_call.region_ref_outer.end

#         if inv_call.region_ref_discovery is not None:
#             if inv_call.region_ref_discovery.pos < ref_pos:
#                 ref_pos = inv_call.region_ref_discovery.pos
#
#             if inv_call.region_ref_discovery.end > ref_end:
#                 ref_end = inv_call.region_ref_discovery.end

        ref_pos = np.max([
            0,
            ref_pos - ref_pad,
        ])

        ref_end = np.min([
            ref_end + ref_pad,
            ref_fai[inv_call.region_ref_outer.chrom]
        ])
        
        # Set regions
        region_tig = pavlib.seq.Region(
            inv_call.region_tig_outer.chrom, tig_pos, tig_end
        )

        region_ref = pavlib.seq.Region(
            inv_call.region_ref_outer.chrom, ref_pos, ref_end
        )

        # Get tig region
        with pysam.FastaFile(input.tig_fa) as fa_file:
            seq_tig = fa_file.fetch(
                region_tig.chrom,
                region_tig.pos,
                region_tig.end
            )

        # Create figure
        fig_dot = pavlib.plot.dotplot_inv_call(
            inv_call, REF_FA, seq_tig=seq_tig,
            region_tig=region_tig, region_ref=region_ref
        )

        # Write
        fig_dot.savefig(output.fig_dot, bbox_inches='tight')
