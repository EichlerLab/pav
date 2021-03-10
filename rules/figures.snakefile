"""
Generate figures from variant discovery.

The rules in this snakefile may not generate variant calls. The pipeline should be run first, then variants can be
plotted with rules in this file. These may be useful for reporting results or troubleshooting.
"""

#
# Inversions
#

# figures_inv_dot_density
#
# Make dot and density plot for an inversion call. May plot with inverted repeat whiskers (whiskers=whisk) or without
# (whiskers=nowhisk). These whiskers show k-mers unique to each inverted repeat flanking the inversion and whether they
# are in reference orientation (pointing up) or inverted orientation (pointing down).
rule figures_inv_dot_density:
    input:
        tig_fa='temp/{asm_name}/align/contigs_{hap}.fa.gz',
        bed_inv='results/{asm_name}/{flag_source}/sv_inv_{hap}.bed.gz',
        bed_density='results/{asm_name}/{flag_source}/density_table/density_{inv_id}_{hap}.tsv.gz'
    output:
        fig_dot='figures/inv/{asm_name}/{hap}/{flag_source}/{inv_id}_dot.{ext}',
        fig_den='figures/inv/{asm_name}/{hap}/{flag_source}/{inv_id}_density.{ext}'
    wildcard_constraints:
        ext='pdf|png'
    run:

        # Read inversion calls and get record
        df_inv = pd.read_csv(input.bed_inv, sep='\t', index_col='ID')

        if wildcards.inv_id not in df_inv.index:
            raise RuntimeError(
                f'Inversion not found in {wildcards.inv_id} callset '
                f'(assembly={wildcards.asm_name}, hap={wildcards.hap})'
            )

        inv_row = df_inv.loc[wildcards.inv_id]

        df_density = pd.read_csv(input.bed_density, sep='\t')

        # Get inversion call
        inv_call = pavlib.inv.get_inv_from_record(inv_row, df_density)

        # Make plots
        fig_dot = pavlib.plot.dotplot_inv_call(
            inv_call, REF_FA, seq_tig=inv_row['SEQ']
        )

        fig_density = pavlib.plot.kmer_density_plot(
            inv_call, hap=wildcards.hap, flank_whiskers=True
        )

        # Write plots
        fig_dot.savefig(output.fig_dot, bbox_inches='tight')
        fig_density.savefig(output.fig_den, bbox_inches='tight')
