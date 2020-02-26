"""
Input functions for Snakemake rules.
"""

### Align ###

def align_input_fasta(wildcards):
    """
    Get input FASTA file from the assembly table given an assembly name and haplotype.
    """

    # Get table entry
    if wildcards.asm_name not in ASM_TABLE.index:
        raise RuntimeError('Cannot get input fasta for {asm_name}: No assembly table entry'.format(**wildcards))

    asm_table_entry = ASM_TABLE.loc[wildcards.asm_name]

    # Get haplotype
    if wildcards.hap not in {'h1', 'h2', 'h0'}:
        raise RuntimeError('Unknown haplotype: {hap}'.format(**wildcards))

    haplotype_col = 'HAP{}'.format(wildcards.hap[1:])

    if haplotype_col not in asm_table_entry:
        raise RuntimeError(
            'No haplotype column "{}" in assembly table for haplotype "{}"'.format(
                haplotype_col, wildcards.hap
            )
        )

    return asm_table_entry[haplotype_col]
