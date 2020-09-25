"""
Input functions for Snakemake rules.
"""

### Align ###

def align_input_fasta(wildcards):
    """
    Get input FASTA file from the assembly table given an assembly name and haplotype.
    """

    if wildcards.asm_name in ASM_TABLE.index:
        # Get table entry

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

    else:
        # Get from pattern
        fa_pattern = config.get('asm_pattern', None)

        if fa_pattern is None:
            raise RuntimeError('Cannot get input fasta for {asm_name}: No assembly table entry and no "fa_pattern" in config'.format(**wildcards))

        if '{asm_name}' not in fa_pattern:
            raise RuntimeError('"Cannot get input fasta for {asm_name}: "asm_pattern" config entry is missing wildcard {{asm_name}}'.format(**wildcards))

        # May have "hap" or "parent" in path (not both). If hap, expect "h1", "h2", etc. If parent, expect "mat" and "pat" (hifiasm-trio)
        if '{hap}' in fa_pattern:
            if '{parent}' in fa_pattern:
                raise RuntimeError('"Cannot get input fasta for {asm_name}: "asm_pattern" contains both {{hap}} and {{parent}} wildcards'.format(**wildcards))

            parent = None

        elif '{parent}' in fa_pattern:

            if wildcards.hap == 'h1':
                parent = 'mat'
            elif wildcards.hap == 'h2':
                parent = 'pat'
            else:
                raise RuntimeError('"Cannot get input fasta for {asm_name} with "{{parent}}" pattern: Haplotype must be "h1" or "h2": Found {hap}'.format(**wildcards))

        else:
            raise RuntimeError('"Cannot get input fasta for {asm_name}: "asm_pattern" config entry is missing wildcard {{hap}}'.format(**wildcards))

        if '{sample}' in fa_pattern:
            re_match = re.match('^([^_]+)_.*', wildcards.asm_name)

            if re_match is None:
                raise RuntimeError('"Cannot get input fasta for {asm_name}: "asm_pattern" contains wildcard {{sample}}, but sample cannot be extracted from the assembly name (expected to be at the start of the assembly name and separated by an underscore)'.format(**wildcards))

            sample = re_match[1]
        else:
            sample = None

        # Create file
        return fa_pattern.format(
            sample=sample,
            asm_name=wildcards.asm_name,
            hap=wildcards.hap,
            parent=parent
        )
