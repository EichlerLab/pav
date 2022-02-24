"""
Input functions for Snakemake rules.
"""

### Align ###

def align_input_list(wildcards):
    """
    Get a list of input files. If there are no input files (empty haplotype), an empty list is returned.
    """

    # Determine if wildcards.hap is present and properly formatted
    if 'hap' in wildcards.keys():

        if not re.match('^h\d+$', wildcards.hap):
            raise RuntimeError('Haplotype wildcard "hap" is malformed: Expected "h" followed by an integer: {hap}'.format(hap=wildcards.hap))

        wildcards_has_hap = True
    else:
        wildcards_has_hap = False

    if 'asm_name' not in wildcards.keys():
        raise RuntimeError('Missing wildcard pattern "asm_name"')

    # Get filename pattern form either the assembly table or config entry (assembly table takes precedent)
    if wildcards.asm_name.strip() in ASM_TABLE.index:
        # Get from assembly table

        record_source = 'table'

        if not wildcards_has_hap:
            raise RuntimeError('Cannot get assembly table record for {asm_name}: Wildcards has no "hap" element'.format(asm_name=wildcards.asm_name))

        # Get table entry
        asm_table_entry = ASM_TABLE.loc[wildcards.asm_name]

        # Get haplotype
        haplotype_col = re.sub('^h', 'HAP', wildcards.hap)

        if haplotype_col not in asm_table_entry:
            raise RuntimeError(
                'No haplotype column "{}" in assembly table for haplotype "{}"'.format(
                    haplotype_col, wildcards.hap
                )
            )

        filename_pattern = asm_table_entry[haplotype_col]

        if not pd.isnull(filename_pattern):
            filename_pattern = filename_pattern.strip()

            if not filename_pattern:
                filename_pattern = None
        else:
            filename_pattern = None

    else:
        # Get from config pattern

        record_source = 'config'

        filename_pattern = config.get('asm_pattern', None)

        if filename_pattern is None:
            raise RuntimeError('Cannot get input file for {asm_name}: No assembly table entry and no "asm_pattern" in config'.format(**wildcards))

    # Return empty list if there are no paths
    if filename_pattern is None:
        return list()

    # Separate semi-colon-delimited paths and fit wildcards
    path_list = list()

    for file_name in filename_pattern.split(';'):

        file_name = file_name.strip()

        if not file_name:
            continue

        # Check for "parent" wildcard (hifiasm-trio)
        if '{parent}' in file_name:

            if wildcards.hap == 'h1':
                parent = 'mat'
            elif wildcards.hap == 'h2':
                parent = 'pat'
            else:
                raise RuntimeError('"Cannot get input file for {asm_name} with "{{parent}}" pattern: Haplotype must be "h1" or "h2": Found {hap}: {file_name}'.format(asm_name=wildcards.asm_name, file_name=file_name, hap=wildcards.hap))

        else:
            parent = None

        # Check for required fields if the record comes from a config entry (must have asm_name and hap or parent)
        if record_source == 'config':

            if '{asm_name}' not in file_name:
                raise RuntimeError('Cannot get input file for {asm_name}: Path from config "asm_name" must have {{asm_name}} wildcard: {file_name}'.format(asm_name=wildcards.asm_name, file_name=file_name))

            if '{hap}' not in file_name and parent is None:
                raise RuntimeError('Cannot get input file for {asm_name}: Path from config "asm_pattern" must have {{hap}} or {{parent}} wildcard: {file_name}'.format(asm_name=wildcards.asm_name, file_name=file_name))

        if '{sample}' in file_name:
            delim = config.get('sample_delimiter', '_')

            if delim == '.':
                delim = '\.'

            elif delim not in {'_', '-', '+', '#'}:
                raise RuntimeError('Sample delimiter in config ("sample_delimiter") must be one of "_", ".", "-", "+", "#": {}'.format(delim))

            re_match = re.match('^([^{delim}]+){delim}.*'.format(delim=delim), wildcards.asm_name)

            if re_match is None:
                raise RuntimeError('"Cannot get input fasta for {asm_name}: "asm_pattern" contains wildcard {{sample}}, but sample cannot be extracted from the assembly name (expected to be at the start of the assembly name and separated by an underscore): {file_name}'.format(asm_name=wildcards.asm_name, file_name=file_name))

            sample = re_match[1]
        else:
            sample = None

        # Create filename
        file_name_parsed = file_name.format(
            sample=sample,
            asm_name=wildcards.asm_name,
            hap=wildcards.hap,
            parent=parent
        )

        if not os.path.isfile(file_name_parsed):
            # Try substituting "hap1" for "h1" and "hap2" for "h2" if h1/h2 is not found.
            # Return this version of the file if it exists, let the pipeline fail on the h1/h2 filename if neither exists.

            if re.match('^h\d+$', wildcards.hap):
                alt_hap = re.sub('^h', 'hap', wildcards.hap)

                file_name_parsed_alt = file_name.format(
                    sample=sample,
                    asm_name=wildcards.asm_name,
                    hap=alt_hap,
                    parent=parent
                )

                if os.path.isfile(file_name_parsed_alt):
                    file_name_parsed = file_name_parsed_alt

        if file_name_parsed is not None:
            path_list.append(file_name_parsed)

        # Return paths
        return path_list
