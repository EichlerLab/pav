# VCF writing procedures

import gzip
import numpy as np
import os
import pandas as pd
import pysam

import Bio.bgzf

import pavlib
import svpoplib


def write_merged_vcf(asm_name, input_dict, output_filename, reference_filename, ref_tsv, symbolic_alt=('sv_inv'), symbolic_seq=None):
    """
    Write a merged VCF file.

    * input_dict
      * key: Tuple of...
        * [0]: varsvtype
        * [1]: "pass" or "fail"
      * Value: Tuple of...
        * [0]: Variant BED file name.
        * [1]: Variant FASTA file name (None if no variant sequences are not used in the VCF).

    :param asm_name: Assembly name. Becomes the name in the VCF column.
    :param input_dict: Input dictionary (see above).
    :param output_filename: Output filename. Will be written as a bgzipped VCF file regardless of file name.
    :param reference_filename: Name of reference file.
    :param ref_tsv: Filename for a TSV of information for each reference sequence (name, length, MD5) or a
        `pd.DataFrame` object of the table already read into memory.
    :param symbolic_alt: Tuple of vartype_svtype values that should get symbolic ALTs in the VCF.
    :param symbolic_seq: Tuple of vartype-svtype values that should get INFO/SEQ if they are also symbolic (no effect
        if the vartype-svtype is not also in symbolic_alt).
    """

    # Check parameters
    if symbolic_alt is not None:
        if isinstance(symbolic_alt, str):
            symbolic_alt = {symbolic_alt}
        else:
            symbolic_alt = set(symbolic_alt)
    else:
        symbolic_alt = {}

    if symbolic_seq is not None:
        if isinstance(symbolic_seq, str):
            symbolic_seq = {symbolic_seq}
        else:
            symbolic_seq = set(symbolic_seq)
    else:
        symbolic_seq = {}

    for key, input_tuple in input_dict.items():
        if not isinstance(key, tuple) or len(key) != 2 or key[0] is None or key[1] is None:
            raise RuntimeError(f'input_dict key is not a tuple (varsvtype, filter): "{key}"')

        if not isinstance(input_tuple, tuple) or len(input_tuple) != 2 or key[0] is None or not isinstance(key[0], str) or not (key[1] is None or isinstance(key[1], str)):
            raise RuntimeError(f'input_dict value for key "{key}" is not a tuple (BED filename, filter): "{input_tuple}"')

        if key[0] in symbolic_alt and key[0] in symbolic_seq and key[1] is None:
            raise RuntimeError(f'input_dict value for key "{key}" is a symbolic ALT with INFO/SEQ (symbolic_alt and symbolic_seq parameters), but no input FASTA file is given')

    # Check assembly name
    if asm_name in {'#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'}:
        raise RuntimeError(f'Assembly name conflicts with a VCF header column name: {asm_name}')

    # Set of known filters
    known_filter_set = set(pavlib.call.FILTER_REASON.keys())

    # Read reference TSV
    if isinstance(ref_tsv, str):
        df_ref = pd.read_csv(ref_tsv, sep='\t')
    elif isinstance(ref_tsv, pd.DataFrame):
        df_ref = ref_tsv
    else:
        raise RuntimeError(f'Unknown type for reference TSV parameter (ref_tsv): {type(ref_tsv)}')

    # Process variant types
    df_list = list()

    any_info_seq = False       # True if the SEQ column is present for any variant
    symbolic_alt_set = set()  # Set of SVTYPEs (uppercase) where symbolic ALTs were written (saved for the header)

    for (varsvtype, filter), input_tuple in input_dict.items():
        vartype, svtype = varsvtype.split('_')

        # Symbolic ALT and INFO/SEQ for symbolic ALTs
        if varsvtype in symbolic_alt:
            is_symbolic_alt = True
            symbolic_alt_set.add(svtype.upper())

            if varsvtype in symbolic_seq:
                is_info_seq = True
                any_info_seq = True

                if input_tuple[1] is None:
                    raise RuntimeError(f'Symbolic ALTs with INFO/SEQ was requested for {vartype}-{svtype}, but no sequence FA file was found.')

            else:
                is_info_seq = False

        else:
            is_symbolic_alt = False
            is_info_seq = False

        if svtype == 'inv' and not is_symbolic_alt:
            raise RuntimeError('INV found without symbolic ALTs set')

        # Read variants
        df = pd.read_csv(input_tuple[0], sep='\t', dtype={'#CHROM': str})

        df.set_index('ID', inplace=True, drop=False)
        df.index.name = 'INDEX'

        if df.shape[0] == 0:
            continue

        # Set FILTER
        if 'FILTER' not in df.columns:
            df['FILTER'] = 'PASS'

        df['FILTER'] = df['FILTER'].apply(lambda val: val.strip() if not pd.isnull(val) and val.strip() != '' else 'PASS')

        uk_filter = {
            val for index, row in df.iterrows() for val in set(row['FILTER'].split(','))
        } - known_filter_set

        if uk_filter:
            uk_list = ', '.join(sorted(uk_filter)[:3]) + ('...' if len(uk_filter) > 3 else '')
            raise RuntimeError(f'Unknown filter(s) found in variant file (var-svtype={varsvtype}, filter={filter}: {uk_list}')

        # Set VARTYPE
        if vartype != 'svindel':
            df['VARTYPE'] = vartype.upper()
        else:
            df['VARTYPE'] = df['SVLEN'].apply(lambda svlen: 'SV' if svlen >= 50 else 'INDEL')

        # Read sequence from FASTA
        if input_tuple[1] is not None:

            # Get Sequence from FASTA and assign to SEQ column (not for SNVs)
            df['SEQ'] = ''

            with gzip.open(input_tuple[1], 'rt') as fa_in:
                for record in Bio.SeqIO.parse(fa_in, 'fasta'):
                    if record.name in df.index:
                        df.loc[record.name, 'SEQ'] = str(record.seq)

            # Check for missing sequences
            df_null_seq = df.loc[pd.isnull(df['SEQ'])]

            if df_null_seq.shape[0] > 0:
                raise RuntimeError(
                    'Missing FASTA sequence for {} variants (vartype={}, svtype={}): {}{}'.format(
                        df_null_seq.shape[0],
                        vartype,
                        svtype,
                        ', '.join([str(val) for val in df_null_seq.iloc[:3]['ID']]),
                        '...' if df_null_seq.shape[0] > 3 else ''
                    )
                )

        # Reformat fields for INFO
        for col in ('HAP', 'HAP_VARIANTS', 'CALL_SOURCE', 'QRY_REGION', 'QRY_STRAND', 'COV_MEAN', 'COV_PROP', 'RGN_REF_INNER', 'RGN_QRY_INNER'):
            if col in df.columns:
                df[col] = df[col].astype(str).apply(lambda val: val.replace(';', ',') if not pd.isnull(val) else val)

        if svtype == 'del':
            df['SVLEN'] = - np.abs(df['SVLEN'])

        # INFO: Base
        df['INFO'] = df.apply(lambda row: 'ID={ID};SVTYPE={SVTYPE}'.format(**row), axis=1)

        # INFO: Add SV/INDEL annotations
        if vartype != 'snv':
            df['INFO'] = df.apply(lambda row: row['INFO'] + ';SVLEN={SVLEN}'.format(**row), axis=1)

        # INFO: Universal INFO columns (all variant types)
        df['INFO'] = df.apply(lambda row: row['INFO'] + ';HAP={HAP};HAP_VARIANTS={HAP_VARIANTS};COV_MEAN={COV_MEAN};COV_PROP={COV_PROP};QRY_REGION={QRY_REGION};QRY_STRAND={QRY_STRAND};CALL_SOURCE={CALL_SOURCE}'.format(**row), axis=1)

        # INFO: Add INV
        if svtype == 'inv':
            df['INFO'] = df.apply(lambda row: row['INFO'] + ';INNER_REF={RGN_REF_INNER};INNER_TIG={RGN_QRY_INNER}'.format(**row), axis=1)

        # INFO: Add COMPOUND
        if 'COMPOUND' in df.columns:
            df['INFO'] = df.apply(lambda row:
                row['INFO'] + (';COMPOUND={COMPOUND}'.format(**row) if not pd.isnull(row['COMPOUND']) else ''),
                axis=1
            )

        # INFO: Add breakpoint homology
        # if svtype in {'ins', 'del'}:
        #     df['INFO'] = df.apply(lambda row: row['INFO'] + ';HOM_REF={HOM_REF};HOM_TIG={HOM_TIG}'.format(**row), axis=1)

        # REF
        if 'REF' not in df.columns:
            if not is_symbolic_alt:
                df['REF'] = list(svpoplib.vcf.ref_base(df, reference_filename))
            else:

                def get_ref_base_by_row(row):
                    pos = max([0, row['POS'] - 1])
                    return ref_file.fetch(row['#CHROM'], pos, pos + 1).upper()

                with pysam.FastaFile(reference_filename) as ref_file:
                    df['REF'] = df.apply(get_ref_base_by_row, axis=1)

        # ALT
        if vartype != 'snv':
            if is_symbolic_alt:

                df['ALT'] = df['SVTYPE'].apply(lambda val: f'<{val}>')

                if is_info_seq:
                    df['INFO'] = df.apply(lambda row: row['INFO'] + ';SEQ={SEQ}'.format(**row))

            else:

                # Check for sequence types that cannot be placed in ALT (may need symbolic ALTs)
                if svtype != 'inv':

                    df['REF'] = df.apply(lambda row:
                        ((row['REF'] + row['SEQ']) if row['POS'] > 0 else (row['SEQ'] + row['REF'])) if row['SVTYPE'] == 'DEL' else row['REF'],
                        axis=1
                    )

                    df['ALT'] = df.apply(lambda row:
                        ((row['REF'] + row['SEQ']) if row['POS'] > 0 else (row['SEQ'] + row['REF'])) if row['SVTYPE'] == 'INS' else row['REF'][0],
                        axis=1
                    )

                else:
                    raise NotImplementedError('No support for writing INV with REF/ALT (must use symbolic ALTs)')

                del df['SEQ']

                df['ALT'] = df['ALT'].apply(lambda val: val.upper())
                df['REF'] = df['REF'].apply(lambda val: val.upper())

        else:
            # Fix position for SNVs (0-based BED to 1-based VCF)
            df['POS'] += 1

            df['ALT'] = df['ALT'].apply(lambda val: val.upper())

        # QUAL
        if 'QUAL' not in df.columns:
            df['QUAL'] = '.'

        # Save columns needed for VCF
        df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'GT']]

        df_list.append(df)

    # Merge
    df = pd.concat(df_list, axis=0)
    df.sort_values(['#CHROM', 'POS'], inplace=True)

    # FILTER headers
    filter_header_list = [
        (filter, reason)
            for filter, reason in pavlib.call.FILTER_REASON.items()
    ]

    # INFO headers
    info_header_list = list()

    info_header_list.append(('ID', '1', 'String', 'Variant ID'))
    info_header_list.append(('SVTYPE', '1', 'String', 'Variant type'))
    info_header_list.append(('SVLEN', '.', 'Integer', 'Variant length'))
    info_header_list.append(('HAP', '.', 'STR', 'List of haplotype names variant was identified in'))
    info_header_list.append(('HAP_VARIANTS', '.', 'STR', 'List of variant IDs identifiying the variant merged in for each haplotype (INFO/HAP order)'))
    info_header_list.append(('COV_MEAN', '.', 'STR', 'Mean coverage for each haplotype under the whole variant (not just breakpoints, INFO/HAP order)'))
    info_header_list.append(('COV_PROP', '.', 'STR', 'Proportion of reference bases under the whole variant with at least one aligned query (assembly sequence) (INFO/HAP order)'))
    info_header_list.append(('QRY_REGION', '.', 'STR', 'Region of the query (assembly sequence) where this variant was found (chrom:pos-end, 1-based, closed coordinates, not BED) (INFO/HAP order)'))
    info_header_list.append(('QRY_STRAND', '.', 'STR', 'Orientation of the aligned query (assembly sequence) at this site (INFO/HAP order)'))
    info_header_list.append(('CALL_SOURCE', '.', 'STR', 'How variant was called - CIGAR, ALIGN_TRUNC, etc (INFO/HAP order)'))

    info_header_list.append(('INNER_REF', '.', 'String', 'Inversion inner breakpoint in reference coordinates (INFO/HAP order)'))
    info_header_list.append(('INNER_TIG', '.', 'String', 'Inversion inner breakpoint in contig coordinates (INFO/HAP order)'))
   # info_header_list.append(('HOM_REF', '.', 'String', 'Perfect breakpoint homology (SV sequence vs reference). Format \'X,Y\' where X homology upstream, and Y is homology downstream. Homology vs reference is often better for DEL.'))
   # info_header_list.append(('HOM_TIG', '.', 'String', 'Perfect breakpoint homology (SV sequence vs contig). Format \'X,Y\' where X homology upstream, and Y is homology downstream. Homology vs contig is often better for INS.'))

    if any_info_seq:
        info_header_list.append(('SEQ', '.', 'String', 'SV or indel sequence'))

    # ALT headers
    alt_header_list = list()

    unknown_alt = symbolic_alt_set - {'INS', 'DEL', 'INV'}

    if unknown_alt:
        raise RuntimeError(f'Unknown symbolic ALTs: {", ".join(sorted(unknown_alt))}')

    if symbolic_alt:
        if 'INS' in symbolic_alt_set:
            alt_header_list.append(('INS', 'Sequence insertion'))

        if 'DEL' in symbolic_alt_set:
            alt_header_list.append(('DEL', 'Sequence deletion'))

        if 'INV' in symbolic_alt_set:
            alt_header_list.append(('INV', 'Inversion'))

    # FORMAT
    df['FORMAT'] = 'GT'

    format_header_list = [
        ('GT', '1', 'String', 'Genotype')
    ]

    # VCF order
    df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'GT']]

    df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', asm_name]

    # Write
    with Bio.bgzf.open(output_filename, 'wt') as out_file:
        for line in svpoplib.vcf.header_list(
            df_ref,
            info_header_list,
            format_header_list,
            alt_header_list,
            filter_header_list,
            variant_source='PAV {}'.format(pavlib.constants.get_version_string()),
            ref_filename=os.path.basename(reference_filename)
        ):
            out_file.write(line)

        out_file.write('\t'.join(df.columns))
        out_file.write('\n')

        for index, row in df.iterrows():
            out_file.write('\t'.join(row.astype(str)))
            out_file.write('\n')
