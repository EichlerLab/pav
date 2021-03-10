"""
Rules for writing VCF output.
"""

_VCF_SVTYPE = {
    'sv': ('ins', 'del', 'inv'),
    'indel': ('ins', 'del'),
    'snv': ('snv', )
}


# vcf_write_vcf
#
# Make VCF headers.
rule vcf_write_vcf:
    input:
        bed_snv_snv='results/{asm_name}/bed/snv_snv.bed.gz',
        bed_indel_ins='results/{asm_name}/bed/indel_ins.bed.gz',
        bed_indel_del='results/{asm_name}/bed/indel_del.bed.gz',
        bed_sv_ins='results/{asm_name}/bed/sv_ins.bed.gz',
        bed_sv_del='results/{asm_name}/bed/sv_del.bed.gz',
        bed_sv_inv='results/{asm_name}/bed/sv_inv.bed.gz',
        fa_indel_ins='results/{asm_name}/bed/fa/indel_ins.fa.gz',
        fa_indel_del='results/{asm_name}/bed/fa/indel_del.fa.gz',
        fa_sv_ins='results/{asm_name}/bed/fa/sv_ins.fa.gz',
        fa_sv_del='results/{asm_name}/bed/fa/sv_del.fa.gz',
        fa_sv_inv='results/{asm_name}/bed/fa/sv_inv.fa.gz',
        ref_tsv='data/ref/contig_info.tsv.gz'
    output:
        vcf='pav_{asm_name}.vcf.gz'
    wildcard_constraints:
        alt_fmt='alt|sym'
    run:

        # Check assembly name
        if wildcards.asm_name in {'#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'}:
            raise RuntimeError(f'Assembly name conflicts with a VCF header column name: {wildcards.asm_name}')

        # Check alt format
        symbolic_alt = False
        # if wildcards.alt_fmt == 'alt':
        #     symbolic_alt = False
        # elif wildcards.alt_fmt == 'sym':
        #     symbolic_alt = True
        # else:
        #     raise RuntimeError(f'Unknown alt format wildcard (alt_fmt): {alt_fmt}')

        # Process variant types
        df_list = list()

        for vartype in ('sv', 'indel', 'snv'):
            for svtype in _VCF_SVTYPE[vartype]:
                print(f'{vartype} - {svtype}')  # DBGTMP

                # Read variants
                bed_file_name = input[f'bed_{vartype}_{svtype}']
                df = pd.read_csv(bed_file_name, sep='\t')

                df['VARTYPE'] = vartype.upper()

                # Read sequence
                if vartype in ('sv', 'indel'):

                    df.set_index('ID', inplace=True)

                    # Get Sequence from FASTA and assign to SEQ column (not for SNVs)
                    fa_file_name = input[f'fa_{vartype}_{svtype}']
                    with gzip.open(fa_file_name, 'rt') as fa_in:
                        df_seq_dict = {
                            record.name: str(record.seq) for record in SeqIO.parse(fa_in, 'fasta')
                        }

                    df['SEQ'] = pd.Series(df_seq_dict)

                    del(df_seq_dict)

                    df.reset_index(inplace=True)

                    # Check for missing sequences
                    df_null_seq = df.loc[pd.isnull(df['SEQ'])]

                    if df_null_seq.shape[0] > 0:
                        id_list = ', '.join(df_null_seq.iloc[:3]['ID'])

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
                for col in ('CALL_SOURCE', 'TIG_REGION', 'QUERY_STRAND', 'RGN_REF_INNER', 'RGN_TIG_INNER'):
                    if col in df.columns:
                        df[col] = df[col].apply(lambda val: val.replace(';', ','))

                if svtype == 'del':
                    df['SVLEN'] = - np.abs(df['SVLEN'])

                # INFO: Base
                df['INFO'] = df.apply(lambda row: 'ID={ID};SVTYPE={SVTYPE}'.format(**row), axis=1)

                # INFO: Add SV/INDEL annotations
                if vartype != 'snv':
                    df['INFO'] = df.apply(lambda row: row['INFO'] + ';SVLEN={SVLEN}'.format(**row), axis=1)

                # INFO: Add contig placement info
                df['INFO'] = df.apply(lambda row: row['INFO'] + ';TIG_REGION={TIG_REGION};QUERY_STRAND={QUERY_STRAND}'.format(**row), axis=1)

                # INFO: Add INV
                if svtype == 'inv':
                    df['INFO'] = df.apply(lambda row: row['INFO'] + ';INNER_REF={RGN_REF_INNER};INNER_TIG={RGN_TIG_INNER}'.format(**row), axis=1)

                # REF
                if 'REF' not in df.columns:
                    df['REF'] = list(svpoplib.vcf.ref_base(df, config['reference']))

                if svtype == 'inv' and not symbolic_alt:
                    df_ref_base = df['REF']
                    df['REF'] = df_ref_base + svpoplib.ref.get_ref_region(df, config['reference']).apply(lambda val: val.upper())

                # ALT
                if vartype != 'snv':
                    if symbolic_alt:
                        df['ALT'] = df['SVTYPE'].apply(lambda val: f'<{val}>')

                        df['INFO'] = df.apply(lambda row: row['INFO'] + ';SEQ={SEQ}'.format(**row))

                    else:

                        # Check for sequence types that cannot be placed in ALT (may need symbolic ALTs)
                        if svtype != 'inv':
                            df['REF'] = df.apply(lambda row: row['REF'] + row['SEQ'] if row['SVTYPE'] == 'DEL' else row['REF'], axis=1)
                            df['ALT'] = df.apply(lambda row: row['REF'] + row['SEQ'] if row['SVTYPE'] == 'INS' else row['REF'][0], axis=1)
                        else:
                            df['ALT'] = df_ref_base + df['SEQ']

                        del df['SEQ']
                else:
                    # Fix position for SNVs (0-based BED to 1-based VCF)
                    df['POS'] += 1

                # Save columns needed for VCF
                df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO', 'GT']]

                df_list.append(df)

        # Merge
        df = pd.concat(df_list, axis=0)
        df.sort_values(['#CHROM', 'POS'], inplace=True)

        # INFO headers
        info_header_list = list()

        info_header_list.append(('ID', '1', 'String', 'Variant ID'))
        info_header_list.append(('SVTYPE', '1', 'String', 'Variant type'))
        info_header_list.append(('SVLEN', '.', 'String', 'Variant length'))
        info_header_list.append(('TIG_REGION', '.', 'String', 'Contig region where variant was found (one per alt with h1 before h2 for homozygous calls)'))
        info_header_list.append(('QUERY_STRAND', '.', 'String', 'Strand of variant in the contig relative to the reference (order follows TIG_REGION)'))
        info_header_list.append(('INNER_REF', '.', 'String', 'Inversion inner breakpoint in reference coordinates (order follows TIG_REGION)'))
        info_header_list.append(('INNER_TIG', '.', 'String', 'Inversion inner breakpoint in contig coordinates (order follows TIG_REGION)'))

        if symbolic_alt:
            info_header_list.append(('SEQ', '.', 'String', 'SV or indel sequence'))

        # ALT headers
        alt_header_list = list()

        if symbolic_alt:
            alt_header_list.append(('INS', 'Sequence insertion'))
            alt_header_list.append(('DEL', 'Sequence deletion'))
            alt_header_list.append(('INV', 'Inversion'))

        # QUAL, FILTER, FORMAT
        if 'QUAL' not in df.columns:
            df['QUAL'] = '.'

        if 'FILTER' in df.columns:
            raise RuntimeError('FILTER is defined in dataframe, but FILTER headers are not yet implemented')

        filter_header_list = list()
        df['FILTER'] = '.'

        df['FORMAT'] = 'GT'

        format_header_list = [
            ('GT', '1', 'String', 'Genotype')
        ]

        # VCF order
        df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'GT']]

        df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', wildcards.asm_name]

        # Write
        df_ref = pd.read_csv(input.ref_tsv, sep='\t')

        with Bio.bgzf.open(output.vcf, 'wt') as out_file:
            for line in svpoplib.vcf.header_list(
                df_ref,
                info_header_list,
                format_header_list,
                alt_header_list,
                filter_header_list,
                variant_source='PAV {}'.format(pavlib.constants.get_version_string()),
                ref_file_name=os.path.basename(config.get('reference'))
            ):
                out_file.write(line)

            out_file.write('\t'.join(df.columns))
            out_file.write('\n')

            for index, row in df.iterrows():
                out_file.write('\t'.join(row.astype(str)))
                out_file.write('\n')

        # Write tabix index if possible
        try:
            shell("""tabix {output.vcf} && touch -r {output.vcf} {output.vcf}.tbi""")
        except:
            pass
