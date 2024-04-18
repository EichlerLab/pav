"""
Call variants from aligned contigs.
"""

import collections
import numpy as np
import os
import pandas as pd
import sys

import Bio.bgzf
import Bio.SeqIO

import pavlib
import svpoplib

global ASM_TABLE
global MERGE_BATCH_COUNT
global REF_FA
global get_config

global expand
global intervaltree
global temp


#
# Merge haplotypes
#

# Generate all BED files
localrules: call_all_bed

rule call_all_bed:
    input:
        bed=lambda wildcards: pavlib.pipeline.expand_pattern(
            'results/{asm_name}/bed_merged/{filter}/{vartype_svtype}.bed.gz', ASM_TABLE,
            vartype_svtype=('svindel_ins', 'svindel_del', 'sv_inv', 'snv_snv'), filter=('pass', 'fail')
        )


# Concatenate variant BED files from batched merges - non-SNV (has variant FASTA).
rule call_merge_haplotypes:
    input:
        bed_batch=lambda wildcards: [
            'temp/{asm_name}/bed/batch/{filter}/{vartype_svtype}/{batch}.bed.gz'.format(
                asm_name=wildcards.asm_name, filter=wildcards.filter, vartype_svtype=wildcards.vartype_svtype, batch=batch
            ) for batch in range(MERGE_BATCH_COUNT)
        ]
    output:
        bed='results/{asm_name}/bed_merged/{filter}/{vartype_svtype}.bed.gz',
        fa='results/{asm_name}/bed_merged/{filter}/fa/{vartype_svtype}.fa.gz'
    wildcard_constraints:
        filter='pass|fail',
        vartype_svtype='svindel_ins|svindel_del|sv_inv'
    run:

        df_list = [pd.read_csv(file_name, sep='\t') for file_name in input.bed_batch if os.stat(file_name).st_size > 0]

        df = pd.concat(
            df_list, axis=0
        ).sort_values(
            ['#CHROM', 'POS', 'END', 'ID']
        )

        with Bio.bgzf.BgzfWriter(output.fa, 'wb') as out_file:
            Bio.SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')

        del df['SEQ']

        df.to_csv(
            output.bed, sep='\t', index=False, compression='gzip'
        )

# Concatenate variant BED files from batched merges - SNV (no variant FASTA).
rule call_merge_haplotypes_snv:
    input:
        bed_batch=lambda wildcards: [
            'temp/{asm_name}/bed/batch/{filter}/snv_snv/{batch}.bed.gz'.format(
                asm_name=wildcards.asm_name, filter=wildcards.filter, batch=batch
            ) for batch in range(MERGE_BATCH_COUNT)
        ]
    output:
        bed='results/{asm_name}/bed_merged/{filter}/snv_snv.bed.gz'
    wildcard_constraints:
        filter='pass|fail'
    run:

        df_list = [pd.read_csv(file_name, sep='\t') for file_name in input.bed_batch if os.stat(file_name).st_size > 0]

        df = pd.concat(
            df_list, axis=0
        ).sort_values(
            ['#CHROM', 'POS', 'END', 'ID']
        ).to_csv(
            output.bed, sep='\t', index=False, compression='gzip'
        )

# Merge by batches.
rule call_merge_haplotypes_batch:
    input:
        tsv='data/ref/merge_batch.tsv.gz',
        bed_var=lambda wildcards: [
            'results/{asm_name}/bed_hap/{filter}/{hap}/{vartype_svtype}.bed.gz'.format(hap=hap, **wildcards)
                for hap in pavlib.pipeline.get_hap_list(wildcards.asm_name, ASM_TABLE)
        ],
        fa_var=lambda wildcards: [
            'results/{asm_name}/bed_hap/{filter}/{hap}/fa/{vartype_svtype}.fa.gz'.format(hap=hap, **wildcards)
                for hap in pavlib.pipeline.get_hap_list(wildcards.asm_name, ASM_TABLE)
        ] if wildcards.vartype_svtype != 'snv_snv' else [],
        bed_callable=lambda wildcards: [
            'results/{asm_name}/callable/callable_regions_{hap}_500.bed.gz'.format(hap=hap, **wildcards)
                for hap in pavlib.pipeline.get_hap_list(wildcards.asm_name, ASM_TABLE)
        ]
    output:
        bed='temp/{asm_name}/bed/batch/{filter}/{vartype_svtype}/{batch}.bed.gz'
    wildcard_constraints:
        filter='pass|fail'
    threads: 12
    run:

        hap_list = pavlib.pipeline.get_hap_list(wildcards.asm_name, ASM_TABLE)

        # Read batch table
        df_batch = pd.read_csv(input.tsv, sep='\t')
        df_batch = df_batch.loc[df_batch['BATCH'] == int(wildcards.batch)]

        subset_chrom = set(df_batch['CHROM'])

        # Get variant type
        var_svtype_list = wildcards.vartype_svtype.split('_')

        if len(var_svtype_list) != 2:
            raise RuntimeError('Wildcard "vartype_svtype" must be two elements separated by an underscore: {}'.format(wildcards.var_svtype))

        # Read variant tables
        bed_list = list()

        for index in range(len(hap_list)):

            # Read variant table
            df = pd.read_csv(input.bed_var[index], sep='\t', dtype={'#CHROM': str}, low_memory=False)
            df = df.loc[df['#CHROM'].isin(subset_chrom)]

            df.set_index('ID', inplace=True, drop=False)
            df.index.name = 'INDEX'

            # Read SEQ
            if df.shape[0] > 0 and len(input.fa_var) > 0:
                df['SEQ'] = pd.Series({
                    record.id: str(record.seq)
                        for record in svpoplib.seq.fa_to_record_iter(
                            input.fa_var[index], record_set=set(df['ID'])
                        )
                })
            else:
                df['SEQ'] = np.nan

            bed_list.append(df)

        # Get configured merge definition
        config_def = pavlib.call.get_merge_params(wildcards, get_config(wildcards))

        print('Merging with def: ' + config_def)
        sys.stdout.flush()

        # Merge
        df = pavlib.call.merge_haplotypes(
            bed_list,
            input.bed_callable,
            hap_list,
            config_def,
            threads=threads
        )

        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


#
# Merge support
#

# Make a table of mappable regions by merging aligned loci with loci covered by alignment-truncating events.
# "flank" parameter is an integer describing how far away records may be to merge (similar to the "bedtools merge"
# "slop" parameter). The flank is not added to the regions that are output.
rule call_callable_regions:
    input:
        bed_align='results/{asm_name}/align/trim-tigref/aligned_tig_{hap}.bed.gz',
        bed_lg_del='temp/{asm_name}/lg_sv/sv_del_{hap}.bed.gz',
        bed_lg_ins='temp/{asm_name}/lg_sv/sv_ins_{hap}.bed.gz',
        bed_lg_inv='temp/{asm_name}/lg_sv/sv_inv_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/callable/callable_regions_{hap}_{flank}.bed.gz'
    run:

        # Get flank param
        try:
            flank = int(wildcards.flank)

        except ValueError:
            raise RuntimeError('Flank parameter is not an integer: {flank}'.format(**wildcards))

        # Merge
        df = pavlib.util.region_merge(
            [
                input.bed_align,
                input.bed_lg_del,
                input.bed_lg_ins,
                input.bed_lg_inv
            ],
            pad=flank
        )

        # Write
        df.to_csv(output.bed, sep='\t', index=False, compression='gzip')


#
# Integrate variant calls from multiple sources (per haplotype, pre-merge)
#

# Run all BED
localrules: call_all_bed_hap

rule call_all_bed_hap:
    input:
        bed=lambda wildcards: pavlib.pipeline.expand_pattern(
            'results/{asm_name}/bed_hap/{filter}/{hap}/{vartype_svtype}.bed.gz', ASM_TABLE,
            filter=('pass', 'fail'), vartype_svtype=('svindel_ins', 'svindel_del', 'sv_inv', 'snv_snv')
        ),
        fa=lambda wildcards: pavlib.pipeline.expand_pattern(
            'results/{asm_name}/bed_hap/fail/{hap}/fa/{vartype_svtype}.fa.gz', ASM_TABLE,
            vartype_svtype=('svindel_ins', 'svindel_del', 'sv_inv')
        ),
        fa_red=lambda wildcards: pavlib.pipeline.expand_pattern(
            'results/{asm_name}/bed_hap/fail/{hap}/redundant/fa/{vartype_svtype}.fa.gz', ASM_TABLE,
            vartype_svtype=('svindel_ins', 'svindel_del', 'sv_inv')
        )

# Create PASS BEDs
localrules: call_all_bed_pass

rule call_all_bed_pass:
    input:
        bed=lambda wildcards: pavlib.pipeline.expand_pattern(
            'results/{asm_name}/bed_hap/pass/{hap}/{vartype_svtype}.bed.gz', ASM_TABLE,
            vartype_svtype=('svindel_ins', 'svindel_del', 'sv_inv', 'snv_snv')
        )

# Write FASTA files for non-PASS variants (not in the redundant set)
rule call_integrate_fail_fa:
    input:
        bed='results/{asm_name}/bed_hap/fail/{hap}/{vartype_svtype}.bed.gz',
        fa='temp/{asm_name}/bed_hap/fail/{hap}/fa/{vartype_svtype}.fa.gz'
    output:
        fa='results/{asm_name}/bed_hap/fail/{hap}/fa/{vartype_svtype}.fa.gz'
    run:

        id_set = set(pd.read_csv(input.bed, sep='\t', usecols=['ID',])['ID'])

        with Bio.bgzf.BgzfWriter(output.fa) as out_file:
            Bio.SeqIO.write(
                svpoplib.seq.fa_to_record_iter(input.fa, record_set=id_set),
                out_file, 'fasta'
            )

# Write FASTA files for non-PASS variants (not in the redundant set)
rule call_integrate_fail_redundant_fa:
    input:
        bed='results/{asm_name}/bed_hap/fail/{hap}/redundant/{vartype_svtype}.bed.gz',
        fa='temp/{asm_name}/bed_hap/fail/{hap}/fa/{vartype_svtype}.fa.gz'
    output:
        fa='results/{asm_name}/bed_hap/fail/{hap}/redundant/fa/{vartype_svtype}.fa.gz'
    run:

        id_set = set(pd.read_csv(input.bed, sep='\t', usecols=['ID',])['ID'])

        with Bio.bgzf.BgzfWriter(output.fa) as out_file:
            Bio.SeqIO.write(
                svpoplib.seq.fa_to_record_iter(input.fa, record_set=id_set),
                out_file, 'fasta'
            )

# Separate multiple calls removed by the TRIM filter at the same site. Calls are separated to the "redundant"
# if they intersect a PASS variant. If there are multiple failed calls that do not intersect a PASS variant, then
# the one from the best alignment (highest QV, then longest alignment, then earliest start position) is chosen.
rule call_integrate_filter_redundant:
    input:
        bed='temp/{asm_name}/bed_hap/fail/{hap}/{vartype_svtype}.bed.gz',
        tsv='results/{asm_name}/bed_hap/fail/{hap}/redundant/intersect_{vartype_svtype}.tsv.gz'
    output:
        bed_nr='results/{asm_name}/bed_hap/fail/{hap}/{vartype_svtype}.bed.gz',
        bed_red='results/{asm_name}/bed_hap/fail/{hap}/redundant/{vartype_svtype}.bed.gz'
    run:

        # Read all FAIL variants
        df = pd.read_csv(input.bed, sep='\t', dtype={'#CHROM': str}, low_memory=False)

        # Initialize IDs acceptned into the NR (nonredundant) set with non-TRIM variants from df
        id_set = set(df.loc[df['FILTER'].apply(lambda val: 'TRIM' not in val.split(',')), 'ID'])

        # Get lead variant IDs for each set
        df_int = pd.read_csv(input.tsv, sep='\t', low_memory=False)

        df_int = df_int.loc[df_int['VARIANTS'].apply(lambda val: len(set(val.split(',')) & id_set) == 0)]  # Remove records that were already accepted into the NR set

        df_int = df_int.loc[df_int['SOURCE'].apply(lambda val: not val.startswith('PASS'))]  # Drop PASS records appearing in the PASS set

        id_set |= set(df_int['VARIANTS'].apply(lambda val: val.split(',')[0]))  # Choose the first of the merged set to represent all the other failed variants from other aligned segments at this site

        # Split variants and write
        df_nr = df.loc[df['ID'].isin(id_set)]

        index_set = set(df_nr.index)
        df_red = df.loc[[index not in index_set for index in df.index]]

        df_nr.to_csv(output.bed_nr, sep='\t', index=False, compression='gzip')
        df_red.to_csv(output.bed_red, sep='\t', index=False, compression='gzip')


# Concatenate variant BED files from batched merges.
rule call_intersect_fail:
    input:
        tsv=lambda wildcards: [
            'temp/{asm_name}/bed_hap/fail/{hap}/intersect/{vartype_svtype}_{batch}.tsv.gz'.format(
                asm_name=wildcards.asm_name, hap=wildcards.hap, vartype_svtype=wildcards.vartype_svtype, batch=batch
            ) for batch in range(MERGE_BATCH_COUNT)
        ]
    output:
        tsv='results/{asm_name}/bed_hap/fail/{hap}/redundant/intersect_{vartype_svtype}.tsv.gz'
    run:

        df_list = [pd.read_csv(file_name, sep='\t') for file_name in input.tsv if os.stat(file_name).st_size > 0]

        df = pd.concat(
            df_list, axis=0
        ).to_csv(
            output.tsv, sep='\t', index=False, compression='gzip'
        )


# Intersect failed calls with FILTER=PASS calls and other failed calls. Used to eliminate redundant failed call
# annotations.
rule call_intersect_fail_batch:
    input:
        bed_pass='results/{asm_name}/bed_hap/pass/{hap}/{vartype_svtype}.bed.gz',
        bed_fail='temp/{asm_name}/bed_hap/fail/{hap}/{vartype_svtype}.bed.gz',
        fa_pass=lambda wildcards: ['results/{asm_name}/bed_hap/pass/{hap}/fa/{vartype_svtype}.fa.gz']
            if wildcards.vartype_svtype != 'snv_snv' else [],
        fa_fail=lambda wildcards: ['temp/{asm_name}/bed_hap/fail/{hap}/fa/{vartype_svtype}.fa.gz']
            if wildcards.vartype_svtype != 'snv_snv' else [],
        bed_align='results/{asm_name}/align/trim-none/aligned_tig_{hap}.bed.gz',
        tsv_batch='data/ref/merge_batch.tsv.gz'
    output:
        tsv=temp('temp/{asm_name}/bed_hap/fail/{hap}/intersect/{vartype_svtype}_{batch}.tsv.gz')
    threads: 1
    run:

        # Get chromosome set for this batch
        df_batch = pd.read_csv(input.tsv_batch, sep='\t', low_memory=False)
        df_batch = df_batch.loc[df_batch['BATCH'] == int(wildcards.batch)]

        subset_chrom = set(df_batch['CHROM'])

        del df_batch

        if len(subset_chrom) > 0:
            # Read
            df_pass = pd.read_csv(
                input.bed_pass, sep='\t', dtype={'#CHROM': str, 'ALIGN_INDEX': str},
                usecols=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'ALIGN_INDEX', 'FILTER'] + (['ALT', 'REF'] if wildcards.vartype_svtype == 'snv_snv' else []),
                low_memory=False
            )

            df_pass = df_pass.loc[df_pass['#CHROM'].isin(subset_chrom)].copy()

            df_fail = pd.read_csv(
                input.bed_fail, sep='\t',  dtype={'#CHROM': str, 'ALIGN_INDEX': str},
                usecols=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'ALIGN_INDEX', 'FILTER'] + (['ALT', 'REF'] if wildcards.vartype_svtype == 'snv_snv' else []),
                low_memory=False
            )

            df_fail = df_fail.loc[df_fail['#CHROM'].isin(subset_chrom)].copy()

            # Add SEQ
            df_pass.set_index('ID', inplace=True, drop=False)
            df_pass.index.name = 'INDEX'

            if df_pass.shape[0] > 0 and len(input.fa_pass) > 0:
                df_pass['SEQ'] = pd.Series({
                    record.id: str(record.seq)
                        for record in svpoplib.seq.fa_to_record_iter(
                            input.fa_pass[0], record_set=set(df_pass['ID'])
                        )
                })
            else:
                df_pass['SEQ'] = np.nan

            df_fail.set_index('ID', inplace=True, drop=False)
            df_fail.index.name = 'INDEX'

            if df_fail.shape[0] > 0 and len(input.fa_fail) > 0:
                df_fail['SEQ'] = pd.Series({
                    record.id: str(record.seq)
                        for record in svpoplib.seq.fa_to_record_iter(
                            input.fa_fail[0], record_set=set(df_fail['ID'])
                        )
                })
            else:
                df_fail['SEQ'] = np.nan


            # Split df_fail into TRIM and non-TRIM, append non-TRIM to df_pass
            df_fail_trim = df_fail.loc[df_fail['FILTER'].apply(lambda filter_val: 'TRIM' in filter_val)]
            index_set = set(df_fail_trim.index)
            df_fail_notrim = df_fail.loc[[index not in index_set for index in df_fail.index]]

            if df_fail_notrim.shape[0] > 0:
                if df_pass.shape[0] > 0:
                    df_pass = pd.concat([df_pass, df_fail_notrim], axis=0)
                else:
                    df_pass = df_fail_notrim

            df_fail = df_fail_trim

            # Prioritize alignments by MAPQ and length
            df_align = pd.read_csv(
                input.bed_align, sep='\t',
                usecols=['INDEX', 'QRY_LEN', 'QRY_POS', 'QRY_END', 'MAPQ'],
                dtype={'INDEX': int}
            )

            index_set = {
                int(val) for val_str in df_fail['ALIGN_INDEX'] for val in str(val_str).split(',')
            }

            df_align = df_align.loc[df_align['INDEX'].isin(index_set)]

            df_align['LEN'] = df_align['QRY_END'] - df_align['QRY_POS']

            index_list = list(index for index in df_align.sort_values(['MAPQ', 'LEN', 'INDEX'])['INDEX'])

            # For variants with multiple alignment indices, choose the least one.
            for index_var in df_fail.index:
                align_index_set = {int(val) for val in df_fail.loc[index_var, 'ALIGN_INDEX'].split(',')}
                df_fail.loc[index_var, 'ALIGN_INDEX'] = [val for val in index_list if val in align_index_set][-1]

            index_set = set(df_fail['ALIGN_INDEX'])
            index_list = [val for val in index_list if val in index_set]

            # Separate failed variants on alignment index
            df_list = [(df_pass, 'PASS')]

            for index in index_list:
                df_list.append((df_fail.loc[df_fail['ALIGN_INDEX'] == index], f'TRIM_{index}'))

            # Intersect
            df_intersect = svpoplib.svmerge.merge_variants(
                bed_list=[val[0] for val in df_list],
                sample_names=[val[1] for val in df_list],
                strategy=pavlib.call.get_merge_params(wildcards, get_config(wildcards)),
                threads=threads
            )

            # Reshape and write
            df_intersect = df_intersect[['ID'] + [col for col in df_intersect.columns if col.startswith('MERGE_')]]

            del df_intersect['MERGE_SRC']
            del df_intersect['MERGE_SRC_ID']

            df_intersect['MERGE_SAMPLES'] = df_intersect['MERGE_SAMPLES'].apply(lambda val_list:
                ','.join([
                    (val[5:] if val.startswith('TRIM_') else val) for val in val_list.split(',')
                ])
            )

            df_intersect.columns = [col[6:] if col.startswith('MERGE_') else col for col in df_intersect.columns]
            df_intersect.columns = ['SOURCE' if col == 'SAMPLES' else col for col in df_intersect.columns]

            df_intersect.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

        else:
            # Write empty file (skipped when batches are concatenated)
            with open(output.tsv, 'wt') as out_file:
                pass

# Filter variants from inside inversions
rule call_integrate_sources:
    input:
        bed_cigar_insdel='temp/{asm_name}/cigar/merged/svindel_insdel_{hap}.bed.gz',
        bed_cigar_snv='temp/{asm_name}/cigar/merged/snv_snv_{hap}.bed.gz',
        bed_lg_ins='temp/{asm_name}/lg_sv/sv_ins_{hap}.bed.gz',
        bed_lg_del='temp/{asm_name}/lg_sv/sv_del_{hap}.bed.gz',
        bed_lg_inv='temp/{asm_name}/lg_sv/sv_inv_{hap}.bed.gz',
        bed_inv='temp/{asm_name}/inv_caller/sv_inv_{hap}.bed.gz',
        bed_cov_tig='results/{asm_name}/align/trim-tig/depth_tig_{hap}.bed.gz',
        bed_filter=lambda wildcards: pavlib.pipeline.get_asm_config(
            wildcards.asm_name, wildcards.hap, ASM_TABLE, config
        )['filter_input']
    output:
        bed_ins_pass='results/{asm_name}/bed_hap/pass/{hap}/svindel_ins.bed.gz',
        bed_del_pass='results/{asm_name}/bed_hap/pass/{hap}/svindel_del.bed.gz',
        bed_inv_pass='results/{asm_name}/bed_hap/pass/{hap}/sv_inv.bed.gz',
        bed_snv_pass='results/{asm_name}/bed_hap/pass/{hap}/snv_snv.bed.gz',
        fa_ins_pass='results/{asm_name}/bed_hap/pass/{hap}/fa/svindel_ins.fa.gz',
        fa_del_pass='results/{asm_name}/bed_hap/pass/{hap}/fa/svindel_del.fa.gz',
        fa_inv_pass='results/{asm_name}/bed_hap/pass/{hap}/fa/sv_inv.fa.gz',
        bed_ins_fail=temp('temp/{asm_name}/bed_hap/fail/{hap}/svindel_ins.bed.gz'),
        bed_del_fail=temp('temp/{asm_name}/bed_hap/fail/{hap}/svindel_del.bed.gz'),
        bed_inv_fail=temp('temp/{asm_name}/bed_hap/fail/{hap}/sv_inv.bed.gz'),
        bed_snv_fail=temp('temp/{asm_name}/bed_hap/fail/{hap}/snv_snv.bed.gz'),
        fa_ins_fail=temp('temp/{asm_name}/bed_hap/fail/{hap}/fa/svindel_ins.fa.gz'),
        fa_del_fail=temp('temp/{asm_name}/bed_hap/fail/{hap}/fa/svindel_del.fa.gz'),
        fa_inv_fail=temp('temp/{asm_name}/bed_hap/fail/{hap}/fa/sv_inv.fa.gz')
    run:

        # Get parameters
        local_config = get_config(wildcards)

        inv_min = local_config.get('inv_min', 0)
        inv_max = local_config.get('inv_max', 1e10)
        inv_inner = local_config.get('inv_inner', False)
        redundant_callset = pavlib.util.as_bool(local_config.get('redundant_callset', False))

        if inv_min is not None and inv_min != 'unlimited':
            inv_min = int(inv_min)
        else:
            inv_min = None

        if inv_max is not None and inv_max != 'unlimited':
            inv_max = int(inv_max)
        else:
            inv_max = None

        # Read tig filter (if present)
        qry_filter_tree = None

        if len(input.bed_filter) > 0:
            qry_filter_tree = collections.defaultdict(intervaltree.IntervalTree)

            for filter_filename in input.bed_filter:
                df_filter = pd.read_csv(filter_filename, sep='\t', header=None, comment='#', usecols=(0, 1, 2))
                df_filter.columns = ['#CHROM', 'POS', 'END']

                for index, row in df_filter.iterrows():
                    qry_filter_tree[row['#CHROM']][row['POS']:row['END']] = True

        # Create large variant filter tree (for filtering small variants inside larger ones)
        compound_filter_tree = collections.defaultdict(intervaltree.IntervalTree)

        # Tuple of (pass, drop) filenames for each variant type
        out_filename_dict = {
            'inv': (output.bed_inv_pass, output.bed_inv_fail, output.fa_inv_pass, output.fa_inv_fail),
            'ins': (output.bed_ins_pass, output.bed_ins_fail, output.fa_ins_pass, output.fa_ins_fail),
            'del': (output.bed_del_pass, output.bed_del_fail, output.fa_del_pass, output.fa_del_fail),
            'snv': (output.bed_snv_pass, output.bed_snv_fail, None, None),
        }

        #
        # Integrate and filter variants
        #

        df_insdel_list = list()  # Collect ins/del tables into this list until they are merged and written

        for vartype in ('inv', 'lg_del', 'lg_ins', 'insdel', 'snv'):

            # Read
            do_write = True      # Write variant call table. Set to False for INS/DEL variants until they are all collected into df_insdel_list.
            is_insdel = False    # Variant is an INS/DEL type, append to df_insdel_list (list is merged and written when both is_insdel and do_write are True).
            is_inv = False       # Variant is an inversion, apply inversion filtering steps.
            add_compound = True  # Add variant regions to the compound filter if True. Does not need to be set for small variants (just consumes CPU time and memory).

            if vartype == 'inv':
                df, filter_dict, compound_dict = pavlib.call.read_variant_table([input.bed_inv, input.bed_lg_inv], True)
                is_inv = True

            elif vartype == 'lg_del':
                df, filter_dict, compound_dict = pavlib.call.read_variant_table(input.bed_lg_del, True)
                do_write = False
                is_insdel = True

            elif vartype == 'lg_ins':
                df, filter_dict, compound_dict = pavlib.call.read_variant_table(input.bed_lg_ins, True)
                do_write = False
                is_insdel = True

            elif vartype == 'insdel':
                df, filter_dict, compound_dict = pavlib.call.read_variant_table(input.bed_cigar_insdel, True)
                do_write = True
                is_insdel = True
                add_compound = False

            elif vartype == 'snv':
                df, filter_dict, compound_dict = pavlib.call.read_variant_table(input.bed_cigar_snv, True)
                add_compound = False

            else:
                assert False, f'vartype in control loop does not match a known value: {vartype}'

            # Apply filters
            if df.shape[0] > 0:

                # Apply query filtered regions
                pavlib.call.apply_qry_filter_tree(df, qry_filter_tree, filter_dict)

                # SVLEN min
                if is_inv and inv_min is not None:
                    for index in df.index[df['SVLEN'] < inv_min]:
                        filter_dict[index].add('SVLEN')

                # SVLEN max
                if is_inv and inv_max is not None:
                    for index in df.index[df['SVLEN'] > inv_max]:
                        filter_dict[index].add('SVLEN')

                # Filter compound
                if not (redundant_callset or is_inv and inv_inner):
                    pavlib.call.apply_compound_filter(df, compound_filter_tree, filter_dict, compound_dict, add_compound)

            # Compound filter
            pavlib.call.update_filter_compound_fields(df, filter_dict, compound_dict)

            del(filter_dict)
            del(compound_dict)

            # Alignment depth
            df['COV_MEAN'] = np.nan
            df['COV_PROP'] = np.nan
            df['COV_QRY'] = ''

            depth_container = pavlib.call.DepthContainer(pd.read_csv(input.bed_cov_tig, sep='\t', dtype={'#CHROM': str}))

            for index, row in df.iterrows():
                df.loc[index, 'COV_MEAN'], df.loc[index, 'COV_PROP'], df.loc[index, 'COV_QRY'] = depth_container.get_depth(row)

            del depth_container

            # Version variant IDs prioritizing PASS over non-PASS (avoid version suffixes on PASS).
            df['ID'] = pavlib.call.version_variant_bed_id(df)

            # Add to df_insdel_list
            if is_insdel:
                df_insdel_list.append(df)

            # Write
            if do_write:
                if is_insdel:
                    # Merge
                    df = pd.concat(df_insdel_list, axis=0).sort_values(['#CHROM', 'POS'])

                    del(df_insdel_list)

                    # Separate INS and DEL, then write
                    for svtype in ('ins', 'del'):

                        # Pass
                        subdf = df.loc[
                            (df['SVTYPE'] == svtype.upper()) & (df['FILTER'] == 'PASS')
                        ]

                        with Bio.bgzf.BgzfWriter(out_filename_dict[svtype][2]) as out_file:
                            Bio.SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(subdf), out_file, 'fasta')

                        del(subdf['SEQ'])

                        subdf.to_csv(out_filename_dict[svtype][0], sep='\t', index=False, compression='gzip')

                        # Fail
                        subdf = df.loc[
                            (df['SVTYPE'] == svtype.upper()) & (df['FILTER'] != 'PASS')
                        ]

                        with Bio.bgzf.BgzfWriter(out_filename_dict[svtype][3]) as out_file:
                            Bio.SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(subdf), out_file, 'fasta')

                        del(subdf['SEQ'])

                        subdf.to_csv(
                            out_filename_dict[svtype][1], sep='\t', index=False, compression='gzip'
                        )

                        del(subdf)

                else:

                    # Pass

                    subdf = df.loc[df['FILTER'] == 'PASS']

                    if out_filename_dict[vartype][2] is not None:
                        with Bio.bgzf.BgzfWriter(out_filename_dict[vartype][2]) as out_file:
                            Bio.SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(subdf), out_file, 'fasta')

                        del(subdf['SEQ'])

                    subdf.to_csv(out_filename_dict[vartype][0], sep='\t', index=False, compression='gzip')

                    # Fail
                    subdf = df.loc[df['FILTER'] != 'PASS']

                    if out_filename_dict[vartype][3] is not None:
                        with Bio.bgzf.BgzfWriter(out_filename_dict[vartype][3]) as out_file:
                            Bio.SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(subdf), out_file, 'fasta')

                        del(subdf['SEQ'])

                    subdf.to_csv(
                        out_filename_dict[vartype][1], sep='\t', index=False, compression='gzip'
                    )

                    del(subdf)

            # Clean
            del(df)



#
# Call from CIGAR
#

# Merge discovery sets from each batch.
rule call_cigar_merge:
    input:
        bed_insdel=expand('temp/{{asm_name}}/cigar/batched/insdel_{{hap}}_{batch}.bed.gz', batch=range(pavlib.cigarcall.CALL_CIGAR_BATCH_COUNT)),
        bed_snv=expand('temp/{{asm_name}}/cigar/batched/snv.bed_{{hap}}_{batch}.gz', batch=range(pavlib.cigarcall.CALL_CIGAR_BATCH_COUNT))
    output:
        bed_insdel=temp('temp/{asm_name}/cigar/merged/svindel_insdel_{hap}.bed.gz'),
        bed_snv=temp('temp/{asm_name}/cigar/merged/snv_snv_{hap}.bed.gz')
    run:

        # INS/DEL
        df_insdel = pd.concat(
            [pd.read_csv(file_name, sep='\t', keep_default_na=False) for file_name in input.bed_insdel],
            axis=0
        ).reset_index(drop=True)

        df_insdel.sort_values(
            ['#CHROM', 'POS', 'END', 'ID']
        ).to_csv(
            output.bed_insdel, sep='\t', index=False, compression='gzip'
        )

        # SNV
        df_snv = pd.concat(
            [pd.read_csv(file_name, sep='\t', keep_default_na=False) for file_name in input.bed_snv],
            axis=0
        ).reset_index(drop=True)

        df_snv.sort_values(
            ['#CHROM', 'POS']
        ).to_csv(
            output.bed_snv, sep='\t', index=False, compression='gzip'
        )


# Call variants by alignment CIGAR parsing.
#
# IDs are not versioned by this rule, versioning must be applied after batches are merged.
rule call_cigar:
    input:
        bed='results/{asm_name}/align/trim-none/aligned_tig_{hap}.bed.gz',
        bed_trim='results/{asm_name}/align/trim-tigref/aligned_tig_{hap}.bed.gz',
        tig_fa_name='temp/{asm_name}/align/contigs_{hap}.fa.gz'
    output:
        bed_insdel=temp('temp/{asm_name}/cigar/batched/insdel_{hap}_{batch}.bed.gz'),
        bed_snv=temp('temp/{asm_name}/cigar/batched/snv.bed_{hap}_{batch}.gz')
    run:

        batch = int(wildcards.batch)

        # Read
        df_align = pd.read_csv(input.bed, sep='\t', dtype={'#CHROM': str}, keep_default_na=False, low_memory=False)

        df_align = df_align.loc[df_align['CALL_BATCH'] == batch]

        # Call
        df_snv, df_insdel = pavlib.cigarcall.make_insdel_snv_calls(df_align, REF_FA, input.tig_fa_name, wildcards.hap, version_id=False)

        # Read trimmed alignments
        df_trim = pd.read_csv(
            input.bed_trim, sep='\t', usecols=['POS', 'END', 'INDEX'],
            index_col='INDEX'
        ).astype(int)

        # Set FILTER - snv
        df_pass = df_trim.reindex(
            list(df_snv['ALIGN_INDEX']), fill_value=-1
        ).set_index(
            df_snv.index, drop=True
        )

        df_snv['FILTER'] = (
                (df_snv['POS'] > df_pass['POS']) & (df_snv['END'] < df_pass['END'])
        ).apply(
            lambda val: 'PASS' if val else 'TRIM'
        )

        # Set FILTER - insdel
        df_pass = df_trim.reindex(
            list(df_insdel['ALIGN_INDEX']), fill_value=-1
        ).set_index(
            df_insdel.index, drop=True
        )

        df_insdel['FILTER'] = (
                (df_insdel['POS'] > df_pass['POS']) & (df_insdel['END'] < df_pass['END'])
        ).apply(
            lambda val: 'PASS' if val else 'TRIM'
        )

        # Write
        df_insdel.to_csv(output.bed_insdel, sep='\t', index=False, compression='gzip')
        df_snv.to_csv(output.bed_snv, sep='\t', index=False, compression='gzip')


#
# Supporting rules
#

# Create a table of merge batch assignments
localrules: call_merge_batch_table

rule call_merge_batch_table:
    input:
        tsv='data/ref/contig_info.tsv.gz'
    output:
        tsv='data/ref/merge_batch.tsv.gz'
    run:

        # Read and sort
        df = pd.read_csv(
            input.tsv, sep='\t', dtype={'CHROM': str, 'LEN': int}
        ).sort_values(
            'LEN', ascending=False
        ).set_index(
            'CHROM'
        )[['LEN']]

        df['BATCH'] = -1

        # Get a list of assignments for each batch
        list_chr = collections.defaultdict(list)
        list_size = collections.Counter()

        def get_smallest():
            """
            Get the next smallest bin.
            """

            min_index = 0

            for i in range(MERGE_BATCH_COUNT):

                if list_size[i] == 0:
                    return i

                if list_size[i] < list_size[min_index]:
                    min_index = i

            return min_index

        for chrom in df.index:
            i = get_smallest()
            df.loc[chrom, 'BATCH'] = i
            list_size[i] += df.loc[chrom, 'LEN']

        # Check
        if np.any(df['BATCH'] < 0):
            raise RuntimeError('Failed to assign all reference contigs to batches (PROGRAM BUG)')

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')
