"""
Data files including reference and data tables for the reference.
"""

import collections
import gzip
import numpy as np
import os
import pandas as pd

import Bio.SeqIO

import pavlib
import svpoplib

global shell


#
# Pre-run targets
#

def data_init_targets(wildcards=None):
    """
    Get a list of input files to be generated before running samples. This target can be used to setup so that each
    sample can be run independently.

    :param wildcards: Ignored. Function signature needed for Snakemake input function.

    :return: List of run targets.
    """

    aligner_set = set()
    part_set = set()

    # Get a set of aligners and partitions
    for asm_name in ASM_TABLE.index:

        local_config = pavlib.config.get_override_config(asm_name, config, ASM_TABLE)

        aligner_set.add(pavlib.config.get_aligner(asm_name, config, ASM_TABLE))

        part_set.add(
            pavlib.config.CONFIG_PARAM_DICT['cigar_partitions'].get_value(
                local_config['cigar_partitions'] if 'cigar_partitions'in local_config else None
            )
        )

        part_set.add(
            pavlib.config.CONFIG_PARAM_DICT['merge_partitions'].get_value(
                local_config['merge_partitions'] if 'cigar_partitions'in local_config else None
            )
        )

    # Construct target list
    target_list = [
        'data/ref/contig_info.tsv.gz',
        'data/ref/n_gap.bed.gz',
        'data/ref/ref.fa.gz',
        'data/ref/ref.fa.gz.fai'
    ]

    if 'lra' in aligner_set:
        target_list += [
            'data/ref/ref.fa.gz.gli',
            'data/ref/ref.fa.gz.mms'
        ]

    for part in part_set:
        target_list.append(
            f'data/ref/partition_{part}.tsv.gz'
        )

    return target_list



#
# Pre-run targets
#

localrules: data_init

# Generate all pre-target runs
rule data_init:
    input:
        files=data_init_targets


#
# Reference contig data table
#

# Contig table.
rule data_ref_contig_table:
    input:
        ref_fa='data/ref/ref.fa.gz'
    output:
        tsv='data/ref/contig_info.tsv.gz'
    run:

        svpoplib.ref.get_ref_info(
            input.ref_fa
        ).to_csv(
            output.tsv, sep='\t', index=True, compression='gzip'
        )


#
# Reference annotation
#

# Find locations of N-gaps.
rule data_align_ref_anno_n_gap:
    input:
        ref_fa='data/ref/ref.fa.gz'
    output:
        bed='data/ref/n_gap.bed.gz'
    run:

        with gzip.open(output.bed, 'wt') as out_file:
            out_file.write('#CHROM\tPOS\tEND\n')

            with gzip.open(input.ref_fa, 'rt') as in_file:
                for record in Bio.SeqIO.parse(in_file, 'fasta'):

                    pos = None
                    end = None

                    enum_list = [i for i, val in enumerate(str(record.seq).upper()) if val == 'N']

                    for index in enum_list:

                        if pos is None:
                            pos = end = index

                        elif index == end + 1:
                            end = index

                        else:
                            out_file.write(f'{record.id}\t{pos}\t{end + 1}\n')
                            pos = end = index

                    if pos is not None:
                        out_file.write(f'{record.id}\t{pos}\t{end + 1}\n')


#
# Reference
#

# Index reference for LRA.
rule data_align_ref_lra_index:
    input:
        fa='data/ref/ref.fa.gz'
    output:
        gli='data/ref/ref.fa.gz.gli',
        mmi='data/ref/ref.fa.gz.mms'
    shell:
        """lra index -CONTIG {input.fa}"""

# Setup reference.
rule data_align_ref:
    input:
        ref_fa=config['reference']
    output:
        ref_fa='data/ref/ref.fa.gz',
        ref_fai='data/ref/ref.fa.gz.fai'
    run:

        is_copy = False

        if os.stat(input.ref_fa).st_size == 0:
            raise RuntimeError(f'Empty input FASTA file: {input.ref_fa}')

        # Link or copy FASTA to FA/GZ
        if os.path.basename(input.ref_fa).endswith('.gz'):
            # Link if gzipped
            os.symlink(os.path.abspath(input.ref_fa), output.ref_fa)

        else:
            # Copy if not gzipped
            pavlib.pipeline.input_tuples_to_fasta([(input.ref_fa, 'fasta')], output.ref_fa)
            is_copy = True

        # Index
        if not is_copy and os.path.isfile(input.ref_fa + '.fai') and os.path.isfile(input.ref_fa + '.gzi'):
            os.symlink(os.path.abspath(input.ref_fa + '.fai'), output.ref_fa + '.fai')
            os.symlink(os.path.abspath(input.ref_fa + '.gzi'), output.ref_fa + '.gzi')

        else:
            shell("""samtools faidx {output.ref_fa}""")

#
# Reference sequence partitions
#

# Create a table of merge batch assignments
localrules: call_merge_partition_table

# Create a table with reference sequences split evenly by length into a set number of partitions (labeled 0 to n - 1).
rule call_merge_partition_table:
    input:
        tsv='data/ref/contig_info.tsv.gz'
    output:
        tsv_part='data/ref/partition_{part_count}.tsv.gz'
    run:

        part_count = int(wildcards.part_count)

        if part_count < 1:
            raise RuntimeError(f'Number of partitions must be at least 1: {part_count}')

        # Read and sort
        df = pd.read_csv(
            input.tsv, sep='\t', dtype={'CHROM': str, 'LEN': int}
        ).sort_values(
            'LEN', ascending=False
        ).set_index(
            'CHROM'
        )[['LEN']]

        df['PARTITION'] = -1

        # Get a list of assignments for each partition
        list_chr = collections.defaultdict(list)
        list_size = collections.Counter()

        def get_smallest():
            """
            Get the next smallest bin.
            """

            min_index = 0

            for i in range(part_count):

                if list_size[i] == 0:
                    return i

                if list_size[i] < list_size[min_index]:
                    min_index = i

            return min_index

        for chrom in df.index:
            i = get_smallest()
            df.loc[chrom, 'PARTITION'] = i
            list_size[i] += df.loc[chrom, 'LEN']

        # Check
        if np.any(df['PARTITION'] < 0):
            raise RuntimeError('Failed to assign all reference contigs to partitions (PROGRAM BUG)')

        # Write
        df.to_csv(output.tsv_part, sep='\t', index=True, compression='gzip')
