"""
Data files including reference and data tables for the reference.
"""

import gzip
import os

import Bio.SeqIO

import pavlib
import svpoplib

global shell


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
