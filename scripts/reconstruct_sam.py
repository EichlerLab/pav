#!/usr/bin/env python3

"""
Reconstruct SAM file from alignment BAM, headers, and contig FASTA. PAV pipes the output into samtools to
create CRAM/BAM files from the alignment BAM files.
"""

import argparse
import gzip
import os
import pandas as pd
import pysam
import sys

import Bio.Seq

# Add PAV libraries and dependencies
PIPELINE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

sys.path.append(PIPELINE_DIR)  # pavlib
sys.path.append(os.path.join(PIPELINE_DIR, 'dep', 'svpop'))  # svpoplib
sys.path.append(os.path.join(PIPELINE_DIR, 'dep'))  # kanapy

import pavlib

if __name__ == '__main__':

    # Get command-line arguments
    parser = argparse.ArgumentParser('Reconstruct SAM from alignment BAM, headers, and contigs')

    parser.add_argument('--bed', help='PAV alignment BED file.')

    parser.add_argument('--headers', help='Original alignment headers.')

    parser.add_argument('--fasta', help='Unaligned contig FASTA.')

    parser.add_argument('--verbose', '-v',
                        default=False, action='store_true',
                        help='Output troubleshooting information to stderr'
                        )

    args = parser.parse_args()

    out_file = sys.stdout


    # Dump headers
    with gzip.open(args.headers, 'rt') as header_file:
        for line in header_file:
            out_file.write(line)

    # Read alignment records
    df = pd.read_csv(args.bed, sep='\t')

    # Open FASTA, write records
    with pysam.FastaFile(args.fasta) as fasta_file:
        for index, row in df.iterrows():

            if args.verbose:
                print('INDEX: {INDEX}'.format(**row), file=sys.stderr)
                sys.stderr.flush()

            # Get coordinates
            ref_bp, tig_bp, clip_h_l, clip_s_l, clip_h_r, clip_s_r = pavlib.align.count_cigar(row)

            if row['REV']:
                clip_r = clip_s_l
                clip_l = clip_s_r
            else:
                clip_r = clip_s_r
                clip_l = clip_s_l

            # Get sequence
            seq = fasta_file.fetch(row['QUERY_ID'], row['QUERY_TIG_POS'] - clip_l, row['QUERY_TIG_END'] + clip_r)

            # Reverse-complement
            if row['REV']:
                seq = str(Bio.Seq.Seq(seq).reverse_complement())

            sys.stdout.write(
                '{query_id:s}\t{flags:s}\t{chrom:s}\t{pos:d}\t{mapq:d}\t{cigar:s}\t{rnext:s}\t{pnext:d}\t{tlen:d}\t{seq:s}\t{qual:s}\n'.format(**{
                    'query_id': row['QUERY_ID'],
                    'flags': row['FLAGS'],
                    'chrom': row['#CHROM'],
                    'pos': row['POS'] + 1,
                    'mapq': row['MAPQ'],
                    'cigar': row['CIGAR'],
                    'rnext': '*',
                    'pnext': 0,
                    'tlen': 0,
                    'seq': seq,
                    'qual': '*'
                })
            )
