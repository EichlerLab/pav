"""
Call variants by CIGAR string.
"""

import Bio
import pandas as pd
import pysam

import pavlib


#
# Definitions
#

# Tag variants called with this source
CALL_SOURCE = 'CIGAR'


def make_insdel_snv_calls(df_align, ref_fa_name, tig_fa_name, hap):
    """
    Parse variants from CIGAR strings.

    :param df_align: Post-cut BED of read alignments.
    :param ref_fa_name: Reference FASTA file name.
    :param tig_fa_name: Contig FASTA file name.
    :param hap: String identifying the haplotype ("h1", "h2").

    :return: A tuple of two dataframes, one for insertions and deletions (SV and indel), and one for SNVs.
    """

    df_insdel_list = list()
    df_snv_list = list()

    seq_ref = None       # Current reference sequence
    seq_ref_name = None  # Current reference contig name

    seq_tig = None       # Current tig sequence
    seq_tig_name = None  # Current tig name
    seq_tig_len = None   # Current tig length
    seq_tig_rev = None   # Aligned contig was reverse-complemented

    # Parse alignment records
    for index, row in df_align.iterrows():

        # Get strand
        is_rev = row['REV']
        strand = '-' if is_rev else '+'
        cluster_match = row['CLUSTER_MATCH']
        align_index = row['INDEX']

        # Load reference and tig sequences
        if seq_ref_name is None or row['#CHROM'] != seq_ref_name:
            with pysam.FastaFile(ref_fa_name) as ref_fa:
                seq_ref_name = row['#CHROM']
                seq_ref = ref_fa.fetch(str(seq_ref_name)).upper()

        if seq_tig_name is None or row['QUERY_ID'] != seq_tig_name or is_rev != seq_tig_rev:
            with pysam.FastaFile(tig_fa_name) as tig_fa:
                seq_tig_name = row['QUERY_ID']
                seq_tig = tig_fa.fetch(str(seq_tig_name)).upper()
                seq_tig_len = len(seq_tig)

                if is_rev:
                    seq_tig = str(Bio.Seq.Seq(seq_tig).reverse_complement())

                seq_tig_rev = is_rev

        # Process CIGAR
        pos_ref = row['POS']
        pos_tig = 0

        cigar_index = 0

        for oplen, op in pavlib.align.cigar_str_to_tuples(row):

            cigar_index += 1

            if op == '=':
                pos_ref += oplen
                pos_tig += oplen

            elif op == 'X':
                # Call SNV(s)

                for i in range(oplen):

                    # Get position and bases
                    pos_ref_snv = pos_ref + i
                    pos_tig_snv = pos_tig + i

                    base_ref = seq_ref[pos_ref_snv]
                    base_tig = seq_tig[pos_tig_snv]

                    # pos_tig_snv to fwd contig if alignment is reversed
                    if is_rev:
                        pos_tig_snv = seq_tig_len - pos_tig_snv - 1

                    # Add variant
                    df_snv_list.append(pd.Series(
                        [
                            seq_ref_name, pos_ref_snv, pos_ref_snv + 1,
                            f'{seq_ref_name}-{pos_ref_snv + 1}-SNV-{base_ref}-{base_tig}', 'SNV', 1,
                            base_ref, base_tig,
                            hap,
                            f'{seq_tig_name}:{pos_tig_snv + 1}-{pos_tig_snv + 1}', strand,
                            0,
                            align_index, cluster_match,
                            CALL_SOURCE
                        ],
                        index=[
                            '#CHROM', 'POS', 'END',
                            'ID', 'SVTYPE', 'SVLEN',
                            'REF', 'ALT',
                            'HAP',
                            'TIG_REGION', 'QUERY_STRAND',
                            'CI',
                            'ALIGN_INDEX', 'CLUSTER_MATCH',
                            'CALL_SOURCE'
                        ]
                    ))

                # Advance
                pos_ref += oplen
                pos_tig += oplen

            elif op == 'I':
                # Call INS

                pos_tig_insdel = pos_tig

                # Get sequence
                seq = seq_tig[pos_tig:(pos_tig + oplen)]

                # pos_tig_insdel to fwd contig if alignment is reversed
                if is_rev:
                    end_tig_insdel = seq_tig_len - pos_tig_insdel
                    pos_tig_insdel = end_tig_insdel - oplen

                else:
                    end_tig_insdel = pos_tig_insdel + oplen

                # Add variant
                df_insdel_list.append(pd.Series(
                    [
                        seq_ref_name, pos_ref, pos_ref + 1,
                        f'{seq_ref_name}-{pos_ref + 1}-INS-{oplen}', 'INS', oplen,
                        hap,
                        f'{seq_tig_name}:{pos_tig_insdel + 1}-{end_tig_insdel}', strand,
                        0,
                        align_index, cluster_match,
                        CALL_SOURCE,
                        seq
                    ],
                    index=[
                        '#CHROM', 'POS', 'END',
                        'ID', 'SVTYPE', 'SVLEN',
                        'HAP',
                        'TIG_REGION', 'QUERY_STRAND',
                        'CI',
                        'ALIGN_INDEX', 'CLUSTER_MATCH',
                        'CALL_SOURCE',
                        'SEQ'
                    ]
                ))

                # Advance
                pos_tig += oplen

                pass

            elif op == 'D':
                # Call DEL

                pos_tig_insdel = pos_tig

                # Get sequence
                seq = seq_ref[pos_ref:(pos_ref + oplen)]

                # pos_tig_insdel to fwd contig if alignment is reversed
                if is_rev:
                    pos_tig_insdel = seq_tig_len - pos_tig_insdel

                # Add variant
                df_insdel_list.append(pd.Series(
                    [
                        seq_ref_name, pos_ref, pos_ref + oplen,
                        f'{seq_ref_name}-{pos_ref + 1}-DEL-{oplen}', 'DEL', oplen,
                        hap,
                        f'{seq_tig_name}:{pos_tig_insdel + 1}-{pos_tig_insdel + 1}', strand,
                        0,
                        align_index, cluster_match,
                        CALL_SOURCE,
                        seq
                    ],
                    index=[
                        '#CHROM', 'POS', 'END',
                        'ID', 'SVTYPE', 'SVLEN',
                        'HAP',
                        'TIG_REGION', 'QUERY_STRAND',
                        'CI',
                        'ALIGN_INDEX', 'CLUSTER_MATCH',
                        'CALL_SOURCE',
                        'SEQ'
                    ]
                ))

                # Advance
                pos_ref += oplen

                pass

            elif op in {'S', 'H'}:
                pos_tig += oplen

            else:
                # Cannot handle CIGAR operation

                if op == 'M':
                    raise RuntimeError((
                        'Illegal operation code in CIGAR string at operation {}: '
                        'Alignments must be generated with =/X (not M): '
                        'opcode={}, subject={}:{}, query={}:{}, align-index={}'
                    ).format(
                        cigar_index, op, seq_ref_name, pos_ref, seq_tig_name, pos_tig, row['INDEX']
                    ))

                else:
                    raise RuntimeError((
                        'Illegal operation code in CIGAR string at operation {}: '
                        'opcode={}, subject={}:{}, query={}:{}, align-index={}'
                    ).format(
                        cigar_index, op, seq_ref_name, pos_ref, seq_tig_name, pos_tig, row['INDEX']
                    ))

    # Merge tables
    if len(df_snv_list) > 0:
        df_snv = pd.concat(df_snv_list, axis=1).T.sort_values(['#CHROM', 'POS'])
    else:
        df_snv = pd.DataFrame(
            [],
            columns=[
                '#CHROM', 'POS', 'END',
                'ID', 'SVTYPE', 'SVLEN',
                'REF', 'ALT',
                'HAP',
                'TIG_REGION', 'QUERY_STRAND',
                'CI',
                'ALIGN_INDEX', 'CLUSTER_MATCH',
                'CALL_SOURCE'
            ]
        )

    if len(df_insdel_list) > 0:
        df_insdel = pd.concat(df_insdel_list, axis=1).T.sort_values(['#CHROM', 'POS'])
    else:
        df_insdel = pd.DataFrame(
            [],
            columns=[
                '#CHROM', 'POS', 'END',
                'ID', 'SVTYPE', 'SVLEN',
                'HAP',
                'TIG_REGION', 'QUERY_STRAND',
                'CI',
                'ALIGN_INDEX', 'CLUSTER_MATCH',
                'CALL_SOURCE',
                'SEQ'
            ]
        )

    # Return tables
    return df_snv, df_insdel
