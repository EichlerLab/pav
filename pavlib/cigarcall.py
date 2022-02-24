"""
Call variants by CIGAR string.
"""

import Bio
import numpy as np
import pandas as pd
import pysam

import pavlib
import svpoplib


#
# Definitions
#

# Tag variants called with this source
CALL_SOURCE = 'CIGAR'

CALL_CIGAR_BATCH_COUNT = 10


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
                seq_ref = ref_fa.fetch(str(seq_ref_name))

        if seq_tig_name is None or row['QUERY_ID'] != seq_tig_name or is_rev != seq_tig_rev:
            with pysam.FastaFile(tig_fa_name) as tig_fa:
                seq_tig_name = row['QUERY_ID']
                seq_tig = tig_fa.fetch(str(seq_tig_name))
                seq_tig_len = len(seq_tig)

                if is_rev:
                    seq_tig = str(Bio.Seq.Seq(seq_tig).reverse_complement())

                seq_tig_rev = is_rev

        seq_ref_upper = seq_ref.upper()
        seq_tig_upper = seq_tig.upper()

        # Process CIGAR
        pos_ref = row['POS']
        pos_tig = 0

        cigar_index = 0

        last_op = None
        last_oplen = 0

        for oplen, op in pavlib.align.cigar_str_to_tuples(row):
            # NOTE: break/continue in this look will not advance last_op and last_oplen (end of loop)

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
                    var_id = f'{seq_ref_name}-{pos_ref_snv + 1}-SNV-{base_ref.upper()}-{base_tig.upper()}'

                    df_snv_list.append(pd.Series(
                        [
                            seq_ref_name, pos_ref_snv, pos_ref_snv + 1,
                            var_id, 'SNV', 1,
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

                # Get sequence
                seq = seq_tig[pos_tig:(pos_tig + oplen)]
                seq_upper = seq.upper()

                # Left shift
                if last_op == '=':
                    left_shift = np.min([
                        last_oplen,
                        pavlib.call.left_homology(pos_ref - 1, seq_ref_upper, seq_upper)  # SV/breakpoint upstream homology
                    ])
                else:
                    left_shift = 0

                sv_pos_ref = pos_ref - left_shift
                sv_end_ref = sv_pos_ref + 1
                sv_pos_tig = pos_tig - left_shift
                sv_end_tig = sv_pos_tig + oplen

                if left_shift != 0:
                    seq = seq_tig[sv_pos_tig:(sv_pos_tig + oplen)]

                # Get positions in the original SV space
                # pos_tig_insdel to fwd contig if alignment is reversed
                if is_rev:
                    end_tig_insdel = seq_tig_len - sv_pos_tig
                    pos_tig_insdel = end_tig_insdel - oplen

                else:
                    pos_tig_insdel = sv_pos_tig
                    end_tig_insdel = pos_tig_insdel + oplen

                # Find breakpoint homology
                seq_upper = seq.upper()

                hom_ref_l = pavlib.call.left_homology(sv_pos_ref - 1, seq_ref_upper, seq_upper)
                hom_ref_r = pavlib.call.right_homology(sv_pos_ref, seq_ref_upper, seq_upper)

                hom_tig_l = pavlib.call.left_homology(sv_pos_tig - 1, seq_tig_upper, seq_upper)
                hom_tig_r = pavlib.call.right_homology(sv_end_tig, seq_tig_upper, seq_upper)

                # Add variant
                var_id = f'{seq_ref_name}-{sv_pos_ref + 1}-INS-{oplen}'

                df_insdel_list.append(pd.Series(
                    [
                        seq_ref_name, sv_pos_ref, sv_end_ref,
                        var_id, 'INS', oplen,
                        hap,
                        f'{seq_tig_name}:{pos_tig_insdel + 1}-{end_tig_insdel}', strand,
                        0,
                        align_index, cluster_match,
                        left_shift, f'{hom_ref_l},{hom_ref_r}', f'{hom_tig_l},{hom_tig_r}',
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
                        'LEFT_SHIFT', 'HOM_REF', 'HOM_TIG',
                        'CALL_SOURCE',
                        'SEQ'
                    ]
                ))

                # Advance
                pos_tig += oplen

                pass

            elif op == 'D':
                # Call DEL

                # Get sequence
                seq = seq_ref[pos_ref:(pos_ref + oplen)]
                seq_upper = seq.upper()

                # Left shift
                if last_op == '=':
                    left_shift = np.min([
                        last_oplen,
                        pavlib.call.left_homology(pos_ref - 1, seq_ref_upper, seq_upper)  # SV/breakpoint upstream homology
                    ])
                else:
                    left_shift = 0

                sv_pos_ref = pos_ref - left_shift
                sv_end_ref = sv_pos_ref + oplen
                sv_pos_tig = pos_tig - left_shift
                sv_end_tig = sv_pos_tig + 1

                # Contig position in original coordinates (translate if - strand)
                pos_tig_insdel = sv_pos_tig

                if is_rev:
                    pos_tig_insdel = seq_tig_len - sv_pos_tig

                # Find breakpoint homology
                seq_upper = seq.upper()

                hom_ref_l = pavlib.call.left_homology(sv_pos_ref - 1, seq_ref_upper, seq_upper)
                hom_ref_r = pavlib.call.right_homology(sv_end_ref, seq_ref_upper, seq_upper)

                hom_tig_l = pavlib.call.left_homology(sv_pos_tig - 1, seq_tig_upper, seq_upper)
                hom_tig_r = pavlib.call.right_homology(sv_pos_tig, seq_tig_upper, seq_upper)

                # Add variant
                var_id = f'{seq_ref_name}-{pos_ref + 1}-DEL-{oplen}'

                df_insdel_list.append(pd.Series(
                    [
                        seq_ref_name, pos_ref, pos_ref + oplen,
                        var_id, 'DEL', oplen,
                        hap,
                        f'{seq_tig_name}:{pos_tig_insdel + 1}-{pos_tig_insdel + 1}', strand,
                        0,
                        align_index, cluster_match,
                        left_shift, f'{hom_ref_l},{hom_ref_r}', f'{hom_tig_l},{hom_tig_r}',
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
                        'LEFT_SHIFT', 'HOM_REF', 'HOM_TIG',
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
                        'opcode={}, subject={}:{} , query={}:{}, align-index={}'
                    ).format(
                        cigar_index, op, seq_ref_name, pos_ref, seq_tig_name, pos_tig, row['INDEX']
                    ))

            # Save last op
            last_op = op
            last_oplen = oplen

    # Merge tables
    if len(df_snv_list) > 0:
        df_snv = pd.concat(df_snv_list, axis=1).T
        df_snv['ID'] = svpoplib.variant.version_id(df_snv['ID'])
        df_snv.sort_values(['#CHROM', 'POS', 'END', 'ID'], inplace=True)

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
        df_insdel = pd.concat(df_insdel_list, axis=1).T
        df_insdel['ID'] = svpoplib.variant.version_id(df_insdel['ID'])
        df_insdel.sort_values(['#CHROM', 'POS', 'END', 'ID'], inplace=True)

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
                'LEFT_SHIFT', 'HOM_REF', 'HOM_TIG',
                'CALL_SOURCE',
                'SEQ'
            ]
        )

    # Return tables
    return df_snv, df_insdel
