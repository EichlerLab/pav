"""
Call large alignment-truncating variants.
"""

import collections
import numpy as np
import pandas as pd

import Bio


# Default parameters
MAX_TIG_DIST_PROP = 1  # Max allowed tig gap as a factor of the minimum alignment length of two records
MAX_REF_DIST_PROP = 3  # Max allowed ref gap as a factor of the minimum alignment length of two records

CALL_SOURCE = 'ALNTRUNC'  # Call source annotation for alignment-truncating events


def scan_for_events(df, df_tig_fai, hap, ref_fa_name, tig_fa_name, max_tig_dist_prop=None, max_ref_dist_prop=None):
    """
    Scan trimmed alignments for alignment-truncating SV events.

    :param df: Dataframe of alignments (trimmed). Output from "align_cut_tig_overlap".
    :param df_tig_fai: Panadas Series keyed by contig names with values as contig lengths (integer).
    :param hap: Haplotype.
    :param ref_fa_name: Reference FASTA file name.
    :param tig_fa_name: Contig FASTA file name.
    :param max_tig_dist_prop: Max allowed tig gap as a factor of the minimum alignment length of two records.
    :param max_ref_dist_prop: Max allowed ref gap as a factor of the minimum alignment length of two records.

    :return:
    """

    # Get parameters
    max_tig_dist_prop = max_tig_dist_prop if max_tig_dist_prop is not None else MAX_TIG_DIST_PROP
    max_ref_dist_prop = max_ref_dist_prop if max_ref_dist_prop is not None else MAX_REF_DIST_PROP

    # Copy df so it can be safely manipulated
    df = df.copy()

    df['QRY_LEN'] = df['END'] - df['POS']

    # Get liftover tool and contig file
    align_lift = asmlib.align.AlignLift(df, df_tig_fai)

    with RefTigManager(ref_fa_name, tig_fa_name) as fa_pair:

        ref_fa = fa_pair[0]
        tig_fa = fa_pair[1]

        # Setup lists
        ins_list = list()
        del_list = list()
        inv_list = list()
        sub_list = list()

        # Get tig/chromosome combinations that appear more than once
        tig_map_count = collections.Counter(df[['#CHROM', 'QUERY_ID']].apply(tuple, axis=1))
        tig_map_count = [(chrom, tig_id) for (chrom, tig_id), count in tig_map_count.items() if count > 1]

        for chrom, tig_id in [(chrom, tig_id) for (chrom, tig_id), count in tig_map_count.items() if count > 1]:

            # Get a list of df indices with alignments for this contig
            tig_index_list = list(df.loc[(df['#CHROM'] == chrom) & (df['QUERY_ID'] == tig_id)].index)
            tig_index_list_len = len(tig_index_list)

            # subindex1 is an index to tig_index_list. Traverse from first element to second from last
            # searching for alignment-truncating SVs
            # * tig_index_list[subindex1] is the index to df
            for subindex1 in range(len(tig_index_list) - 1):
                subindex2 = subindex1 + 1

                row1 = df.loc[tig_index_list[subindex1]]

                is_rev = row1['REV']

                while subindex2 < tig_index_list_len:

                    row2 = df.loc[tig_index_list[subindex2]]

                    # Report
                    #            print('Break: {}:{:,d}/{:,d} (dtig={}, dref={}) ({}, index={},{} (sub={},{}))'.format(
                    #                row1['#CHROM'], row1['END'], row2['POS'],
                    #                dist_tig, dist_ref,
                    #                row1['QUERY_ID'],
                    #                tig_index_list[subindex1],
                    #                tig_index_list[subindex2],
                    #                subindex1, subindex2
                    #            ))

                    if row1['REV'] == row2['REV']:
                        # Scan for INS/DEL

                        query_pos = row2['QUERY_TIG_END'] if is_rev else row1['QUERY_TIG_END']
                        query_end = row1['QUERY_TIG_POS'] if is_rev else row2['QUERY_TIG_POS']

                        dist_tig = query_end - query_pos
                        dist_ref = row2['POS'] - row1['END']

                        min_aln_len = np.min([row1['QRY_LEN'], row2['QRY_LEN']])

                        # Stop if there is large gap between the aligned bases (tig or ref)
                        if (
                                np.abs(dist_tig) / min_aln_len > max_tig_dist_prop
                        ) or (
                                np.abs(dist_ref) / min_aln_len > max_ref_dist_prop
                        ):
                            subindex2 += 1
                            continue

                        if dist_tig >= 50 or dist_ref >= 50:

                            if dist_tig < 50:
                                # Call DEL

                                sv_id = '{}-{}-DEL-{}'.format(chrom, row1['END'], dist_ref)

                                seq = ref_fa.fetch(chrom, row1['END'], row2['POS'])

                                # Append
                                del_list.append(pd.Series(
                                    [
                                        chrom, row1['END'], row2['POS'], sv_id, 'DEL', dist_ref,
                                        hap,
                                        tig_id, query_pos, query_pos + 1, '-' if row1['REV'] else '+',
                                        CALL_SOURCE, row1['CLUSTER_MATCH'],
                                        seq
                                    ],
                                    index=[
                                        '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
                                        'HAP',
                                        'QUERY_ID', 'QUERY_POS', 'QUERY_END', 'QUERY_STRAND',
                                        'CALL_SOURCE', 'CLUSTER_MATCH',
                                        'SEQ'
                                    ]
                                ))

                                # Done with this contig break
                                break

                            elif dist_ref < 50:
                                # Call INS

                                seq = tig_fa.fetch(row1['QUERY_ID'], query_pos, query_end)

                                sv_id = '{}-{}-INS-{}'.format(chrom, row1['END'], dist_tig)

                                if is_rev:
                                    seq = str(Bio.Seq.Seq(seq).reverse_complement())

                                ins_list.append(pd.Series(
                                    [
                                        chrom, row1['END'], row1['END'] + 1, sv_id, 'INS', dist_tig,
                                        hap,
                                        tig_id,
                                        query_pos, query_end,
                                        '-' if row1['REV'] else '+',
                                        CALL_SOURCE,
                                        seq
                                    ],
                                    index=[
                                        '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
                                        'HAP',
                                        'QUERY_ID',
                                        'QUERY_POS', 'QUERY_END',
                                        'QUERY_STRAND',
                                        'CALL_SOURCE',
                                        'SEQ'
                                    ]
                                ))

                                # Done with this contig break
                                break

                            else:
                                # Call SUB
                                pass

                        # Next record
                        subindex2 += 1

                    else:
                        # subindex1 and subindex3 are in opposite directions

                        subindex3 = subindex2 + 1

                        row3 = df.loc[tig_index_list[subindex3]]

                        if row3['REV'] == row1['REV']:
                            # Find inversion

                            if INVERSION-FOUND:
                                # Call inversion


                                # Save record
                                inv_list.append(pd.Series(
                                    [
                                        chrom, row1['END'], row1['END'] + 1, sv_id, 'INS', dist_tig,
                                        hap,
                                        tig_id,
                                        query_pos, query_end,
                                        '-' if row1['REV'] else '+',
                                        CALL_SOURCE,
                                        seq
                                    ],
                                    index=[
                                        '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
                                        'HAP',
                                        'QUERY_ID',
                                        'QUERY_POS', 'QUERY_END',
                                        'QUERY_STRAND',
                                        'CALL_SOURCE',
                                        'SEQ'
                                    ]
                                ))

                                # Consume breakpoints
                                subindex1 = subindex3
                                break

                    # Next subindex2
                    subindex2 += 1

            # Advance
            subindex1 += 1

    # Return
    return df_ins, df_del, df_inv, df_sub


class RefTigManager:
    """
    Manage a pair of open FASTA files, reference and tig.
    """

    def __init__(self, ref_fa_name, tig_fa_name):

        self.ref_fa = pysam.FastaFile(ref_fa_name)
        self.tig_fa = pysam.FastaFile(tig_fa_name)

    def __enter__(self):
        return self.ref_fa, self.tig_fa

    def __exit__(self, type, value, traceback):
        self.ref_fa.close()
        self.tig_fa.close()
