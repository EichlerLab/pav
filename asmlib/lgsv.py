"""
Call large alignment-truncating variants.
"""

import collections
import numpy as np
import os
import pandas as pd
import sys

import asmlib.seq
import asmlib.inv
import kanapy.util.kmer

# Default parameters
MAX_TIG_DIST_PROP = 1  # Max allowed tig gap as a factor of the minimum alignment length of two records
MAX_REF_DIST_PROP = 3  # Max allowed ref gap as a factor of the minimum alignment length of two records

CALL_SOURCE = 'ALNTRUNC'  # Call source annotation for alignment-truncating events


def scan_for_events(df, df_tig_fai, hap, ref_fa_name, tig_fa_name, k_size, n_tree=None,
                    threads=1, log=sys.stdout, density_out_dir=None, max_tig_dist_prop=None, max_ref_dist_prop=None,
    ):
    """
    Scan trimmed alignments for alignment-truncating SV events.

    :param df: Dataframe of alignments (trimmed). Output from "align_cut_tig_overlap".
    :param df_tig_fai: Panadas Series keyed by contig names with values as contig lengths (integer).
    :param hap: Haplotype.
    :param ref_fa_name: Reference FASTA file name.
    :param tig_fa_name: Contig FASTA file name.
    :param k_size: K-mer size for inversion detection.
    :param n_tree: Locations of n-base regions in the reference to ignore regions with N's or `None` if no N-filtering
        should be attempted. If a region is expanded into N's, then it is discarded. If defined, this object should be
        dict-like keyed by chromosome names where each value is an IntervalTree of N bases.
    :param threads: Number of concurrent threads for inversion calling.
    :param log: Open file to write logs. Defaults to stdout.
    :param density_out_dir: Write density tables to this directory if not `None`. Directory must exist before this
        function is called.
    :param max_tig_dist_prop: Max allowed tig gap as a factor of the minimum alignment length of two records.
    :param max_ref_dist_prop: Max allowed ref gap as a factor of the minimum alignment length of two records.

    :return: A tuple of dataframes for SV calls: (INS, DEL, INV).
    """

    # Get parameters
    max_tig_dist_prop = max_tig_dist_prop if max_tig_dist_prop is not None else MAX_TIG_DIST_PROP
    max_ref_dist_prop = max_ref_dist_prop if max_ref_dist_prop is not None else MAX_REF_DIST_PROP

    # Copy df so it can be safely manipulated
    df = df.copy()

    df['QRY_LEN'] = df['END'] - df['POS']

    # Get liftover tool and k-mer util
    align_lift = asmlib.align.AlignLift(df, df_tig_fai)
    k_util = kanapy.util.kmer.KmerUtil(k_size)

    # with RefTigManager(ref_fa_name, tig_fa_name) as fa_pair:
    #
    #     ref_fa = fa_pair[0]
    #     tig_fa = fa_pair[0]

    # Setup lists
    ins_list = list()
    del_list = list()
    inv_list = list()

    # Get tig/chromosome combinations that appear more than once
    tig_map_count = collections.Counter(df[['#CHROM', 'QUERY_ID']].apply(tuple, axis=1))
    tig_map_count = [(chrom, tig_id) for (chrom, tig_id), count in tig_map_count.items() if count > 1]

    for chrom, tig_id in tig_map_count:

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

                if row2['REV'] == is_rev:
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

                            log.write('DEL: {}\n'.format(sv_id))
                            log.flush()

                            seq = asmlib.seq.region_seq_fasta(
                                asmlib.seq.Region(chrom, row1['END'], row2['POS']),
                                ref_fa_name
                            )

                            # Append
                            del_list.append(pd.Series(
                                [
                                    chrom, row1['END'], row2['POS'], sv_id, 'DEL', dist_ref,
                                    hap,
                                    tig_id, query_pos, query_pos + 1, '-' if row1['REV'] else '+',
                                    CALL_SOURCE, row1['CLUSTER_MATCH'],
                                    '-' if is_rev else '+',
                                    dist_tig,
                                    seq
                                ],
                                index=[
                                    '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
                                    'HAP',
                                    'QUERY_ID', 'QUERY_POS', 'QUERY_END', 'QUERY_STRAND',
                                    'CALL_SOURCE', 'CLUSTER_MATCH',
                                    'CI',
                                    'STRAND',
                                    'SEQ'
                                ]
                            ))

                            # Done with this contig break
                            break

                        elif dist_ref < 50:
                            # Call INS

                            seq = asmlib.seq.region_seq_fasta(
                                asmlib.seq.Region(tig_id, query_pos, query_end),
                                tig_fa_name,
                                rev_compl=is_rev
                            )

                            sv_id = '{}-{}-INS-{}'.format(chrom, row1['END'], dist_tig)

                            log.write('INS: {}\n'.format(sv_id))
                            log.flush()

                            ins_list.append(pd.Series(
                                [
                                    chrom, row1['END'], row1['END'] + 1, sv_id, 'INS', dist_tig,
                                    hap,
                                    tig_id,
                                    query_pos, query_end,
                                    '-' if row1['REV'] else '+',
                                    CALL_SOURCE, row1['CLUSTER_MATCH'],
                                    '-' if is_rev else '+',
                                    dist_ref,
                                    seq
                                ],
                                index=[
                                    '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
                                    'HAP',
                                    'QUERY_ID', 'QUERY_POS', 'QUERY_END', 'QUERY_STRAND',
                                    'CALL_SOURCE', 'CLUSTER_MATCH',
                                    'STRAND',
                                    'CI',
                                    'SEQ'
                                ]
                            ))

                            # Done with this contig break
                            break

                        else:

                            region_flag = asmlib.seq.Region(chrom, row1['END'], row2['POS'], is_rev=row1['REV'])

                            # Try INV
                            inv_call = asmlib.inv.scan_for_inv(
                                region_flag,
                                ref_fa_name, tig_fa_name,
                                align_lift, k_util,
                                threads=threads,
                                n_tree=n_tree,
                                log=log
                            )

                            if inv_call is not None:
                                log.write('INV (2-tig): {}\n'.format(inv_call))
                                log.flush()

                                seq = asmlib.seq.region_seq_fasta(
                                    inv_call.region_tig_outer,
                                    tig_fa_name,
                                    rev_compl=is_rev
                                )

                                inv_list.append(pd.Series(
                                    [
                                        inv_call.region_ref_outer.chrom,
                                        inv_call.region_ref_outer.pos,
                                        inv_call.region_ref_outer.end,
                                        inv_call.id,
                                        'INV',
                                        inv_call.svlen,

                                        inv_call.region_ref_outer.to_base1_string(),
                                        inv_call.region_ref_inner.to_base1_string(),

                                        inv_call.region_tig_outer.to_base1_string(),
                                        inv_call.region_tig_inner.to_base1_string(),

                                        inv_call.region_ref_discovery.to_base1_string(),
                                        inv_call.region_tig_discovery.to_base1_string(),
                                        inv_call.region_flag.to_base1_string(),
                                        inv_call.max_inv_den_diff,
                                        CALL_SOURCE,
                                        '-' if is_rev else '+',
                                        seq
                                    ],
                                    index=[
                                        '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
                                        'RGN_REF_OUTER', 'RGN_REF_INNER',
                                        'RGN_TIG_OUTER', 'RGN_TIG_INNER',
                                        'RGN_REF_DISC', 'RGN_TIG_DISC', 'RGN_REF_FLAG',
                                        'MAX_DENSITY_DIFF',
                                        'CALL_SOURCE',
                                        'STRAND',
                                        'SEQ'
                                    ]
                                ))

                                # Save density table
                                if density_out_dir is not None:
                                    inv_call.df.to_csv(
                                        os.path.join(density_out_dir,
                                                     'density_{}_{}.tsv.gz'.format(inv_call.id, hap)),
                                        sep='\t', index=False, compression='gzip'
                                    )

                                # Done with this contig break
                                break

                            # Else: Call SUB?

                    # Next record
                    subindex2 += 1

                elif subindex2 + 1 < tig_index_list_len:
                    # subindex1 and subindex2 are in opposite directions

                    subindex3 = subindex2 + 1

                    row3 = df.loc[tig_index_list[subindex3]]

                    if row3['REV'] == row1['REV']:  # Already known: row1['REV'] != row2['REV']

                        # Find inversion in 3-part alignment with middle in opposite orientation as flanks
                        region_flag = asmlib.seq.Region(chrom, row1['END'], row3['POS'], is_rev=row1['REV'])

                        # Try INV
                        inv_call = asmlib.inv.scan_for_inv(
                            region_flag,
                            ref_fa_name, tig_fa_name,
                            align_lift, k_util,
                            threads=threads,
                            n_tree=n_tree,
                            log=log
                        )

                        if inv_call is not None:
                            log.write('INV (3-tig): {}\n'.format(inv_call))
                            log.flush()

                            seq = asmlib.seq.region_seq_fasta(
                                inv_call.region_tig_outer,
                                tig_fa_name,
                                rev_compl=is_rev
                            )

                            inv_list.append(pd.Series(
                                [
                                    inv_call.region_ref_outer.chrom,
                                    inv_call.region_ref_outer.pos,
                                    inv_call.region_ref_outer.end,
                                    inv_call.id,
                                    'INV',
                                    inv_call.svlen,

                                    inv_call.region_ref_outer.to_base1_string(),
                                    inv_call.region_ref_inner.to_base1_string(),

                                    inv_call.region_tig_outer.to_base1_string(),
                                    inv_call.region_tig_inner.to_base1_string(),

                                    inv_call.region_ref_discovery.to_base1_string(),
                                    inv_call.region_tig_discovery.to_base1_string(),
                                    inv_call.region_flag.to_base1_string(),
                                    inv_call.max_inv_den_diff,
                                    CALL_SOURCE,
                                    '-' if is_rev else '+',
                                    seq
                                ],
                                index=[
                                    '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
                                    'RGN_REF_OUTER', 'RGN_REF_INNER',
                                    'RGN_TIG_OUTER', 'RGN_TIG_INNER',
                                    'RGN_REF_DISC', 'RGN_TIG_DISC', 'RGN_REF_FLAG',
                                    'MAX_DENSITY_DIFF',
                                    'CALL_SOURCE',
                                    'STRAND',
                                    'SEQ'
                                ]
                            ))

                            # Save density table
                            if density_out_dir is not None:
                                inv_call.df.to_csv(
                                    os.path.join(density_out_dir,
                                                 'density_{}_{}.tsv.gz'.format(inv_call.id, hap)),
                                    sep='\t', index=False, compression='gzip'
                                )

                            # Done with this contig break
                            break

                # Next subindex2
                subindex2 += 1

        # Advance
        subindex1 += 1

    # Concat records
    if len(ins_list) > 0:
        df_ins = pd.concat(ins_list, axis=1).T.sort_values(['#CHROM', 'POS', 'END'])
    else:
        df_ins = pd.DataFrame([], columns=[
            '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
            'HAP',
            'QUERY_ID', 'QUERY_POS', 'QUERY_END', 'QUERY_STRAND',
            'CALL_SOURCE', 'CLUSTER_MATCH',
            'STRAND',
            'CI',
            'SEQ'
        ])

    if len(del_list) > 0:
        df_del = pd.concat(del_list, axis=1).T.sort_values(['#CHROM', 'POS', 'END'])
    else:
        df_del = pd.DataFrame([], columns=[
            '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
            'HAP',
            'QUERY_ID', 'QUERY_POS', 'QUERY_END', 'QUERY_STRAND',
            'CALL_SOURCE', 'CLUSTER_MATCH',
            'STRAND',
            'CI',
            'SEQ'
        ])

    if len(inv_list) > 0:
        df_inv = pd.concat(inv_list, axis=1).T.sort_values(['#CHROM', 'POS', 'END'])
    else:
        df_inv = pd.DataFrame([], columns=[
            '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
            'RGN_REF_OUTER', 'RGN_REF_INNER',
            'RGN_TIG_OUTER', 'RGN_TIG_INNER',
            'RGN_REF_DISC', 'RGN_TIG_DISC', 'RGN_REF_FLAG',
            'MAX_DENSITY_DIFF',
            'CALL_SOURCE',
            'STRAND',
            'SEQ'
        ])

    # Return records
    return ((df_ins, df_del, df_inv))


# class RefTigManager:
#     """
#     Manage a pair of open FASTA files, reference and tig.
#     """
#
#     def __init__(self, ref_fa_name, tig_fa_name):
#
#         self.ref_fa = pysam.FastaFile(ref_fa_name)
#         self.tig_fa = pysam.FastaFile(tig_fa_name)
#
#     def __enter__(self):
#         return self.ref_fa, self.tig_fa
#
#     def __exit__(self, type, value, traceback):
#         self.ref_fa.close()
#         self.tig_fa.close()
