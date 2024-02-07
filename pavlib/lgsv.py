"""
Call large alignment-truncating variants.
"""

import collections
import numpy as np
import os
import pandas as pd
import sys

import pavlib.seq
import pavlib.inv
import kanapy.util.kmer
import svpoplib

import Bio.Seq

# Default parameters
MAX_TIG_DIST_PROP = 1  # Max allowed tig gap as a factor of the minimum alignment length of two records
MAX_REF_DIST_PROP = 3  # Max allowed ref gap as a factor of the minimum alignment length of two records

DIST_PROP_LEN_MAPQ = (20000, 40)  # If a flanking alignment is at least as long as the first element and MAPQ is at
# least the second element, ignore MAX_TIG_DIST_PROP and MAX_REF_DIST_PROP and accept the alignments.

CALL_SOURCE = 'ALNTRUNC'  # Call source annotation for alignment-truncating events

CALL_SOURCE_INV_DENSITY = 'ALNTRUNC-DEN'
CALL_SOURCE_INV_NO_DENSITY = 'ALNTRUNC-NODEN'


def scan_for_events(df, df_tig_fai, hap, ref_fa_name, tig_fa_name, k_size, n_tree=None,
                    threads=1, log=sys.stdout, density_out_dir=None, max_tig_dist_prop=None, max_ref_dist_prop=None,
                    srs_tree=None, max_region_size=None
    ):
    """
    Scan trimmed alignments for alignment-truncating SV events.

    :param df: Dataframe of alignments (trimmed). Output from "align_cut_tig_overlap".
    :param df_tig_fai: Pandas Series keyed by contig names with values as contig lengths (integer).
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
    :param srs_tree: Inversion density "--staterunsmooth" parameters (for `density.py`). May be a tree, a list of
        limits, or `None` to use the default for all sizes. See `inv.scan_for_inv()` for details.
    :param max_region_size: Max region size for inversion scanning. Value 0 disables the limit, and `None` sets the
        default limit, `inv.MAX_REGION_SIZE`. See `inv.scan_for_inv()` for details.

    :return: A tuple of dataframes for SV calls: (INS, DEL, INV).
    """

    # Get parameters
    max_tig_dist_prop = max_tig_dist_prop if max_tig_dist_prop is not None else MAX_TIG_DIST_PROP
    max_ref_dist_prop = max_ref_dist_prop if max_ref_dist_prop is not None else MAX_REF_DIST_PROP

    # Copy df so it can be safely manipulated
    df = df.copy()

    df['QRY_LEN'] = df['END'] - df['POS']

    # Get liftover tool and k-mer util
    align_lift = pavlib.align.AlignLift(df, df_tig_fai)
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

    # Cache reference sequence (upper-case for homology searches)
    seq_cache_ref = SeqCache(ref_fa_name, uppercase=True)
    seq_cache_tig = SeqCache(tig_fa_name, uppercase=True)

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

                    # Determine if contigs are in reference orientation or reverse orientation
                    # In reference orientation, the left-most alignment record on the reference (row1) is also
                    # left-most in the contig sequence. If the contig is reversed, then row1 comes before row2 on the
                    # reference, but row1 comes after row2 on the contig.
                    if row1['QUERY_TIG_POS'] < row2['QUERY_TIG_POS']:
                        if row2['QUERY_TIG_POS'] < row1['QUERY_TIG_END']:
                            raise RuntimeError(
                                'Contig ranges overlap for two alignment records (should not occur after alignment trimming): '
                                f'Index {row1["INDEX"]} ({row1["QUERY_ID"]}:{row1["QUERY_TIG_POS"]}-{row1["QUERY_TIG_END"]}) and '
                                f'Index {row2["INDEX"]} ({row2["QUERY_ID"]}:{row2["QUERY_TIG_POS"]}-{row2["QUERY_TIG_END"]})'
                            )

                        query_pos = row1['QUERY_TIG_END']
                        query_end = row2['QUERY_TIG_POS']

                    else:
                        if row1['QUERY_TIG_POS'] < row2['QUERY_TIG_END']:
                            raise RuntimeError(
                                'Contig ranges overlap for two alignment records (should not occur after alignment trimming): '
                                f'Index {row1["INDEX"]} ({row1["QUERY_ID"]}:{row1["QUERY_TIG_POS"]}-{row1["QUERY_TIG_END"]}) and '
                                f'Index {row2["INDEX"]} ({row2["QUERY_ID"]}:{row2["QUERY_TIG_POS"]}-{row2["QUERY_TIG_END"]})'
                            )

                        query_pos = row2['QUERY_TIG_END']
                        query_end = row1['QUERY_TIG_POS']

                    dist_tig = query_end - query_pos
                    dist_ref = row2['POS'] - row1['END']

                    if dist_tig < 0:
                        raise RuntimeError(
                            f'Contig query positions are out of order (program bug): Contig distance is negative ({dist_tig}): '
                            f'Index {row1["INDEX"]} ({row1["QUERY_ID"]}:{row1["QUERY_TIG_POS"]}-{row1["QUERY_TIG_END"]}) and '
                            f'Index {row2["INDEX"]} ({row2["QUERY_ID"]}:{row2["QUERY_TIG_POS"]}-{row2["QUERY_TIG_END"]})'
                        )

                    min_aln_len = np.min([row1['QRY_LEN'], row2['QRY_LEN']])
                    min_mapq = np.min([row1['MAPQ'], row2['MAPQ']])

                    # Stop at large gap between the aligned bases where flanking alignments are not well supported
                    if min_aln_len < DIST_PROP_LEN_MAPQ[0] or min_mapq < DIST_PROP_LEN_MAPQ[1]:

                        # Rescue if alignments are sufficently large
                        if (
                                np.abs(dist_tig) / min_aln_len > max_tig_dist_prop
                        ) or (
                                np.abs(dist_ref) / min_aln_len > max_ref_dist_prop
                        ):
                            subindex2 += 1
                            continue

                    # Call DEL
                    if dist_ref >= 50 and dist_tig < 50:

                        # Call DEL
                        svlen = dist_ref

                        pos_ref = row1['END']
                        end_ref = row2['POS']
                        pos_tig = query_pos
                        end_tig = pos_tig + 1

                        ref_region = pavlib.seq.Region(chrom, pos_ref, end_ref)

                        # Get reference and tig sequences (left shift and homology search)
                        seq_ref = seq_cache_ref.get(chrom, False)
                        seq_tig = seq_cache_tig.get(tig_id, is_rev)

                        # Get SV sequence
                        seq = pavlib.seq.region_seq_fasta(ref_region, ref_fa_name)

                        # Left-shift through matching bases and get position
                        left_shift = np.min([
                            pavlib.align.match_bp(row1, True),  # Do not shift through anything but matched bases
                            pavlib.call.left_homology(pos_ref - 1, seq_ref, seq.upper())  # SV/breakpoint upstream homology
                        ])

                        if left_shift > 0:
                            pos_ref -= left_shift
                            end_ref -= left_shift
                            pos_tig -= left_shift
                            end_tig -= left_shift

                            ref_region = pavlib.seq.Region(chrom, pos_ref, end_ref)
                            seq = pavlib.seq.region_seq_fasta(ref_region, ref_fa_name)

                        # Get ID and log
                        sv_id = '{}-{}-DEL-{}'.format(chrom, pos_ref, svlen)

                        log.write('DEL: {}\n'.format(sv_id))
                        log.flush()

                        # Find breakpoint homology
                        seq_upper = seq.upper()

                        hom_ref_l = pavlib.call.left_homology(pos_ref - 1, seq_ref, seq_upper)
                        hom_ref_r = pavlib.call.right_homology(end_ref, seq_ref, seq_upper)

                        hom_tig_l = pavlib.call.left_homology(pos_tig - 1, seq_tig, seq_upper)
                        hom_tig_r = pavlib.call.right_homology(pos_tig, seq_tig, seq_upper)

                        # Append
                        del_list.append(pd.Series(
                            [
                                chrom, pos_ref, end_ref,
                                sv_id, 'DEL', svlen,
                                hap,
                                f'{tig_id}:{pos_tig + 1}-{end_tig}', '-' if row1['REV'] else '+',
                                dist_tig,
                                '{},{}'.format(row1['INDEX'], row2['INDEX']), row1['CLUSTER_MATCH'],
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

                        # Done with this contig break
                        break

                    # Call INS
                    elif dist_ref < 50 and dist_tig >= 50:

                        # Call INS
                        pos_ref = row1['END']
                        end_ref = pos_ref + 1
                        pos_tig = query_pos
                        end_tig = query_end

                        svlen = dist_tig

                        tig_region = pavlib.seq.Region(tig_id, pos_tig, end_tig, is_rev=is_rev)

                        # Get reference and tig sequences (left shift and homology search)
                        seq_ref = seq_cache_ref.get(chrom, False)
                        seq_tig = seq_cache_tig.get(tig_id, is_rev)

                        # SV sequence
                        seq = pavlib.seq.region_seq_fasta(
                            tig_region,
                            tig_fa_name,
                            rev_compl=is_rev
                        )

                        # Left-shift through matching bases and get position
                        left_shift = np.min([
                            pavlib.align.match_bp(row1, True),  # Do not shift through anything but matched bases
                            pavlib.call.left_homology(pos_ref - 1, seq_ref, seq.upper())  # SV/breakpoint upstream homology
                        ])

                        if left_shift > 0:
                            pos_ref -= left_shift
                            end_ref -= left_shift
                            pos_tig -= left_shift
                            end_tig -= left_shift

                            tig_region = pavlib.seq.Region(tig_id, pos_tig, end_tig, is_rev=is_rev)

                            seq = pavlib.seq.region_seq_fasta(
                                tig_region,
                                tig_fa_name,
                                rev_compl=is_rev
                            )

                        # Get ID and log
                        sv_id = '{}-{}-INS-{}'.format(chrom, pos_ref, svlen)

                        log.write('INS: {}\n'.format(sv_id))
                        log.flush()

                        # Find breakpoint homology
                        seq_upper = seq.upper()

                        hom_ref_l = pavlib.call.left_homology(pos_ref - 1, seq_ref, seq_upper)
                        hom_ref_r = pavlib.call.right_homology(pos_ref, seq_ref, seq_upper)

                        hom_tig_l = pavlib.call.left_homology(pos_tig - 1, seq_tig, seq_upper)
                        hom_tig_r = pavlib.call.right_homology(end_tig, seq_tig, seq_upper)

                        # APPEND
                        ins_list.append(pd.Series(
                            [
                                chrom, pos_ref, end_ref,
                                sv_id, 'INS', svlen,
                                hap,
                                tig_region.to_base1_string(), '-' if is_rev else '+',
                                dist_ref,
                                '{},{}'.format(row1['INDEX'], row2['INDEX']), row1['CLUSTER_MATCH'],
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

                        # Done with this contig break
                        break

                    # INV (2 contig)
                    elif dist_ref >= 50 and dist_tig >= 50:

                        region_flag = pavlib.seq.Region(chrom, row1['END'], row2['POS'], is_rev=row1['REV'])

                        # Try INV
                        inv_call = pavlib.inv.scan_for_inv(
                            region_flag,
                            ref_fa_name, tig_fa_name,
                            align_lift, k_util,
                            max_region_size=max_region_size,
                            threads=threads,
                            n_tree=n_tree,
                            srs_tree=srs_tree,
                            log=log,
                            min_exp_count=1  # If alignment truncation does not contain inverted k-mers at sufficient density, stop searching
                        )

                        if inv_call is not None:
                            log.write('INV (2-tig): {}\n'.format(inv_call))
                            log.flush()

                            seq = pavlib.seq.region_seq_fasta(
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

                                    hap,

                                    inv_call.region_tig_outer.to_base1_string(),
                                    '-' if is_rev else '+',

                                    0,

                                    inv_call.region_ref_inner.to_base1_string(),
                                    inv_call.region_tig_inner.to_base1_string(),

                                    inv_call.region_ref_discovery.to_base1_string(),
                                    inv_call.region_tig_discovery.to_base1_string(),

                                    inv_call.region_flag.region_id(),
                                    'ALNTRUNC',

                                    '{},{}'.format(row1['INDEX'], row2['INDEX']),
                                    row1['CLUSTER_MATCH'],

                                    CALL_SOURCE_INV_DENSITY,

                                    seq
                                ],
                                index=[
                                    '#CHROM', 'POS', 'END',
                                    'ID', 'SVTYPE', 'SVLEN',
                                    'HAP',
                                    'TIG_REGION', 'QUERY_STRAND',
                                    'CI',
                                    'RGN_REF_INNER', 'RGN_TIG_INNER',
                                    'RGN_REF_DISC', 'RGN_TIG_DISC',
                                    'FLAG_ID', 'FLAG_TYPE',
                                    'ALIGN_INDEX', 'CLUSTER_MATCH',
                                    'CALL_SOURCE',
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

                # INV (3 contig)
                elif subindex2 + 1 < tig_index_list_len:

                    subindex3 = subindex2 + 1

                    row3 = df.loc[tig_index_list[subindex3]]

                    # Determine if the three alignment records have an inversion signature
                    if (
                        row3['REV'] == row1['REV']
                    ) and ( # Middle alignment record must be between the flanking alignment records on tnhe contig
                        (
                            not row1['REV'] and (row1['QUERY_TIG_END'] < (row2['QUERY_TIG_POS'] + row2['QUERY_TIG_END']) // 2 < row3['QUERY_TIG_POS'])
                        ) or (
                            row1['REV'] and (row3['QUERY_TIG_POS'] < (row2['QUERY_TIG_POS'] + row2['QUERY_TIG_END']) // 2 < row1['QUERY_TIG_END'])
                        )
                    ):

                        # Classic INV signature: +,-,+ or -,+,-
                        # subindex1 and subindex2 are in opposite directions (already tested)

                        # Find inversion in 3-part alignment with middle in opposite orientation as flanks
                        region_flag = pavlib.seq.Region(chrom, row1['END'], row3['POS'], is_rev=row1['REV'])

                        # Try INV
                        inv_call = pavlib.inv.scan_for_inv(
                            region_flag,
                            ref_fa_name, tig_fa_name,
                            align_lift, k_util,
                            max_region_size=max_region_size,
                            threads=threads,
                            n_tree=n_tree,
                            srs_tree=srs_tree,
                            log=log,
                            min_exp_count=1  # If alignment truncation does not contain inverted k-mers at sufficient density, stop searching
                        )

                        # Recover inversion if alignment supports and density fails.
                        if inv_call is None and subindex2 == subindex1 + 1 and subindex3 == subindex1 + 2:
                            region_ref = pavlib.seq.Region(chrom, row2['POS'], row2['END'])
                            region_tig = pavlib.seq.Region(row2['QUERY_ID'], row2['QUERY_TIG_POS'], row2['QUERY_TIG_END'])

                            inv_call = pavlib.inv.InvCall(
                                region_ref, region_ref,
                                region_tig, region_tig,
                                region_ref, region_tig,
                                region_ref,
                                None
                            )

                            call_source = CALL_SOURCE_INV_NO_DENSITY

                        else:
                            call_source = CALL_SOURCE_INV_DENSITY

                        if inv_call is not None:
                            log.write('INV (3-tig): {}\n'.format(inv_call))
                            log.flush()

                            seq = pavlib.seq.region_seq_fasta(
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

                                    hap,

                                    inv_call.region_tig_outer.to_base1_string(),
                                    '-' if is_rev else '+',

                                    0,

                                    inv_call.region_ref_inner.to_base1_string(),
                                    inv_call.region_tig_inner.to_base1_string(),

                                    inv_call.region_ref_discovery.to_base1_string(),
                                    inv_call.region_tig_discovery.to_base1_string(),

                                    inv_call.region_flag.region_id(),
                                    'ALNTRUNC',

                                    '{},{},{}'.format(row1['INDEX'], row2['INDEX'], row3['INDEX']),
                                    row1['CLUSTER_MATCH'],

                                    call_source,

                                    seq
                                ],
                                index=[
                                    '#CHROM', 'POS', 'END',
                                    'ID', 'SVTYPE', 'SVLEN',
                                    'HAP',
                                    'TIG_REGION', 'QUERY_STRAND',
                                    'CI',
                                    'RGN_REF_INNER', 'RGN_TIG_INNER',
                                    'RGN_REF_DISC', 'RGN_TIG_DISC',
                                    'FLAG_ID', 'FLAG_TYPE',
                                    'ALIGN_INDEX', 'CLUSTER_MATCH',
                                    'CALL_SOURCE',
                                    'SEQ'
                                ]
                            ))

                            # Save density table
                            if density_out_dir is not None and inv_call.df is not None:
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
        df_ins = pd.concat(ins_list, axis=1).T
        df_ins['ID'] = svpoplib.variant.version_id(df_ins['ID'])
        df_ins.sort_values(['#CHROM', 'POS', 'END', 'ID'], inplace=True)

    else:
        df_ins = pd.DataFrame([], columns=[
            '#CHROM', 'POS', 'END',
            'ID', 'SVTYPE', 'SVLEN',
            'HAP',
            'TIG_REGION', 'QUERY_STRAND',
            'CI',
            'ALIGN_INDEX', 'CLUSTER_MATCH',
            'LEFT_SHIFT', 'HOM_REF', 'HOM_TIG',
            'CALL_SOURCE',
            'SEQ'
        ])

    if len(del_list) > 0:
        df_del = pd.concat(del_list, axis=1).T
        df_del['ID'] = svpoplib.variant.version_id(df_del['ID'])
        df_del.sort_values(['#CHROM', 'POS', 'END', 'ID'], inplace=True)

    else:
        df_del = pd.DataFrame([], columns=[
            '#CHROM', 'POS', 'END',
            'ID', 'SVTYPE', 'SVLEN',
            'HAP',
            'TIG_REGION', 'QUERY_STRAND',
            'CI',
            'ALIGN_INDEX', 'CLUSTER_MATCH',
            'LEFT_SHIFT', 'HOM_REF', 'HOM_TIG',
            'CALL_SOURCE',
            'SEQ'
        ])

    if len(inv_list) > 0:
        df_inv = pd.concat(inv_list, axis=1).T
        df_inv['ID'] = svpoplib.variant.version_id(df_inv['ID'])
        df_inv.sort_values(['#CHROM', 'POS', 'END', 'ID'], inplace=True)

    else:
        df_inv = pd.DataFrame([], columns=[
            '#CHROM', 'POS', 'END',
            'ID', 'SVTYPE', 'SVLEN',
            'HAP',
            'TIG_REGION', 'QUERY_STRAND',
            'CI',
            'RGN_REF_INNER', 'RGN_TIG_INNER',
            'RGN_REF_DISC', 'RGN_TIG_DISC',
            'FLAG_ID', 'FLAG_TYPE',
            'ALIGN_INDEX', 'CLUSTER_MATCH',
            'CALL_SOURCE',
            'SEQ'
        ])

    # Return records
    return df_ins, df_del, df_inv


class SeqCache:
    """
    Keep a cache of a sequence string in upper-case. Stores the last instance of the sequence and the ID. When a
    new ID is requested, the old sequnece is discarded and the new one is loaded.
    """

    def __init__(self, fa_filename, uppercase=True):
        """
        Create a cache object to read from indexed FASTA file `fa_filename`.

        :param fa_filename: Indexed FASTA file name.
        :param uppercase: `True` if sequences should be made upper-case, otherwise, preserve case.
        """

        self.fa_filename = fa_filename
        self.uppercase = uppercase

        self.id = None
        self.seq = None
        self.is_rev = None

    def get(self, sequence_id, is_rev):
        """
        Get a sequence. Returns the cached version if ID matches the current ID, otherwise, the correct sequence is
        retrieved, cached, and returned.

        :param sequence_id: Sequence ID string or Region.
        :param is_rev: `True` if the sequence is reverse-complemented. Retrieving the same sequence ID as the cached
            sequence with `is_rev` mismatch will reload the sequence in the requested orientation.

        :return: Sequence.
        """

        if self.id != sequence_id or self.is_rev != is_rev:
            new_seq = pavlib.seq.region_seq_fasta(
                sequence_id, self.fa_filename
            )

            if self.uppercase:
                new_seq = new_seq.upper()

            if is_rev:
                new_seq = str(Bio.Seq.Seq(new_seq).reverse_complement())

            self.seq = new_seq

            self.id = sequence_id
            self.is_rev = is_rev

        return self.seq
