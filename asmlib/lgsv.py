"""
Call large alignment-truncating variants.
"""

# Default parameters
MAX_TIG_DIST_PROP = 1  # Max allowed tig gap as a factor of the minimum alignment length of two records
MAX_REF_DIST_PROP = 3  # Max allowed ref gap as a factor of the minimum alignment length of two records

CALL_SOURCE = 'ALNTRUNC'


def scan_for_events(df, hap, max_tig_dist_prop=None, max_ref_dist_prop=None):
    """
    Scan trimmed alignments for alignment-truncating SV events.

    :param df: Dataframe of alignments (trimmed). Output from "align_cut_tig_overlap".
    :param hap: Haplotype.
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

            while subindex2 < tig_index_list_len:
                row1 = df.loc[tig_index_list[subindex1]]
                row2 = df.loc[tig_index_list[subindex2]]

                dist_tig = row2['QUERY_POS'] - row1['QUERY_END']
                dist_ref = row2['POS'] - row1['END']
                min_aln_len = np.min([row1['QRY_LEN'], row2['QRY_LEN']])

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

                    # Stop if there is large gap between the aligned bases (tig or ref)
                    if (
                            np.abs(dist_tig) / min_aln_len > max_tig_dist_prop
                    ) or (
                            np.abs(dist_ref) / min_aln_len > max_ref_dist_prop
                    ):
                        subindex2 += 1
                        continue

                    if dist_tig >= 50 or dist_ref >= 50:

                        if dist_tig == 0:
                            # DEL

                            # Append
                            del_list.append(pd.Series(
                                [
                                    chrom, pos, end, id, 'SV', 'DEL', dist_ref,
                                    hap,
                                    tig_id, tig_pos, tig_end, tig_strand,
                                    CALLSOURCE, row1['CLUSTER_MATCH'],
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

                        elif dist_ref == 0:

                            #tig_pos = row1['QUERY_END']  # Correct for rev alignment
                            #tig_end = row2['QUERY_POS']
                            #seq =

                            ins_list.append(pd.Series(
                                [
                                    chrom, row1['END'], row1['END'] + 1, id, 'SV', 'INS', dist_tig,
                                    hap,
                                    tig_id, tig_pos, tig_end, '-' if row1['REV'] else '+',
                                    CALL_SOURCE,
                                    seq
                                ],
                                index=[
                                    '#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN',
                                    'HAP',
                                    'QUERY_ID', 'QUERY_POS', 'QUERY_END', 'QUERY_STRAND',
                                    'CALL_SOURCE',
                                    'SEQ'
                                ]
                            ))

                    # Next record
                    subindex2 += 1

                else:

                    subindex3 = subindex2 + 1

                    while subindex3 < tig_index_list_len:
                        row3 = df.loc[tig_index_list[subindex3]]

                        if row3['REV'] == row1['REV']:
                            #                        print('\t* Try INV: {}:{}: {}, {} ({}, index={},{},{})'.format(
                            #                            row1['#CHROM'], row1['END'], dist_tig, dist_ref, row1['QUERY_ID'],
                            #                            tig_index_list[subindex1],
                            #                            tig_index_list[subindex2],
                            #                            tig_index_list[subindex3]
                            #                        ))

                            # Save record
                            flag_list.append(pd.Series(  # INV flag
                                [
                                    'INV',
                                    row1['#CHROM'], row1['QUERY_ID'],
                                    row1['END'], row2['POS'],
                                    dist_tig, dist_ref,
                                    tig_index_list[subindex1], tig_index_list[subindex2], tig_index_list[subindex3],
                                    subindex1, subindex2, subindex3,
                                    '{}:{:d}-{:d}'.format(row1['#CHROM'], row1['END'] - 10000, row3['POS'] + 10000)
                                ],
                                index=[
                                    'TYPE',
                                    '#CHROM', 'QUERY_ID',
                                    'REF_ST_UP', 'REF_ST_DN',
                                    'DIST_TIG', 'DIST_REF',
                                    'TIG_IDX_1', 'TIG_IDX_2', 'TIG_IDX_3',
                                    'SUB_IDX_1', 'SUB_IDX_2', 'SUB_IDX_3',
                                    'WIN_10K'
                                ]
                            ))

                            # TODO: Call INV event

                        # Next subindex3
                        subindex3 += 1

                    # Next subindex2
                    subindex2 += 1
