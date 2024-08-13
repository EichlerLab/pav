"""
Variant anchor.
"""

import numpy as np
import pandas as pd
import pavlib

SEGMENT_TABLE_FIELDS = [  # Fields for creating an empty table
    '#CHROM', 'POS', 'END',
    'IS_ANCHOR', 'IS_ALIGNED', 'IS_REV',
    'QRY_ID', 'QRY_POS', 'QRY_END',
    'LEN_REF', 'LEN_QRY',
    'GAP_REF',
    'INDEX', 'SEG_ORDER',
    'SCORE',
    'TRIM_REF_L', 'TRIM_REF_R', 'TRIM_QRY_L', 'TRIM_QRY_R'
]

class AnchoredInterval:
    """
    Represents one interval of alignments spanning from one anchor to another. This object is fed into call methods.

    :param chain_node: A ChainNode object describing the start and end index positions of the chain anchors in query-
        sorted order.
    :param df_align: DataFrame of alignment records with the 'SCORE' column added.
    :param score_model: Model for scoring alignments and template switches.

    Attributes:
        * chain_node: AnchorChainNode object used to create this interval.
        * is_rev: True If both anchors are reverse-complemented and False if both are not (error if they don't match).
        * qry_id: Query ID.
        * chrom: Reference chromosome name.
        * df_segment: Dataframe of segments with the anchors in the first and last positions.
        * seg_n: Number of segments between the anchors (length of df_segment without anchors).
        * len_ref: Distance between anchors on the reference. If negative, then the switch from anchor-to-anchor goes
            backwards in reference space. "ref_region" has the absolute value of the distance.
        * len_qry: Distance between anchors on the query.
        * ref_region: Reference region (pavlib.seq.Region).
        * qry_region: Query region (pavlib.seq.Region).
    """

    def __init__(self, chain_node, df_align, score_model):
        self.chain_node = chain_node

        segment_list = list()  # List of query segments broken by template switches and untemplated insertions

        if chain_node.end_index <= chain_node.start_index:
            raise RuntimeError(f'Chain node end index {chain_node.end_index} must be greater than the start index {chain_node.start_index}')

        head_row = df_align.loc[chain_node.start_index]
        tail_row = df_align.loc[chain_node.end_index]

        self.is_rev = head_row['REV']
        self.qry_id = head_row['QRY_ID']
        self.chrom = head_row['#CHROM']

        if self.is_rev != tail_row['REV'] or self.qry_id != tail_row['QRY_ID'] or self.chrom != tail_row['#CHROM']:
            raise RuntimeError(f'Anchors {chain_node.start_index} (index={head_row["INDEX"]}) and {chain_node.end_index} (index={tail_row["INDEX"]}): Non-matching REV, QRY_ID, or CHROM')

        # Create templated and untemplated insertions for each aligned segment and transition.
        # Save row_l outside the loop, row_r may be transformed inside the loop and should be preserved for the next round
        # as row_r.
        index_l = chain_node.start_index
        row_l = df_align.loc[index_l]

        index_r = index_l + 1

        seg_order = 0

        # Add head anchor
        segment_list.append(
            pd.Series(
                [
                    head_row['#CHROM'], head_row['POS'], head_row['END'],
                    True, True, head_row['REV'],
                    self.qry_id, head_row['QRY_POS'], head_row['QRY_END'],
                    head_row['END'] - head_row['POS'], head_row['QRY_END'] - head_row['QRY_POS'],
                    0,
                    head_row['INDEX'], seg_order,
                    head_row['SCORE'],
                    head_row['TRIM_REF_L'], head_row['TRIM_REF_R'], head_row['TRIM_QRY_L'], head_row['TRIM_QRY_R']
                ],
                index=[
                    '#CHROM', 'POS', 'END',
                    'IS_ANCHOR', 'IS_ALIGNED', 'IS_REV',
                    'QRY_ID', 'QRY_POS', 'QRY_END',
                    'LEN_REF', 'LEN_QRY',
                    'GAP_REF',
                    'INDEX', 'SEG_ORDER',
                    'SCORE',
                    'TRIM_REF_L', 'TRIM_REF_R', 'TRIM_QRY_L', 'TRIM_QRY_R'
                ]
            )
        )

        seg_order += 1

        # Traverse alignment records setting segments
        while index_r <= chain_node.end_index:

            # Get right row
            row_r = df_align.loc[index_r]

            if row_r['QRY_ID'] != self.qry_id:
                raise RuntimeError(f'Query ID in alignment record {row_r["INDEX"]} ("{row_r["QRY_ID"]}") does not match the query ID for this chain "{self.qry_id}"')

            # Get query and reference positions
            gap_qry_pos = row_l['QRY_END']
            gap_qry_end = row_r['QRY_POS']

            if self.is_rev:  # Get rows in reference order
                row_l_ref = row_r
                row_r_ref = row_l
            else:
                row_l_ref = row_l
                row_r_ref = row_r

            gap_ref_pos = row_l_ref['END']
            gap_ref_end = row_r_ref['POS']

            gap_len_qry = gap_qry_end - gap_qry_pos
            gap_len_ref = gap_ref_end - gap_ref_pos

            if gap_len_qry < 0:
                raise RuntimeError(f'Got negative gap between alignment records {index_l} and {index_r}: {gap_len_qry}')

            if gap_len_qry > 0:
                # Add gap between template switches as an untemplated insertion
                segment_list.append(
                    pd.Series(
                        [
                            np.nan, -1, -1,
                            False, False, False,
                            self.qry_id, gap_qry_pos, gap_qry_end,
                            0, gap_len_qry,
                            0,
                            -1, seg_order,
                            score_model.gap(gap_len_qry),
                            0, 0, 0, 0
                        ],
                        index=[
                            '#CHROM', 'POS', 'END',
                            'IS_ANCHOR', 'IS_ALIGNED', 'IS_REV',
                            'QRY_ID', 'QRY_POS', 'QRY_END',
                            'LEN_REF', 'LEN_QRY',
                            'GAP_REF',
                            'INDEX', 'SEG_ORDER',
                            'SCORE',
                            'TRIM_REF_L', 'TRIM_REF_R', 'TRIM_QRY_L', 'TRIM_QRY_R'
                        ]
                    )
                )

                seg_order += 1

            # Add templated segment
            if index_r != chain_node.end_index:
                segment_list.append(
                    pd.Series(
                        [
                            row_r['#CHROM'], row_r['POS'], row_r['END'],
                            False, True, row_r['REV'],
                            self.qry_id, row_r['QRY_POS'], row_r['QRY_END'],
                            row_r['END'] - row_r['POS'], row_r['QRY_END'] - row_r['QRY_POS'],
                            gap_len_ref,
                            row_r['INDEX'], seg_order,
                            row_r['SCORE'],
                            row_r['TRIM_REF_L'], row_r['TRIM_REF_R'], row_r['TRIM_QRY_L'], row_r['TRIM_QRY_R']
                        ],
                        index=[
                            '#CHROM', 'POS', 'END',
                            'IS_ANCHOR', 'IS_ALIGNED', 'IS_REV',
                            'QRY_ID', 'QRY_POS', 'QRY_END',
                            'LEN_REF', 'LEN_QRY',
                            'GAP_REF',
                            'INDEX', 'SEG_ORDER',
                            'SCORE',
                            'TRIM_REF_L', 'TRIM_REF_R', 'TRIM_QRY_L', 'TRIM_QRY_R'
                        ]
                    )
                )

            seg_order += 1

            # Next aligned segment
            index_l = index_r
            index_r = index_l + 1

            row_l = row_r

        # Add tail anchor
        segment_list.append(
            pd.Series(
                [
                    tail_row['#CHROM'], tail_row['POS'], tail_row['END'],
                    True, True, tail_row['REV'],
                    self.qry_id, tail_row['QRY_POS'], tail_row['QRY_END'],
                    tail_row['END'] - tail_row['POS'], tail_row['QRY_END'] - tail_row['QRY_POS'],
                    0,
                    tail_row['INDEX'], seg_order,
                    tail_row['SCORE'],
                    tail_row['TRIM_REF_L'], tail_row['TRIM_REF_R'], tail_row['TRIM_QRY_L'], tail_row['TRIM_QRY_R']
                ],
                index=[
                    '#CHROM', 'POS', 'END',
                    'IS_ANCHOR', 'IS_ALIGNED', 'IS_REV',
                    'QRY_ID', 'QRY_POS', 'QRY_END',
                    'LEN_REF', 'LEN_QRY',
                    'GAP_REF',
                    'INDEX', 'SEG_ORDER',
                    'SCORE',
                    'TRIM_REF_L', 'TRIM_REF_R', 'TRIM_QRY_L', 'TRIM_QRY_R'
                ]
            )
        )

        seg_order += 1

        # Concat segments, order by query.
        # Set complex reference position start/end by the anchoring alignment records.
        self.df_segment = pd.concat(
            segment_list, axis=1
        ).T.astype({
            '#CHROM': str, 'POS': int, 'END': int,
            'IS_ANCHOR': bool, 'IS_ALIGNED': bool, 'IS_REV': bool,
            'QRY_ID': str, 'QRY_POS': int, 'QRY_END': int,
            'LEN_REF': int, 'LEN_QRY': int,
            'GAP_REF': int,
            'INDEX': int, 'SEG_ORDER': int,
            'SCORE': float
        })

        self.df_segment['STRAND'] = self.df_segment.apply(lambda row:
            '*' if not row['IS_ALIGNED'] else ({True: '-', False: '+'}[row['IS_REV']]),
            axis=1
        )

        # Get query and reference regions
        qry_pos = self.df_segment.iloc[0]['QRY_END']
        qry_end = self.df_segment.iloc[-1]['QRY_POS']

        if self.is_rev:
            self.df_segment = self.df_segment.iloc[::-1]

        ref_pos = self.df_segment.iloc[0]['END']
        ref_end = self.df_segment.iloc[-1]['POS']

        self.seg_n = self.df_segment.iloc[1:-1].shape[0]

        self.len_ref = ref_end - ref_pos

        if ref_end - ref_pos < 0:
            self.ref_region = pavlib.seq.Region(self.chrom, ref_end, ref_pos)
            # ref_pos, ref_end = ref_end, ref_pos
        else:
            self.ref_region = pavlib.seq.Region(self.chrom, ref_pos, ref_end)

        # Get query region
        self.qry_region = pavlib.seq.Region(self.qry_id, qry_pos, qry_end, self.is_rev)
        self.ref_region = pavlib.seq.Region(self.chrom, ref_pos, ref_end)

        self.len_qry = len(self.qry_region)

        # Test invariants
        if self.len_qry < 0:
            raise ValueError(f'Qry region length is negative: {self.qry_region}')

        if self.seg_n != self.df_segment.shape[0] - 2:
            raise ValueError(f'Segment count mismatch: seg_n={self.seg_n} != segments - 2 ({self.df_segment.shape[0] - 2})')

    def __repr__(self):
        return f'AnchoredInterval(start={self.chain_node.start_index}, end={self.chain_node.end_index}, chr={self.ref_region.chrom}, query={self.qry_region.chrom}, seg_n={self.seg_n}, is_rev={self.is_rev})'
