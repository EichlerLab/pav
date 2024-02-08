# Alignment lift

import collections
import intervaltree
import numpy as np

import pavlib

from .align import *


class AlignLift:
    """
    Create an alignment liftover object for translating between reference and contig alignments (both directions). Build
    liftover from alignment data in a DataFrame (requires #CHROM, POS, END, QRY_ID, QUERY_POS, and QUERY_END). The
    DataFrame index must not have repeated values. The data-frame is also expected not to change, so feed this object
    a copy if it might be altered while this object is in use.
    """

    def __init__(self, df, df_fai, cache_align=10):
        """
        Create AlignLift instance.

        :param df: Alignment BED file (post-cut) as DataFrame.
        :param df_fai: Contig FAI file as Series.
        :param cache_align: Number of alignment records to cache.
        """

        self.df = df
        self.df_fai = df_fai
        self.cache_align = cache_align

        # Check df
        if len(set(df.index)) != df.shape[0]:
            raise RuntimeError('Cannot create AlignLift object with duplicate index values')

        # Build a reference-coordinate and a tig-coordinate tree
        self.ref_tree = collections.defaultdict(intervaltree.IntervalTree)
        self.tig_tree = collections.defaultdict(intervaltree.IntervalTree)

        for index, row in df.iterrows():
            self.ref_tree[row['#CHROM']][row['POS']:row['END']] = index
            self.tig_tree[row['QRY_ID']][row['QRY_POS']:row['QRY_END']] = index

        # Build alignment caching structures
        self.cache_queue = collections.deque()

        self.ref_cache = dict()
        self.tig_cache = dict()

    def lift_to_sub(self, query_id, coord, gap=False):
        """
        Lift coordinates from query (tig) to subject (reference).

        Returns tuple(s) of:
            0: Query ID.
            1: Position.
            2: Is reverse complemented. `True` if the alignment was lifted through a reverse-complemented
                alignment record.
            3: Minimum position.
            4: Maximum position.
            5: Align record index.

        If the lift hits an alignment gap, then the function can try to resolve the lift to the best subject position
        based on the alignment using a strategy outlined by the `gap` parameter.

        :param query_id: Query record ID.
        :param coord: Query coordinates. May be a single value, list, or tuple.
        :param gap: Interpolate into an alignment gap between two contigs if `True`.

        :return: Single coordinate tuple or list of coordinate tuples if `coord` is a list or tuple. Returns `None`
            for coordinates that cannot be lifted.
        """

        # Determine type
        if issubclass(coord.__class__, list) or issubclass(coord.__class__, tuple):
            ret_list = True
        else:
            ret_list = False
            coord = (coord,)

        # Do lift
        lift_coord_list = list()

        for pos in coord:
            pos_org = pos  # Pre-reverse position

            # Find matching records
            match_set = self.tig_tree[query_id][pos:(pos + 1)]

            if len(match_set) == 1:
                index = list(match_set)[0].data

            elif len(match_set) == 0 and gap:
                lift_coord_list.append(self._get_subject_gap(query_id, pos))
                continue

            else:
                lift_coord_list.append(None)
                continue

            # Get lift tree
            if index not in self.tig_cache.keys():
                self._add_align(index)

            lift_tree = self.tig_cache[index]

            # Get row
            row = self.df.loc[index]

            # Reverse coordinates of pos if the alignment is reverse-complemented. Translates from QRY_POS-space
            # to QUERY_POS-space.
            if row['REV']:
                pos = self.df_fai[query_id] - pos

            # Get match record
            match_set = lift_tree[pos]

            if len(match_set) == 1:
                match_interval = list(match_set)[0]

            if len(match_set) == 0:
                # Allow queries to match if they end exactly at the alignment end
                match_set = lift_tree[pos - 1]

                match_interval = list(match_set)[0] if len(match_set) == 1 else None

                if not match_interval or match_interval.end != pos:
                    lift_coord_list.append(None)

                    raise RuntimeError(
                        (
                            'Found no matches in a lift-tree for a record within a '
                            'global to-subject tree: {}:{} (index={}, gap={})'
                        ).format(query_id, pos_org, index, gap)
                    )

            elif len(match_set) > 1:
                lift_coord_list.append(None)

                raise RuntimeError(
                    (
                        'Found multiple matches in a lift-tree for a record within a '
                        'global to-subject tree: {}:{} (index={}, gap={})'
                    ).format(query_id, pos_org, index,  gap)
                )

            # Interpolate coordinates
            if match_interval.data[1] - match_interval.data[0] > 1:
                lift_pos = match_interval.data[0] + (pos - match_interval.begin)

                lift_coord_list.append((
                    row['#CHROM'],
                    lift_pos,
                    row['REV'],
                    lift_pos,
                    lift_pos,
                    (row['INDEX'],)
                ))

            else:  # Lift from missing bases on the target (insertion or deletion)
                lift_coord_list.append((
                    row['#CHROM'],
                    match_interval.data[1],
                    row['REV'],
                    match_interval.data[1],
                    match_interval.data[1],
                    (row['INDEX'],)
                ))

        # Return coordinates
        if ret_list:
            return lift_coord_list
        else:
            return lift_coord_list[0]

    def lift_to_qry(self, subject_id, coord):
        """
        Lift coordinates from subject (reference) to query (tig).

        Returns tuple(s) of:
            0: Query ID.
            1: Position.
            2: Is reverse complemented. `True` if the alignment was lifted through a reverse-complemented
                alignment record.
            3: Minimum position.
            4: Maximum position.
            5: Align record index.

        :param subject_id: Subject ID.
        :param coord: Subject coordinates. May be a single value, list, or tuple.

        :return: Single coordinate tuple or list of coordinate tuples if `coord` is a list or tuple. Returns `None`
            for coordinates that cannot be lifted.
        """

        # Determine type
        if issubclass(coord.__class__, list) or issubclass(coord.__class__, tuple):
            ret_list = True
        else:
            ret_list = False
            coord = (coord,)

        # Do lift
        lift_coord_list = list()

        for pos in coord:
            match_set = self.ref_tree[subject_id][pos:(pos + 1)]

            # Check coordinates
            if len(match_set) == 0:
                lift_coord_list.append(None)
                continue

                # raise ValueError('Subject region {}:{} has no to-query lift records'.format(subject_id, pos))

            if len(match_set) > 1:
                lift_coord_list.append(None)
                continue

                # raise ValueError(
                #     'Subject region {}:{} has {} to-query lift records'.format(subject_id, pos, len(match_set))
                # )

            # Get lift tree
            index = list(match_set)[0].data

            if index not in self.ref_cache.keys():
                self._add_align(index)

            lift_tree = self.ref_cache[index]

            # Save row
            row = self.df.loc[index]

            # Get match record
            match_set = lift_tree[pos:(pos + 1)]

            if len(match_set) != 1:
                raise RuntimeError(
                    (
                        'Program bug: Found no matches in a lift-tree for a record withing a '
                        'global to-query tree: {}:{} (index={})'
                    ).format(subject_id, pos, index)
                )

            match_interval = list(match_set)[0]

            # Interpolate coordinates
            if match_interval.data[1] - match_interval.data[0] > 1:
                qry_pos = match_interval.data[0] + (pos - match_interval.begin)

            else:  # Lift from missing bases on the target (insertion or deletion)
                qry_pos = match_interval.data[1]

            if row['REV']:
                qry_pos = self.df_fai[row['QRY_ID']] - qry_pos

            lift_coord_list.append((
                row['QRY_ID'],
                qry_pos,
                row['REV'],
                qry_pos,
                qry_pos,
                (row['INDEX'],)
            ))

        # Return coordinates
        if ret_list:
            return lift_coord_list
        else:
            return lift_coord_list[0]

    def lift_region_to_sub(self, region, gap=False):
        """
        Lift region to subject.

        :param region: Query region.
        :param gap: Interpolate within gap if `True`.

        :return: Subject region or `None` if it could not be lifted.
        """

        # Lift
        sub_pos, sub_end = self.lift_to_sub(region.chrom, (region.pos, region.end), gap)

        # Check lift: Must lift both ends to the same subject ID
        if sub_pos is None or sub_end is None:
            return None

        if sub_pos[0] != sub_end[0] or (sub_pos[2] is not None and sub_end[2] is not None and sub_pos[2] != sub_end[2]):
            return None

        # Return
        return pavlib.seq.Region(
            sub_pos[0], sub_pos[1], sub_end[1],
            is_rev=False,
            pos_min=sub_pos[3], pos_max=sub_pos[4],
            end_min=sub_end[3], end_max=sub_end[4],
            pos_aln_index=(sub_pos[5],),
            end_aln_index=(sub_end[5],)
        )

    def lift_region_to_qry(self, region):
        """
        Lift region to query.

        :param region: Subject region.

        :return: Query region or `None` if it could not be lifted.
        """

        # Lift
        query_pos, query_end = self.lift_to_qry(region.chrom, (region.pos, region.end))

        # Check lift: Must lift both ends to the same query ID
        if query_pos is None or query_end is None:
            return None

        if query_pos[0] != query_end[0] or query_pos[2] != query_end[2]:
            return None

        # Return
        return pavlib.seq.Region(
            query_pos[0], query_pos[1], query_end[1],
            is_rev=query_pos[2],
            pos_min=query_pos[3], pos_max=query_pos[4],
            end_min=query_end[3], end_max=query_end[4],
            pos_aln_index=(query_pos[5],),
            end_aln_index=(query_end[5],)
        )

    def _get_subject_gap(self, query_id, pos):
        """
        Interpolate lift coordinates to an alignment gap.

        :param pos: Position on the contig.

        :return: A tuple of (subject, pos, rev, min, max).
        """

        # Check arguments
        if pos is None:
            return None

        # Get alignment records for this contig
        subdf = self.df.loc[self.df['QRY_ID'] == query_id]

        # Must be flanked by two contigs on either side
        if not np.any(subdf['QRY_END'] < pos):
            return None

        if not np.any(subdf['QRY_POS'] > pos):
            return None

        # Get left and right rows for the alignment record flanking this position
        row_l = subdf.loc[
            subdf.loc[subdf['QRY_END'] < pos, 'QRY_END'].sort_values().index[-1]
        ]

        row_r = subdf.loc[
            subdf.loc[subdf['QRY_POS'] > pos, 'QRY_POS'].sort_values().index[0]
        ]

        # Rows must be mapped to the same subject
        if row_l['#CHROM'] != row_r['#CHROM']:
            return None

        return (
            (
                row_l['#CHROM'],
                int((row_l['QRY_END'] + row_r['QRY_POS']) / 2),
                row_l['REV'] if row_l['REV'] == row_r['REV'] else None,
                row_l['QRY_END'],
                row_r['QRY_POS'],
                (row_l['INDEX'], row_r['INDEX'])
            )
        )

    def _add_align(self, index):
        """
        Add an alignment from DataFrame index `index`.

        :param index: DataFrame index.
        """

        # No alignment to add if it's already cached.
        if index in self.ref_cache.keys():

            # Append to end (last used) and skip adding alignment
            while index in self.cache_queue:
                self.cache_queue.remove(index)

            self.cache_queue.appendleft(index)

            return

        # Make space for this alignment
        self._check_and_clear()

        # Get row
        row = self.df.loc[index]

        # Build lift trees
        sub_bp = row['POS']
        qry_bp = 0

        itree_ref = intervaltree.IntervalTree()
        itree_qry = intervaltree.IntervalTree()

        # Get CIGAR and check query start position
        cigar_op_list = list(cigar_str_to_tuples(self.df.loc[index]))

        clipped_bp = 0

        cigar_index = 0

        while cigar_index < len(cigar_op_list) and cigar_op_list[cigar_index][1] in {'S', 'H'}:
            clipped_bp += cigar_op_list[cigar_index][0]
            cigar_index += 1

        # qry_map_rgn = pavlib.seq.region_from_string(row['QRY_MAP_RGN'])
        #
        # if qry_map_rgn.pos != clipped_bp:
        #     raise RuntimeError(
        #         'Number of clipped bases ({}) does not match query start position {}: Alignment {}:{} ({}:{})'.format(
        #             clipped_bp, qry_map_rgn.pos, row['#CHROM'], row['POS'], row['QRY_ID'], qry_map_rgn.pos
        #         )
        #     )

        # Build trees
        for cigar_len, cigar_op in cigar_op_list:

            if cigar_op in {'=', 'X', 'M'}:

                itree_ref[sub_bp:(sub_bp + cigar_len)] = (qry_bp, qry_bp + cigar_len)
                itree_qry[qry_bp:(qry_bp + cigar_len)] = (sub_bp, sub_bp + cigar_len)

                sub_bp += cigar_len
                qry_bp += cigar_len

            elif cigar_op == 'I':

                itree_qry[qry_bp:(qry_bp + cigar_len)] = (sub_bp, sub_bp + 1)

                qry_bp += cigar_len

            elif cigar_op == 'D':

                itree_ref[sub_bp:(sub_bp + cigar_len)] = (qry_bp, qry_bp + 1)

                sub_bp += cigar_len

            elif cigar_op in {'S', 'H'}:

                qry_bp += cigar_len

                clipped_bp += cigar_len

            else:
                row = self.df.loc[index]

                raise RuntimeError('Unhandled CIGAR operation: {}: Alignment {}:{} ({})'.format(
                    cigar_op, row['#CHROM'], row['POS'], row['QRY_ID']
                ))

                # raise RuntimeError('Unhandled CIGAR operation: {}: Alignment {}:{} ({}:{})'.format(
                #     cigar_op, row['#CHROM'], row['POS'], row['QRY_ID'], qry_map_rgn.pos
                # ))

        # Cache trees
        self.ref_cache[index] = itree_ref
        self.tig_cache[index] = itree_qry

        # Add index to end of queue
        self.cache_queue.appendleft(index)

    def _check_and_clear(self):
        """
        Check alignment cache and clear if necessary to make space.
        """

        while len(self.cache_queue) >= self.cache_align:
            index = self.cache_queue.pop()

            del(self.ref_cache[index])
            del(self.tig_cache[index])
