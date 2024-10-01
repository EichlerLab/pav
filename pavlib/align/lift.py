# Alignment lift

import collections
import intervaltree
import numpy as np

import pavlib


class AlignLift:
    """
    Create an alignment liftover object for translating between reference and query alignments (both directions). Build
    liftover from alignment data in a DataFrame (requires #CHROM, POS, END, QRY_ID, QUERY_POS, and QUERY_END). The
    DataFrame index must not have repeated values. The data-frame is also expected not to change, so feed this object
    a copy if it might be altered while this object is in use.
    """

    def __init__(self, df, df_fai, cache_align=10):
        """
        Create AlignLift instance.

        :param df: Alignment BED file (post-cut) as DataFrame.
        :param df_fai: Query FAI file as Series.
        :param cache_align: Number of alignment records to cache.
        """

        self.df = df
        self.df_fai = df_fai
        self.cache_align = cache_align

        # Check df
        if len(set(df.index)) != df.shape[0]:
            raise RuntimeError('Cannot create AlignLift object with duplicate index values')

        # Build a reference-coordinate and a query-coordinate tree
        self.ref_tree = collections.defaultdict(intervaltree.IntervalTree)
        self.qry_tree = collections.defaultdict(intervaltree.IntervalTree)

        for index, row in df.iterrows():
            self.ref_tree[row['#CHROM']][row['POS']:row['END']] = index
            self.qry_tree[row['QRY_ID']][row['QRY_POS']:row['QRY_END']] = index

        # Build alignment caching structures
        self.cache_queue = collections.deque()

        self.ref_cache = dict()
        self.qry_cache = dict()

    def lift_to_ref(self, query_id, coord, gap=False):
        """
        Lift coordinates from query to reference.

        Returns tuple(s) of:
            0: Query ID.
            1: Position.
            2: Is reverse complemented. `True` if the alignment was lifted through a reverse-complemented
                alignment record.
            3: Minimum position.
            4: Maximum position.
            5: Align record index.
            6: Alignment record score.

        If the lift hits an alignment gap, then the function can try to resolve the lift to the best reference position
        based on the alignment using a strategy outlined by the `gap` parameter.

        :param query_id: Query record ID.
        :param coord: Query coordinates. May be a single value, list, or tuple.
        :param gap: Interpolate into an alignment gap between two queries if `True`.

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
            pos_org = pos  # Pre-reverse position for error reporting

            # Find matching records
            match_set = self.qry_tree[query_id][pos:(pos + 1)]

            if len(match_set) == 0:
                if gap:
                    lift_coord_list.append(self._get_ref_gap(query_id, pos))
                else:
                    lift_coord_list.append(None)
                continue

            # Get list of lift matches
            match_list = list()

            for match_element in match_set:
                index = match_element.data

                # Get lift tree
                if index not in self.qry_cache.keys():
                    self._add_align(index)

                lift_tree = self.qry_cache[index]

                # Get row
                row = self.df.loc[index]

                # Reverse coordinates of pos if the alignment is reverse-complemented. Translates from QRY_POS-space
                # to QUERY_POS-space.
                if row['REV']:
                    pos = self.df_fai[query_id] - pos

                # Get match record
                lift_set = lift_tree[pos]

                if len(lift_set) == 1:
                    match_interval = list(lift_set)[0]

                elif len(lift_set) == 0:
                    # Allow queries to match if they end exactly at the alignment end
                    lift_set = lift_tree[pos - 1]

                    match_interval = list(lift_set)[0] if len(lift_set) == 1 else None

                    if not match_interval or match_interval.end != pos:
                        lift_coord_list.append(None)

                        raise RuntimeError(
                            (
                                'Found no matches in a lift-tree for a record within a '
                                'global to-reference tree: {}:{} (index={}, gap={})'
                            ).format(query_id, pos_org, index, gap)
                        )

                else:
                    raise RuntimeError(
                        (
                            'Found multiple matches in a lift-tree for a record within a '
                            'global to-reference tree: {}:{} (index={}, gap={})'
                        ).format(query_id, pos_org, index,  gap)
                    )

                # Interpolate coordinates
                if match_interval.data[1] - match_interval.data[0] > 1:
                    lift_pos = match_interval.data[0] + (pos - match_interval.begin)

                    match_list.append((
                        row['#CHROM'],
                        lift_pos,
                        row['REV'],
                        lift_pos,
                        lift_pos,
                        row['INDEX'],
                        row['SCORE']
                    ))

                else:  # Lift from missing bases on the target (insertion or deletion)
                    match_list.append((
                        row['#CHROM'],
                        match_interval.data[1],
                        row['REV'],
                        match_interval.data[1],
                        match_interval.data[1],
                        row['INDEX'],
                        row['SCORE']
                    ))

            # Append
            lift_coord_list.append(match_list)

        # Return coordinates
        if ret_list:
            return lift_coord_list
        else:
            return lift_coord_list[0]

    def lift_to_qry(self, reference_id, coord):
        """
        Lift coordinates from reference to query.

        For each element in `coord`, returns a list of tuples in the form:
            0: Query ID.
            1: Position.
            2: Is reverse complemented. `True` if the alignment was lifted through a reverse-complemented
                alignment record.
            3: Minimum position.
            4: Maximum position.
            5: Align record index.
            6: Alignment record score.

        If `coord` is a `list` or `tuple`, then a list for each position (where each element is a list of tuples, one
        for each alignment intersect). If `coord` is not a `list` or `tuple`, then there is no outer list and the
        returned results is a list of tuples for that position. (i.e. same as `lift_to_qry('chrom', (x,))[0]`). For
        each `coord` element that has no elements, `None` is returned instead of an empty list (a single `None` if
        `coord` is a single position not in a `list` or `tuple` and has no lift).

        :param reference_id: Subject ID.
        :param coord: Subject coordinates. May be a single value, list, or tuple.

        :return: Single coordinate tuple or list of coordinate tuples if `coord` is a list or tuple. Returns `None`
            for coordinates that cannot be lifted.
        """

        # Determine type
        if isinstance(coord, (list, tuple)):
            ret_list = True
        else:
            ret_list = False
            coord = (coord,)

        # Do lift
        lift_coord_list = list()

        for pos in coord:
            match_set = self.ref_tree[reference_id][pos:(pos + 1)]

            # Check coordinates
            if len(match_set) == 0:
                lift_coord_list.append(None)
                continue

            # Get list of lift matches
            match_list = list()

            for match_element in match_set:

                # Get lift tree
                index = match_element.data

                if index not in self.ref_cache.keys():
                    self._add_align(index)

                lift_tree = self.ref_cache[index]

                # Save row
                row = self.df.loc[index]

                # Get match record
                lift_set = lift_tree[pos:(pos + 1)]

                if len(lift_set) != 1:
                    raise RuntimeError(
                        (
                            'Program bug: Found no matches in a lift-tree for a record withing a '
                            'global to-query tree: {}:{} (index={})'
                        ).format(reference_id, pos, index)
                    )

                lift_interval = list(lift_set)[0]

                # Interpolate coordinates
                if lift_interval.data[1] - lift_interval.data[0] > 1:
                    qry_pos = lift_interval.data[0] + (pos - lift_interval.begin)

                else:  # Lift from missing bases on the target (insertion or deletion)
                    qry_pos = lift_interval.data[1]

                if row['REV']:
                    qry_pos = self.df_fai[row['QRY_ID']] - qry_pos

                match_list.append((
                    row['QRY_ID'],
                    qry_pos,
                    row['REV'],
                    qry_pos,
                    qry_pos,
                    row['INDEX'],
                    row['SCORE']
                ))

            # Add records for this position
            lift_coord_list.append(match_list)

        # Return coordinates
        if ret_list:
            return lift_coord_list
        else:
            return lift_coord_list[0]

    def lift_region_to_ref(self, region, gap=False, check_orientation=True):
        """
        Lift region to reference.

        :param region: Query region.
        :param gap: Interpolate within gap if `True`.

        :return: Subject region or `None` if it could not be lifted.
        """

        # Lift
        ref_pos, ref_end = self.lift_to_ref(region.chrom, (region.pos, region.end), gap)

        # Check lift: Must lift both ends to the same reference ID
        if ref_pos is None or ref_end is None:
            return None

        # Single position
        if len(ref_pos) != 1 or len(ref_end) != 1:
            return None

        ref_pos = ref_pos[0]
        ref_end = ref_end[0]

        if ref_pos[0] != ref_end[0]:
            return None

        # Check orientation
        is_rev = None

        if ref_pos[2] is not None and ref_end[2] is not None:

            if ref_pos[2] == ref_end[2]:
                is_rev = ref_pos[2]
            else:
                if check_orientation:
                    return None

        # Return region
        return pavlib.seq.Region(
            ref_pos[0], ref_pos[1], ref_end[1],
            is_rev=is_rev,
            pos_min=ref_pos[3], pos_max=ref_pos[4],
            end_min=ref_end[3], end_max=ref_end[4],
            pos_aln_index=(ref_pos[5],),
            end_aln_index=(ref_end[5],)
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
        if query_pos is None or query_end is None or len(query_pos) != 1 or len(query_end) != 1:
            return None

        query_pos = query_pos[0]
        query_end = query_end[0]

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

    def _get_ref_gap(self, query_id, pos):
        """
        Interpolate lift coordinates to an alignment gap.

        :param pos: Position on the query.

        :return: A tuple of (reference, pos, rev, min, max).
        """

        # Check arguments
        if pos is None:
            return None

        # Get alignment records for this query
        subdf = self.df.loc[self.df['QRY_ID'] == query_id]

        # Must be flanked by two query on either side
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

        # Rows must be mapped to the same reference
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
        ref_bp = row['POS']
        qry_bp = 0

        itree_ref = intervaltree.IntervalTree()
        itree_qry = intervaltree.IntervalTree()

        # Get CIGAR and check query start position
        cigar_op_list = list(pavlib.align.util.cigar_str_to_tuples(self.df.loc[index]))

        clipped_bp = 0

        cigar_index = 0

        while cigar_index < len(cigar_op_list) and cigar_op_list[cigar_index][1] in {'S', 'H'}:
            clipped_bp += cigar_op_list[cigar_index][0]
            cigar_index += 1

        # Build trees
        for cigar_len, cigar_op in cigar_op_list:

            if cigar_op in {'=', 'X', 'M'}:

                itree_ref[ref_bp:(ref_bp + cigar_len)] = (qry_bp, qry_bp + cigar_len)
                itree_qry[qry_bp:(qry_bp + cigar_len)] = (ref_bp, ref_bp + cigar_len)

                ref_bp += cigar_len
                qry_bp += cigar_len

            elif cigar_op == 'I':

                itree_qry[qry_bp:(qry_bp + cigar_len)] = (ref_bp, ref_bp + 1)

                qry_bp += cigar_len

            elif cigar_op == 'D':

                itree_ref[ref_bp:(ref_bp + cigar_len)] = (qry_bp, qry_bp + 1)

                ref_bp += cigar_len

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
        self.qry_cache[index] = itree_qry

        # Add index to end of queue
        self.cache_queue.appendleft(index)

    def _check_and_clear(self):
        """
        Check alignment cache and clear if necessary to make space.
        """

        while len(self.cache_queue) >= self.cache_align:
            index = self.cache_queue.pop()

            del(self.ref_cache[index])
            del(self.qry_cache[index])
