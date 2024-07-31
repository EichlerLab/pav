"""
Alignment chainging functions
"""

import collections


class AnchorChainNode:
    """
    Describes an SV (simple or complex) anchored by aligned contig bases on each side (start_index and end_index).

    :param start_index: Index of the left-most aligned segment in contig coordinate order.
    :param start_index: Index of the right-most aligned segment in contig coordinate order.
    """

    def __init__(self, start_index, end_index):
        self.start_index = start_index
        self.end_index = end_index

        self.next_node_list = list()

    def __repr__(self):
        """
        :return: String representation of a choin object.
        """

        return f'AnchorChainNode(start={self.start_index}, end={self.end_index}, children={self.n_children()})'

    def n_children(self):
        """
        Get the number of child nodes from this node.

        :return: Number if child nodes.
        """
        return len(self.next_node_list)


class AnchorChainContainer:
    """
    A container of anchor chains. Allows chains to be constructed
    """

    def __init__(self):

        # Chain end list. Key: end_index. Value: list of anchor chain nodes with "end_index" equal to the key index.
        # Needed for connecting new nodes to the chain graphs ending on new node start indexes.
        self.chain_end_list = collections.defaultdict(list)

        # List of nodes with no parents.
        self.source_node_list = list()

        # A dictionary of index of all anchor chains. Key: (start_index, end_index). Value: Anchor chain.
        self.chain_dict = dict()

    def add_anchor(self, start_index, end_index):
        """
        Add a pair of anchors to this graph.

        :param start_index: First aligned segment.
        :param end_index: Last aligned segment.
        """

        if start_index >= end_index:
            raise RuntimeError(f'end_index ({end_index}) must be greater that start_index ({start_index})')

        # Skip if this anchor already exists
        if (start_index, end_index) in self.chain_dict:
            return

        # Create node
        new_node = AnchorChainNode(start_index, end_index)

        self.chain_dict[(start_index, end_index)] = new_node

        # Find parent anchors
        parent_node_list = self.chain_end_list[start_index]

        if len(parent_node_list) > 0:
            for parent_node in parent_node_list:
                parent_node.next_node_list.append(new_node)

        else:
            self.source_node_list.append(new_node)

        # Append to chain end list
        self.chain_end_list[end_index].append(new_node)


def can_anchor(row_a, row_b, score_model, min_score=1000, gap_scale=1):
    """
    Determine if two alignment rows can anchor a rearrangement. Requires the "SCORE" column is added to the alignment
    rows.

    Both rows are records from an alignment BED file representing a candidate poir of alignments for anchoring an SV
    (simple or complex) between them. If either row is `None`, they are not aligned to the same reference contig, or
    if they are not aligned in the same orientation, `False` is returned.

    Each row has an alignment score ("SCORE" column), which is computed if the "SCORE" column is absent using
    `score_model`. If either anchor's score is less than `min_score`, `False` is returned. The minimum value of the pair
    alignment scores (from row_a and row_b) is compared to the gap between them. The gap between the alignments is found
    by adding the number of contig bases skipped between the alignment records and the number of reference bases skipped
    moving from `row_a` to `row_b` (orientation does not matter, may skip forward (DEL) or backward (DUP)).

    The gap is scored as one gap with `score_model`. If the gap penalty exceeds the minimum alignment score of the
    anchors, then `False` is returned. Otherwise, `True` is returned and these alignments can function as anchors for
    an SV between them.

    :param row_a: Row earlier in the alignment chain (contig position order).
    :param row_b: Row later in the alignment chain (contig position order).
    :param score_model: Alignment scoring model used to score the gap between two rows (ref gap and tig gap).
    :param min_score: Cannot be an anchor if either anchor's alignment score is less than this value.
    :param gap_scale: Scale gap score by this factor. A value of less than 1 reduces the gap penalty (e.g. 0.5 halves
        it), and a value greater than 1 increases the gap penalty (e.g. 2.0 doubles it).

    :return: `True` if both rows are not `None` and are collinear in contig and reference space.
    """

    # Both rows should be present and in the same orientation
    if row_a is None or row_b is None or row_a['#CHROM'] != row_b['#CHROM']:
        return False

    is_rev = row_a['REV']

    if is_rev != row_b['REV']:
        return False

    anchor_score = min([row_a['SCORE'], row_b['SCORE']])

    if anchor_score < min_score:
        return False

    # Check reference contiguity
    if is_rev:
        ref_l_end = row_b['END']
        ref_r_pos = row_a['POS']
    else:
        ref_l_end = row_a['END']
        ref_r_pos = row_b['POS']

    # Score gap
    gap_len = row_b['QRY_POS'] - row_a['QRY_END']

    if gap_len < 0:
        raise RuntimeError(f'Alignment rows are out of order: Negative distance {gap_len}: row_a index "{row_a.name}", row_b index "{row_b.name}"')

    gap_len += abs(ref_r_pos - ref_l_end)

    if gap_len == 0:
        return True

    return score_model.gap(gap_len) * gap_scale < anchor_score


def can_reach_anchor(row_l, row_r, score_model):
    """
    Determine if a left-most anchor can reach as far as a right-most anchor. This function only tells the traversal
    algorithm when to stop searching for further right-most anchor candidates, it does not determine if the two
    alignments are anchors (does not consider the score of the right-most alignment or the reference position or
    orientation).

    :param row_l: Left-most archor in contig coordinates.
    :param row_r: Right-most archor in contig coordinates.
    :param score_model: Model for determining a gap score.

    :return: `True` if the anchor in alignment `row_l` can reach as far as the start of `row_r` based on the alignment
        score of `row_l` and the distance in contig coordinates between the end of `row_l` and the start of `row_r`.
    """

    # Check rows
    if row_l is None or row_r is None:
        raise RuntimeError('Cannot score tig distance for records "None"')

    if row_l['QRY_ID'] != row_r['QRY_ID']:
        raise RuntimeError(f'Cannot score tig distance for mismatching queries: "{row_l["QRY_ID"]}" and "{row_r["QRY_ID"]}"')

    tig_dist = row_r['QRY_POS'] - row_l['QRY_END']

    if tig_dist < 0:
        raise RuntimeError(f'Cannot score tig distance for out-of-order (by tig coordinates): indexes "{row_l["INDEX"]}" and "{row_r["INDEX"]}"')

    # Cannot anchor if query distance is too large
    return tig_dist == 0 or row_l['SCORE'] - score_model.gap(tig_dist) > 0
