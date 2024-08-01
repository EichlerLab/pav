"""
Alignment routines.
"""

import numpy as np

import pavlib


class AffineScoreModel:
    """
    Affine score model with default values modeled on minimap2 (2.26) default parameters.

    :param match: Match score [2].
    :param mismatch: Mismatch penalty [4].
    :param affine_gap: A list of tuples of two elements (gap open, gap extend) penalty.
    :param template_switch: Template switch penalty. If none, defaults to 2x the penalty of a 50 bp gap.
    """

    def __init__(self, match=2, mismatch=4, affine_gap=((4, 2), (24, 1)), template_switch=None):
        self.match = np.abs(match)
        self.mismatch = np.abs(mismatch)
        self.affine_gap = (
            (np.abs(gap_open), np.abs(gap_extend))
                for gap_open, gap_extend in affine_gap
        )

        if template_switch is None:
            self.template_switch = - np.abs(2 * self.gap(50))
        else:
            self.template_switch = - np.abs(template_switch)

    def match(self, n=1):
        """
        Score match.

        :param n: Number of matching bases.
        """

        return self.match * n

    def mismatch(self, n=1):
        """
        Score mismatch.

        :param n: Number of mismatched bases.
        """

        return -self.mismatch * n

    def gap(self, n=1):
        """
        Score gap (insertion or deletion). Compute all affine aligment gap scores and return the lowest.

        :param n: Size of gap.
        """

        if n == 0:
            return 0

        return - np.min(
            [
                gap_open + (gap_extend * n)
                    for gap_open, gap_extend in self.affine_gap
            ]
        )

    def score(self, op_len, op_code):
        """
        Score a CIGAR operation and return 0 if the operation is not scored (S, H, N, and P CIGAR operations).

        :param op_len: CIGAR operation length.
        :param op_code: CIGAR operation code, must be capitalized.
        """

        if op_code == '=':
            return self.match * op_len

        elif op_code == 'X':
            return -self.mismatch * op_len

        elif op_code in {'I', 'D'}:
            return -self.gap(op_len)

        elif op_code in {'S', 'H'}:
            return 0

        elif op_code == 'M':
            raise RuntimeError('Cannot score alignments with match ("M") in CIGAR string (requires "=" and "X")')

        elif op_code in {'N', 'P'}:
            return 0

        else:
            raise RuntimeError('Unrecognized CIGAR op code: {op_code}')


def score_align(row, score_model=None):
    """
    Add a SCORE column based on a scoring model for alignment events.

    :param row: Alignment row.
    :param score_model: Alignment model or None to use the default affine-gap model.

    :return: Alignment score.
    """

    if score_model is None:
        score_model = AffineScoreModel()

    return sum([
        score_model.score(op_len, op_code)
            for op_len, op_code in pavlib.align.util.cigar_str_to_tuples(row)
    ])
