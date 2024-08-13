"""
Alignment routines.
"""

import abc
import numpy as np
import re

import pavlib

AFFINE_SCORE_MATCH = 2
AFFINE_SCORE_MISMATCH = 4
AFFINE_SCORE_GAP = ((4, 2), (24, 1))
AFFINE_SCORE_TS = None

DEFAULT_ALIGN_SCORE_MODEL = f'affine::match={AFFINE_SCORE_MATCH},mismatch={AFFINE_SCORE_MISMATCH},gap={";".join([f"{gap_open}:{gap_extend}" for gap_open, gap_extend in AFFINE_SCORE_GAP])}'

class ScoreModel(object, metaclass=abc.ABCMeta):
    """
    Score model interface.
    """

    @abc.abstractmethod
    def match(self, n=1):
        """
        Score match.

        :param n: Number of matching bases.
        """

        raise NotImplementedError

    @abc.abstractmethod
    def mismatch(self, n=1):
        """
        Score mismatch.

        :param n: Number of mismatched bases.
        """

        raise NotImplementedError

    @abc.abstractmethod
    def gap(self, n=1):
        """
        Score gap (insertion or deletion). Compute all affine aligment gap scores and return the lowest.

        :param n: Size of gap.
        """

        raise NotImplementedError

    @abc.abstractmethod
    def template_switch(self):
        """
        Score a template switch.

        :return: Template switch score.
        """

        raise NotImplementedError

    def score(self, op_len, op_code):
        """
        Score a CIGAR operation and return 0 if the operation is not scored (S, H, N, and P CIGAR operations). The
        CIGAR operations may be a capital letter or symbol (e.g. "=", "X", "I", etc.) on the numeric CIGAR operation
        code (defined in `pavlib.align.util`).

        :param op_len: CIGAR operation length.
        :param op_code: CIGAR operation code, must be capitalized.
        """

        if op_code in {'=', pavlib.align.util.CIGAR_EQ}:
            return self.match(op_len)

        elif op_code in {'X', pavlib.align.util.CIGAR_X}:
            return self.mismatch(op_len)

        elif op_code in {'I', 'D', pavlib.align.util.CIGAR_I, pavlib.align.util.CIGAR_D}:
            return self.gap(op_len)

        elif op_code in {'S', 'H', pavlib.align.util.CIGAR_S, pavlib.align.util.CIGAR_H}:
            return 0

        elif op_code in {'M', pavlib.align.util.CIGAR_M}:
            raise RuntimeError('Cannot score alignments with match ("M") in CIGAR string (requires "=" and "X")')

        elif op_code in {'N', 'P', pavlib.align.util.CIGAR_N, pavlib.align.util.CIGAR_P}:
            return 0

        else:
            raise RuntimeError(f'Unrecognized CIGAR op code: {op_code}')

    def score_cigar_tuples(self, cigar_tuples, rev=False):
        """
        Sum the scores for each operation over a list of CIGAR operations. The CIGAR operations may be a capital
        letter or symbol (e.g. "=", "X", "I", etc.) on the numeric CIGAR operation code (defined in
        `pavlib.align.util`).

        :param cigar_tuples: List of tuples of (op_len, op_code). If string, then assume this is a list of CIGAR
            operations and convert to tuples.
        :param rev: If True, each tuple is "(CIGAR op, CIGAR length)" (pysam) instead of the default order "(CIGAR
            length, CIGAR op)" (CIGAR-string order).

        :return: Sum of scores for each CIGAR operation.
        """

        if cigar_tuples is None:
            raise RuntimeError('CIGAR list is None')

        if isinstance(cigar_tuples, str):
            if len(cigar_tuples.strip()) == 0:
                raise RuntimeError('CIGAR string is empty')

            cigar_tuples = pavlib.align.util.cigar_str_to_tuples(cigar_tuples)

        if not rev:
            return np.sum([self.score(op_len, op_code) for op_len, op_code in cigar_tuples])
        else:
            return np.sum([self.score(op_len, op_code) for op_code, op_len in cigar_tuples])

    @abc.abstractmethod
    def __repr__(self):
        return 'ScoreModel(Interface)'


class AffineScoreModel(ScoreModel):
    """
    Affine score model with default values modeled on minimap2 (2.26) default parameters.

    :param match: Match score [2].
    :param mismatch: Mismatch penalty [4].
    :param affine_gap: A list of tuples of two elements (gap open, gap extend) penalty.
    :param template_switch: Template switch penalty. If none, defaults to 2x the penalty of a 50 bp gap.
    """

    def __init__(self, match=AFFINE_SCORE_MATCH, mismatch=AFFINE_SCORE_MISMATCH, affine_gap=AFFINE_SCORE_GAP, template_switch=AFFINE_SCORE_TS):
        self.score_match = np.abs(match)
        self.score_mismatch = -np.abs(mismatch)
        self.score_affine_gap = tuple((
            (-np.abs(gap_open), -np.abs(gap_extend))
                for gap_open, gap_extend in affine_gap
        ))

        if template_switch is None:
            self.score_template_switch = 2 * self.gap(50)
        else:
            try:
                self.score_template_switch = - np.abs(float(template_switch))
            except ValueError:
                raise ValueError(f'template_switch parameter is not numeric: {template_switch}')

    def match(self, n=1):
        """
        Score match.

        :param n: Number of matching bases.
        """

        return self.score_match * n

    def mismatch(self, n=1):
        """
        Score mismatch.

        :param n: Number of mismatched bases.
        """

        return self.score_mismatch * n

    def gap(self, n=1):
        """
        Score gap (insertion or deletion). Compute all affine alignment gap scores and return the lowest penalty.

        :param n: Size of gap.
        """

        if n == 0:
            return 0

        return np.max(
            [
                gap_open + (gap_extend * n)
                    for gap_open, gap_extend in self.score_affine_gap
            ]
        )

    def template_switch(self):
        """
        Score a template switch.

        :return: Template switch score.
        """
        return self.score_template_switch

    def __repr__(self):
        gap_str = ';'.join([f'{abs(gap_open)}:{abs(gap_extend)}' for gap_open, gap_extend in self.score_affine_gap])

        return f'AffineScoreModel(match={self.score_match},mismatch={-self.score_mismatch},gap={gap_str},ts={-self.score_template_switch})'


def get_score_model(param_string=None):
    """
    Get score model from a string of alignment parameters.

    :param param_string: Parameter string. May be None or an empty string (default model is used). If the string is
        an instance of `ScoreModel`, then the `ScoreModel` object is returned. Otherwise, the string is parsed and
        a score model object is returned.

    :return: A `ScoreModel` object.
    """

    if isinstance(param_string, ScoreModel):
        return param_string

    if param_string is not None:
        param_string = param_string.strip()

    if param_string is None or len(param_string) == 0:
        param_string = DEFAULT_ALIGN_SCORE_MODEL

    if '::' in param_string:
        model_type, model_params = re.split(r'::', param_string, maxsplit=1)

    else:
        model_type, model_params = 'affine', param_string

    if model_type == 'affine':
        return get_affine_by_params(model_params)

    raise RuntimeError(f'Unrecognized score model type: {model_type}')


def get_affine_by_params(param_string):
    """
    Parse a string to get alignment parameters from it.

    :param param_string: Parameter string.

    :return: A configured AffineScoreModel.
    """

    # Set defaults
    params = [
        AFFINE_SCORE_MATCH,
        AFFINE_SCORE_MISMATCH,
        AFFINE_SCORE_GAP,
        AFFINE_SCORE_TS
    ]

    keys = ['match', 'mismatch', 'gap', 'ts']

    # Sanitize parameter string
    if param_string is not None:
        param_string = param_string.strip()

        if len(param_string) == 0:
            param_string = None

    if param_string is None:
        param_string = DEFAULT_ALIGN_SCORE_MODEL

    # Parse param string
    param_pos = 0

    for tok in param_string.split(','):
        tok = tok.strip()

        if len(tok) == 0:
            param_pos += 1
            continue  # Accept default for missing

        if '=' in tok:
            param_pos = None  # Do not allow positional parameters after named ones

            key, val = tok.split('=')

            key = key.strip()
            val = val.strip()

        else:
            if param_pos is None:
                raise RuntimeError(f'Named parameters (with "=") must be specified after positional parameters (no "="): {param_string}')

            key = keys[param_pos]
            val = tok

            param_pos += 1

        if key in {'match', 'mismatch', 'ts'}:
            try:
                val = abs(float(val))
            except ValueError:
                raise RuntimeError(f'Unrecognized alignment parameter: {key} (allowed: match, mismatch, gap, ts)')

        if key == 'match':
            params[0] = val

        elif key == 'mismatch':
            params[1] = val

        elif key == 'gap':

            gap_list = list()

            for gap_pair in val.split(';'):
                gap_tok = gap_pair.split(':')

                if len(gap_tok) != 2:
                    raise RuntimeError(f'Unrecognized gap format: Expected "open:extend" for each element (multiple pairs separated by ";"): "{gap_pair}" in {val}')

                try:
                    gap_open = abs(float(gap_tok[0].strip()))
                    gap_extend = abs(float(gap_tok[1].strip()))

                except ValueError:
                    raise RuntimeError(f'Unrecognized gap format: Expected integer values for "open:extend" for each gap cost element: {val}')

                gap_list.append((gap_open, gap_extend))

            params[2] = tuple(gap_list)

        elif key == 'ts':
            params[3] = val

        else:
            raise RuntimeError(f'Unrecognized alignment parameter: {key} (allowed: match, mismatch, gap, ts)')

    # Return alignment object
    return AffineScoreModel(*params)
