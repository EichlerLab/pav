"""
Routines for aligned contigs.
"""

import Bio.SeqIO
import collections
import numpy as np
import pandas as pd
import pysam
import re

import kanapy


class Region:
    """
    Represents a region (chromosome, pos, and end) in 0-based half-open coordinates (BED).

    Tracks orientation as `is_rev` for other code to use (not used by this object). By default, it assumes that `is_rev`
    is `True` if `pos` > `end`.

    If breakpoints are uncertain, minimum and maximum values can be set for pos and end.

    If the region is associated with alignment records, the alignment record index for pos and end can be also be
    tracked with this object.
    """

    def __init__(
            self,
            chrom, pos, end,
            is_rev=None,
            pos_min=None, pos_max=None,
            end_min=None, end_max=None,
            pos_aln_index=None,
            end_aln_index=None
    ):
        self.chrom = chrom
        self.pos = int(pos)
        self.end = int(end)

        self.pos_min = self.pos if pos_min is None else int(pos_min)
        self.pos_max = self.pos if pos_max is None else int(pos_max)
        self.end_min = self.end if end_min is None else int(end_min)
        self.end_max = self.end if end_max is None else int(end_max)

        self.pos_aln_index = pos_aln_index
        self.end_aln_index = end_aln_index

        # Swap and set orientation to reverse if pos and end are reversed
        if self.pos > self.end:
            tmp = self.pos
            self.pos = self.end
            self.end = tmp

            self.end_min = self.pos if pos_min is None else int(pos_min)
            self.end_max = self.pos if pos_max is None else int(pos_max)
            self.pos_min = self.end if end_min is None else int(end_min)
            self.pos_max = self.end if end_max is None else int(end_max)

            tmp = self.pos_aln_index
            self.pos_aln_index = self.end_aln_index
            self.end_aln_index = tmp

            if is_rev is None:
                is_rev = True

        if is_rev is None:
            is_rev = False

        self.is_rev = is_rev

    def __repr__(self):
        """
        :return: Coordinate string in 0-based half-open coordinates (BED).
        """
        return '{}:{}-{}'.format(self.chrom, self.pos, self.end)

    def to_base1_string(self):
        """
        :return: Coordinate string in 1-based closed coordinates (Samtools, UCSC browser).
        """
        return '{}:{}-{}'.format(self.chrom, self.pos + 1, self.end)

    def __len__(self):
        """
        Get length of region.

        :return: Length of region.
        """
        return self.end - self.pos

    def region_id(self):
        """
        Get a region ID using SV-Pop notation for variants and regions.

        :return: Region ID.
        """
        return '{}-{}-RGN-{}'.format(self.chrom, self.pos, self.end - self.pos)

    def expand(self, expand_bp, min_pos=0, max_end=None, shift=True, balance=0.5):
        """
        Expand this region by `expand_bp`.

        Do not expand beyond `min_pos` and `max_end`, if defined. If position limits are reached and `shift` is `True`,
        then set the boundary at the limit and expand in the free position up to its limit, if defined. The region is
        never expanded beyond `expand_bp`, but if limits are reached, it may expand less.

        Expand can be asymetrical if `balance` is defined and not 0.5 (must be between 0 and 1). `pos` is expanded by
        `int(expand_bp * balance)`, and `end` is expanded by `int(expand_bp * (1 - balance))`.

        :param expand_bp: Expand this number of bases. May be negative (untested).
        :param min_pos: Lower limit on the new position.
        :param max_end: Upper limit on the new end.
        :param shift: If `True`, shift the expand window if one end exceeds a limit. Tries to preserve `expand_bp` while
            honoring the limits.
        :param balance: Shift `pos` by `int(expand_bp * balance)` and `end` by `int(expand_bp * (1 - balance))`. If
            `0.5`, both positions are shifted equally.
        """

        # Set balance
        if balance is None:
            balance = 0.5

        try:
            if not (0 <= balance <= 1):
                raise RuntimeError('balance must be in range [0, 1]: {}'.format(balance))

        except ValueError:
            raise RuntimeError('balance is not numeric: {}'.format(balance))

        # Get new positions
        expand_pos = int(expand_bp * balance)
        expand_end = np.max([0, expand_bp - expand_pos])

        new_pos = int(self.pos - expand_pos)
        new_end = int(self.end + expand_end)

        # Shift pos
        if min_pos is not None and new_pos < min_pos:
            min_diff = min_pos - new_pos

            if shift:
                new_end += min_diff

            new_pos = min_pos

        # Shift end
        if max_end is not None:
            if max_end.__class__ == pd.core.series.Series and self.chrom in max_end.index:
                max_end = max_end[self.chrom]
            else:
                max_end = None

        if max_end is not None and new_end > max_end:
            max_diff = new_end - max_end

            if shift:
                new_pos -= max_diff

                if new_pos < min_pos:
                    new_pos = min_pos

            new_end = max_end

        # Check for over-contraction (if expand_bp was negative)
        if new_end < new_pos:
            new_end = new_pos = (new_end + new_pos) // 2

        # Assign new coordinates
        self.pos = new_pos
        self.end = new_end

        self.pos_min = self.pos
        self.pos_max = self.pos
        self.end_min = self.end
        self.end_max = self.end

    def __getitem__(self, key):
        """
        Handle `self[key]`.

        :param key: Key (chrom, pos, or end).

        :return: Value.
        """

        if key not in {'chrom', 'pos', 'pos1', 'end'}:
            raise IndexError('No key in Region: {}'.format(key))

        if key == 'pos1':
            return self.pos + 1

        return self.__dict__[key]

    def __eq__(self, other):
        """
        Handle self == other.

        :param other: Other object.

        :return: `True` if `other` is a `Region` with the same chromosome and coordinates.
        """
        if not issubclass(other.__class__, Region):
            return False

        return (
            self.chrom == other.chrom and
            self.pos == other.pos and
            self.end == other.end
        )

    def __add__(self, other):

        if not np.issubdtype(other.__class__, np.integer):
            raise RuntimeError('Region addition expected integer, receieved: {:s}'.format(other.__class__.__name__))

        self.pos += other
        self.end += other

    def __sub__(self, other):

        if not np.issubdtype(other.__class__, np.integer):
            raise RuntimeError('Region addition expected integer, receieved: {:s}'.format(other.__class__.__name__))

        self.pos -= other
        self.end -= other

    def copy(self):
        """
        Deep copy this object.

        :return: New copy.
        """
        return Region(
            self.chrom, self.pos, self.end, self.is_rev, self.pos_min, self.pos_max, self.end_min, self.end_max
        )


def region_from_string(rgn_str, is_rev=None, base0half=False):
    """
    Get a region object from a string (e.g. "chr1:123456-234567").

    :param rgn_str: Region string.
    :param is_rev: Region is reverse-complemented if `True`. May be used by some code when extracting sequence
        (e.g. inversion caller). If `None`, determine automatically if positions are reversed (pos > end).
    :param base0half: If `False` (default), then string is in base-1 closed coordinates (the first base in the
        chromosome is 1, and the positions are inclusive). If `True`, expect BED coordinates.

    :return: Region object.
    """

    match_obj = re.match(r'^([^:]+):(\d+)-(\d+)$', rgn_str)

    if match_obj is None:
        raise RuntimeError('Region is not in expected format (chrom:pos-end): {}'.format(rgn_str))

    pos = int(match_obj[2])
    end = int(match_obj[3])

    if not base0half:
        pos -= 1

    return Region(match_obj[1], pos, end, is_rev=is_rev)


def region_from_id(region_id):
    """
    Take a region ID (CHROM-POS-SVTYPE-LEN) and translate to a Region object. Typically used to translate "RGN" IDs.

    :param region_id: Region ID.

    :return: Region object.
    """

    tok = region_id.split('-')

    if len(tok) != 4:
        raise RuntimeError('Unrecognized region ID: {}'.format(region_id))

    return Region(tok[0], int(tok[1]) - 1, int(tok[1]) - 1 + int(tok[3]))


def ref_kmers(region, fa_file_name, k_util):

    ref_seq = region_seq_fasta(region, fa_file_name, False)

    ### Get reference k-mer counts ###

    ref_mer_count = collections.Counter()

    for kmer in kanapy.util.kmer.stream(ref_seq, k_util):
        ref_mer_count[kmer] += 1

    return ref_mer_count


def region_seq_fasta(region, fa_file_name, rev_compl=None):
    """
    Get sequence from an indexed FASTA file. FASTA must have ".fai" index.

    :param region: Region object.
    :param fa_file_name: FASTA file name.
    :param rev_compl: Reverse-complement sequence is `True`. If `None`, reverse-complement if `region.is_rev`.

    :return: String sequence.
    """

    with pysam.FastaFile(fa_file_name) as fa_file:
        sequence = fa_file.fetch(region.chrom, region.pos, region.end)

        if rev_compl is None:
            if region.is_rev:
                return str(Bio.Seq.Seq(sequence).reverse_complement())
        else:
            if rev_compl:
                return str(Bio.Seq.Seq(sequence).reverse_complement())

        return sequence


# Some pysam implementations do not d o well with repeated calls against an open pysam.FastaFile, so use
# region_seq_fasta() for all calls.
# def region_seq_pysam(region, fa_file, rev_compl=None):
#     """
#     Get a region sequence from an open pysam file.
#
#     :param region: Region object.
#     :param fa_file: Open pysam FastaFile
#     :param rev_compl: Reverse-complement sequence is `True`. If `None`, reverse-complement if `region.is_rev`.
#
#     :return: Sequence from `region`
#     """
#
#     sequence = fa_file.fetch(region.chrom, region.pos, region.end)
#
#     if rev_compl is None:
#         if region.is_rev:
#             return str(Bio.Seq.Seq(sequence).reverse_complement())
#     else:
#         if rev_compl:
#             return str(Bio.Seq.Seq(sequence).reverse_complement())
#
#     return sequence
