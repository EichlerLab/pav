"""
Routines for aligned contigs.
"""

import Bio.SeqIO
import Bio.Seq
import collections
import gzip
import numpy as np
import os
import pandas as pd
import pysam
import re
import shutil

import svpoplib
import kanapy

import pavlib


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
        self.chrom = str(chrom)
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
        :return: Coordinate string in 1-based closed coordinates (Samtools, UCSC browser).
        """

        return self.to_base1_string()

    def to_base1_string(self):
        """
        :return: Coordinate string in 1-based closed coordinates (Samtools, UCSC browser).
        """
        return '{}:{}-{}'.format(self.chrom, self.pos + 1, self.end)

    def to_bed_string(self):
        """
        :return: Coordinate string in 0-based closed coordinates (BED).
        """
        return '{}\t{}\t{}'.format(self.chrom, self.pos + 1, self.end)

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

        return (
            self.chrom == other.chrom and
            self.pos == other.pos and
            self.end == other.end
        )

    def __lt__(self, other):
        """
        Determine if this object is less than other.

        :param other: Other object.

        :return: `True` if this object is less than `other`.
        """

        return (self.chrom, self.pos, self.end) < (other.chrom, other.pos, other.end)

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

    def contains(self, other):
        if other is None:
            return False

        return other.chrom == self.chrom and other.pos >= self.pos and other.end <= self.end

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

    rgn_str_no_comma = rgn_str.replace(',', '')
    match_obj = re.match(r'^([^:]+):(\d+)-(\d+)$', rgn_str_no_comma)

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
    """
    Get a counter keyed by k-mers.

    :param region: Region to extract sequence from.
    :param fa_file_name: FASTA file to extract sequence from.
    :param k_util: K-mer utility for k-merizing sequence.

    :return: A collections.Counter object key k-mer keys and counts.
    """

    ref_seq = region_seq_fasta(region, fa_file_name, False)

    ### Get reference k-mer counts ###

    ref_mer_count = collections.Counter()

    for kmer in kanapy.util.kmer.stream(ref_seq, k_util):
        ref_mer_count[kmer] += 1

    return ref_mer_count


def region_seq_fasta(region, fa_file_name, rev_compl=None):
    """
    Get sequence from an indexed FASTA file. FASTA must have ".fai" index.

    :param region: Region object to extract a region, or a string with the record ID to extract a whole record.
    :param fa_file_name: FASTA file name.
    :param rev_compl: Reverse-complement sequence is `True`. If `None`, reverse-complement if `region.is_rev`.

    :return: String sequence.
    """

    with pysam.FastaFile(fa_file_name) as fa_file:

        if region.__class__ == str:
            is_region = False
        elif region.__class__ == Region:
            is_region = True
        else:
            raise RuntimeError('Unrecognized region type: {}: Expected Region (pavlib.seq) or str'.format(str(region.__class__.__name__)))

        if is_region:
            sequence = fa_file.fetch(region.chrom, region.pos, region.end)
        else:
            sequence = fa_file.fetch(region, None, None)

        if rev_compl is None:
            if is_region and region.is_rev:
                return str(Bio.Seq.Seq(sequence).reverse_complement())
        else:
            if rev_compl:
                return str(Bio.Seq.Seq(sequence).reverse_complement())

        return sequence

def variant_seq_from_region(df, fa_file_name, region_col='QRY_REGION', strand_col='QRY_STRAND', id_col='ID', seq_upper=False):
    """
    Get sequence from an indexed FASTA file. FASTA must have ".fai" index.

    :param df: Variant DataFrame.
    :param fa_file_name: FASTA file name.
    :param region_col: Region column name. If column is "-" or `True`, the sequence is reverse complemented. Must be
        "+", "-", `True`, or `False`. The boolean values allow "IS_REV" column to be used where the sequence is on
        the negative strand if `True`. If `None`, no reverse complementing is performed.
    :param strand_col: Strand column name.
    :param id_col: ID column name.
    :param seq_upper: Upper-case sequences if `True` (removes soft-mask annotations in the sequence string).

    :return: Iterator over Bio.Seq sequence objects with sequences in reference orientation.
    """

    if df is None:
        raise RuntimeError('Variant DataFrame is missing')

    fa_file_name = fa_file_name.strip() if fa_file_name is not None else None

    if fa_file_name is None or len(fa_file_name) == 0:
        raise RuntimeError('FASTA file name is missing')

    if not os.path.isfile(fa_file_name):
        raise RuntimeError(f'FASTA file does not exist or is not a regular file: {fa_file_name}')

    if region_col not in df.columns:
        raise RuntimeError(f'Region column not found in DataFrame: "{region_col}"')

    if strand_col is not None and strand_col not in df.columns:
        raise RuntimeError(f'Strand column not found in DataFrame: "{strand_col}"')

    if id_col not in  df.columns:
        raise RuntimeError(f'ID column not found in DataFrame: "{id_col}"')

    n_blank = np.sum(df[id_col].apply(len) == 0)

    if n_blank > 0:
        raise RuntimeError(f'Found {n_blank} empty IDs in DataFrame: "{id_col}"')

    dup_ids = sorted([id for id, count in collections.Counter(df[id_col]).items() if count > 1])

    if dup_ids:
        n_dup = len(dup_ids)
        dup_ids = ', '.join(dup_ids[3:]) + (', ...' if n_dup > 3 else '')
        raise RuntimeError(f'Found {n_dup} duplicate IDs in DataFrame: "{id_col}": {dup_ids}')

    with pavlib.io.FastaReader(fa_file_name) as fa_file:

        for index, row in df.iterrows():

            region = region_from_string(row[region_col])

            seq = Bio.Seq.Seq(fa_file.fetch(region.chrom, region.pos, region.end))

            if seq_upper:
                seq = seq.upper()

            if strand_col is not None:
                if row[strand_col] in {'-', True}:
                    seq = seq.reverse_complement()
                elif row[strand_col] not in {'+', True}:
                    raise RuntimeError(f'Unrecognized strand in record {row[id_col]}: "{row[strand_col]}"')

            yield Bio.SeqRecord.SeqRecord(
                seq, id=row[id_col], description=''
            )
