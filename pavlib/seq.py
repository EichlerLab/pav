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

    :param region: Region object to extract a region, or a string with the recrord ID to extract a whole record.
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


def expand_input(file_name_list):
    """
    Expand input to a list of tuples:
    * [0]: File name.
    * [1]: File type ("fasta", "fastq", or "gfa")

    File type does not change if the file is gzipped, the downstream input function will resolve that.

    This function traverses FOFN files (file of file names) recursively until FASTA, FASTQ, or GFA files are found.

    :param file_name_list: List of input file names.
    """

    # Check arguments
    if file_name_list is None:
        raise RuntimeError('Cannot create input FASTA: Input name list is None')

    # Check files
    if issubclass(file_name_list.__class__, str):
        file_name_list = [file_name_list]

    elif issubclass(file_name_list.__class__, tuple):
        file_name_list = list(file_name_list)

    elif issubclass(file_name_list.__class__, set):
        file_name_list = sorted(file_name_list)

    elif not issubclass(file_name_list.__class__, list):
        raise RuntimeError(f'Unrecognized type for input file name list: {file_name_list.__class__}')


    # Generate a list of files traversing into FOFN files
    file_name_tuples = list()
    fofn_set = list()  # Set of visited FOFN files (prevent recursive traversal)

    while len(file_name_list) > 0:

        # Next file name
        file_name = file_name_list[0].strip()
        file_name_list = file_name_list[1:]

        if not file_name:
            continue

        # Get file extension
        file_name_lower = file_name.lower()

        if file_name_lower.endswith('.gz'):  # Strip GZ, the downstream input functions will detect file type
            file_name_lower = file_name_lower.rsplit('.', 1)[0]

        if '.' not in file_name_lower:
            raise RuntimeError(f'No recognizable extension in file name: {file_name}')

        file_name_ext = file_name_lower.rsplit('.', 1)[1]

        # Expand FOFN files
        if file_name_ext == 'fofn':

            # Check for recursive FOFN traversal
            file_name_real = os.path.realpat(file_name)

            if file_name_real in fofn_set:
                raise RuntimeWarning(f'Detected recursive FOFN traversal, ignoring redundant entry: {file_name}')
                continue

            fofn_set.add(file_name_real)

            # Append FOFN entries to the input file list
            with svpoplib.seq.PlainOrGzReader(file_name) as in_file:
                for line in in_file:
                    line = line.strip()

                    if line:
                        file_name_list.append(line.lower())

        # FASTA files
        if file_name_ext in {'fasta', 'fa', 'fn'}:
            file_name_tuples.append((file_name, 'fasta'))

        elif file_name_ext in {'fastq', 'fq'}:
            file_name_tuples.append((file_name, 'fastq'))

        elif file_name_ext == 'gfa':
            file_name_tuples.append((file_name, 'gfa'))

        else:
            raise RuntimeError(f'Unrecognized file extension {file_name_ext}: {file_name}')

        # Return tuples
        return file_name_tuples


def input_tuples_to_fasta(file_name_tuples, out_file_name):
    """
    Convert a list of input files to a single FASTA entry. Input files may be FASTA, FASTQ, or GFA.

    :param file_name_tuples: List of tuples for each input entry ([0]: File name, [1]: File format). The file format
        must be "fasta", "fastq", or "gfa" (case sensitive).
    :param out_file_name: Output file. If `file_name_tuples` or contains only empty files, the output file is also
        empty (0 bytes) indicating to PAV that data for this haplotype are missing.
    """

    # Check input files, fail early
    has_data = False

    if file_name_tuples is None:
        file_name_tuples = []

    for index in range(len(file_name_tuples)):
        file_name, file_format = file_name_tuples[index]

        if file_format not in {'fasta', 'fastq', 'gfa'}:
            raise RuntimeError(f'Unrecognized file format "{file_format}": {file_name}')

        if not os.path.isfile(file_name):
            raise RuntimeError(f'Input file does not exist or is not a regular file: {file_name}')

        if os.stat(file_name).st_size > 0:
            has_data = True
        else:
            file_name_tuples[index][1] = 'empty'  # Mark file as empty

    # Stop if there are no records. An empty file to signal downstream steps that there is no data for this haplotype.
    if not has_data:
        with open(out_file_name, 'w') as out_file:
            pass

        return

    # Record iterator
    def input_record_iter():

        record_id_set = set()

        for file_name, file_format in file_name_tuples:

            if file_format in {'fasta', 'fastq'}:
                for record in svpoplib.seq.fa_to_record_iter(file_name, input_format=file_format):

                    if record.id in record_id_set:
                        raise RuntimeError(f'Duplicate record ID in input: {record_id}')

                    record_id_set.add(record.id)

                    yield record

            elif file_format == 'gfa':
                for record in svpoplib.seq.gfa_to_record_iter(file_name):

                    if record.id in record_id_set:
                        raise RuntimeError(f'Duplicate record ID in input: {record_id}')

                    record_id_set.add(record.id)

                    yield record

            elif file_format not in {'skip', 'empty'}:
                raise RuntimeError(f'Program bug: Unrecognized file type "{file_format}" after checking input.')

    # Traverse entries
    with Bio.bgzf.open(out_file_name, 'wb') as out_file:
        Bio.SeqIO.write(input_record_iter(), out_file, 'fasta')

    return


# def copy_fa_to_gz(in_file_name, out_file_name):
#     """
#     Copy FASTA file to gzipped FASTA if the file is not already gzipped. If the input file is empty, then write an
#     empty file.
#
#     :param in_file_name: Input FASTA file.
#     :param out_file_name: Output gzipped FASTA file.
#     """
#
#     # Write empty files in FASTA is empty
#     if os.stat(in_file_name).st_size == 0:
#         with open(out_file_name, 'w') as out_file:
#             pass
#
#         return
#
#     # Determine if file is BGZF compressed
#     is_bgzf = False
#
#     try:
#         with Bio.bgzf.open(in_file_name, 'r') as in_file_test:
#             is_bgzf = True
#
#     except ValueError:
#         pass
#
#     # Copy or compress
#     if is_bgzf:
#
#         # Copy file if already compressed
#         shutil.copyfile(in_file_name, out_file_name)
#
#     else:
#         # Compress to BGZF
#
#         is_gz = False
#
#         try:
#             with gzip.open(in_file_name, 'r') as in_file_test:
#
#                 line = next(in_file_test)
#
#                 is_gz = True
#
#         except OSError:
#             pass
#
#         if is_gz:
#             # Re-compress to BGZF
#
#             with gzip.open(in_file_name, 'rb') as in_file:
#                 with Bio.bgzf.open(out_file_name, 'wb') as out_file:
#                     for line in in_file:
#                         out_file.write(line)
#
#         else:
#             # Compress plain text
#
#             with open(in_file_name, 'r') as in_file:
#                 with Bio.bgzf.open(out_file_name, 'wb') as out_file:
#                     for line in in_file:
#                         out_file.write(line)

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
