"""
Routines for aligned contigs.
"""

import Bio.SeqIO
import collections
import intervaltree
import io
import numpy as np
import pandas as pd
import pysam
import re
import subprocess
import sys

import kanapy

class Region:
    """
    Represents a region (chromosome, pos, and end) in 0-based half-open coordinates (BED).
    """

    def __init__(self, chrom, pos, end):
        self.chrom = chrom
        self.pos = int(pos)
        self.end = int(end)

        if self.pos > self.end:
            tmp= self.pos
            self.pos = self.end
            self.end = tmp

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
        return Region(self.chrom, self.pos, self.end)


def region_from_string(rgn_str, base0half=False):
    """
    Get a region object from a string (e.g. "chr1:123456-234567").

    :param rgn_str: Region string.
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

    return Region(match_obj[1], pos, end)

def ref_kmers(region, fa_file_name, k_util):

    ref_seq = fa_region(region, fa_file_name)

    ### Get reference k-mer counts ###

    ref_mer_count = collections.Counter()

    for kmer in kanapy.util.kmer.stream(ref_seq, k_util):
        ref_mer_count[kmer] += 1

    return ref_mer_count

def tig_mer_stream(region, aln_file_name, subseq_exe, k_util, index=True):


    subseq_generator = subseq_region(region, aln_file_name, subseq_exe)

    # Get first aligned contig region
    try:
        sub_region, tig_seq = next(subseq_generator)

    except StopIteration:
        return None

    # Abort if there is more than one aligned contig
    try:
        sub_region, tig_seq = next(subseq_generator)

        return None  # Found another record, stop

    except StopIteration:
        pass  # Only one record


    ## Get tig k-mers ##

    return (
        sub_region,
        list(kanapy.util.kmer.stream(tig_seq, k_util, index=index))
    )

def fa_region(region, fa_file_name):
    """
    Get sequence from an indexed FASTA file. FASTA must have ".fai" index.

    :param region: Region object.
    :param fa_file_name: FASTA file name.

    :return: String sequence.
    """

    with pysam.FastaFile(fa_file_name) as fa_file:
        return fa_file.fetch(region.chrom, region.pos, region.end)

def subseq_region(region, aln_file_name, subseq_exe):
    """
    Given reference coordinates, find corresponding coordinates and sequence on aligned contigs.

    :param region: Reference region.
    :param aln_file_name: Alignment file name (CRAM/BAM).
    :param subseq_exe: Path to subseq executable.

    :return: An iterator for each contig aligned to the region. Each element of the iterator is a tuple of the
        contig region (Region object, contig context) and the sequence from the contig over that region.
    """

    # Open subseq
    with subprocess.Popen(
        """{subseq} -r {region} {tig_align}""".format(
            subseq=subseq_exe,
            region=str(region),
            tig_align=aln_file_name
        ).split(' '),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    ) as proc:

        outs, errs = proc.communicate()

        if proc.returncode != 0:
            print(errs.decode(), file=sys.stderr)
            raise RuntimeError(
                'Subseq failed extracting region from contig alignments {}: Return code = {}'.format(str(region), proc.returncode))

    # Parse subseq output (one or more FASTA records)
    with io.StringIO(outs.decode()) as in_file:
        for record in Bio.SeqIO.parse(in_file, 'fasta'):

            # Parse record ID as region
            region_match = re.match(r'^(.*):(\d+)-(\d+)$', record.id)

            if region_match is None:
                raise RuntimeError('Subseq return unrecognized region format in record ID: {}'.format(record.id))

            # Return record
            yield (
                Region(region_match[1], region_match[2], region_match[3]),
                str(record.seq)
            )

def cigar_lift_to_subject(subject_region, query_region, aln_file_name, ref_fa, query_start=None):
    """
    Get a list of coordinates through an alignment that match a query sequence to the subject (reference) it is aligned
    to.

    :param subject_region: Subject (reference) region the query is aligned to.
    :param query_region: Query region to characterize. All coordinates within this region are returned.
    :param aln_file_name: Alignment file name.
    :param ref_fa: Reference file name.
    :param query_start: Choose an alignment only if the first query (contig) base mapped is this value. Setting this
        value to the known start position of a query record allows overlapping alignments to be handled if the
        contig start is known ahead of time. If `None`, choose any alignment and generate an error if more than two
        contigs with the same name (`query_region.chrom`) are mapped to the same `subject_region` region.

    :return: A list of tuples describing the coordinates of `subject_region` and `query_region` through an alignment
        within those regions. Each tuple item includes (query start, query end, subject start, subject end, cigar
        op code, cigar op length).
    """

    lift_list = None

    # Open alignment file
    with pysam.AlignmentFile(aln_file_name, 'rb', reference_filename=ref_fa) as aln_file:
        for record in aln_file.fetch(subject_region.chrom, subject_region.pos, subject_region.end):
            if (
                    record.query_name == query_region.chrom and
                    record.query_alignment_start <= query_region.pos and
                    record.query_alignment_end >= query_region.end
            ):

                # Skip if query_start_pos does not match
                if query_start is not None and record.query_alignment_start != query_start:
                    continue

                # Check for overlapping alignment records
                if lift_list is not None:
                    raise RuntimeError('More than one matching query sequence "{}" aligned within "{}"'.format(
                        query_region, subject_region
                    ))

                # Build a list of records
                lift_list = lift_list_from_alignment_record(record, query_region)

    # Return list
    return lift_list

def get_query_sequence(subject_region, query_name, aln_file_name, ref_fa):
    """
    Get query sequence from an alignment record.

    :param subject_region: Subject (reference) region of the aligned contig. Contig must align to start and end
        positions.
    :param query_name: Contig name.
    :param aln_file_name: Alignment file name.
    :param ref_fa: Reference FASTA.

    :return: Sequence or None if no sequence was found.
    """

    with pysam.AlignmentFile(aln_file_name, 'rb', reference_filename=ref_fa) as aln_file:
        for record in aln_file.fetch(subject_region.chrom, subject_region.pos, subject_region.end):
            if (
                    record.query_name == query_name and
                    record.reference_start == subject_region.pos and
                    record.reference_end == subject_region.end
            ):
                return record.seq

    # No records found
    return None

def lift_list_from_alignment_record(record, query_region):
    """
    Generatea list of match reference/subject coordinates. See `cigar_lift_to_subject()`.

    :param record: Alignment record.
    :param query_region: Region to generate a lift-list for.

    :return: List of coordinates. See `cigar_lift_to_subject()`.
    """
    lift_list = list()

    query_pos = 0  # Code reads clipped bases and adds, do not start at record.query_alignment_start
    subject_pos = record.reference_start

    query_pos_start = query_region.pos
    query_pos_end = query_region.end

    cigar_index = 0

    for cigar_op, cigar_len in record.cigar:

        # Stop processing if this record start beyond the final position in query_region
        if query_pos > query_pos_end:
            break

        last_query_pos = query_pos
        last_subject_pos = subject_pos

        # Update positions on this CIGAR op
        if cigar_op in {0, 7, 8}:  # {M, = , X}
            # Consume both
            query_pos += cigar_len
            subject_pos += cigar_len

        elif cigar_op in {1, 4, 5}:  # {I, S, H}
            # Consume query only
            query_pos += cigar_len

        elif cigar_op in {2, 3}:  # {D, N}
            # Consume subject only
            subject_pos += cigar_len

        elif cigar_op != 6:  # P

            # Unrecognized opcode
            raise RuntimeError(
                'Unrecognized opcode in CIGAR string for alignment {}:{}-{} ({}:{}-{}) at CIGAR index {}: {} (oplen = {})'.format(
                    record.query_name, record.query_alignment_start, record.query_alignment_end,
                    record.reference_name, record.reference_start, record.reference_end,
                    cigar_index, cigar_op, cigar_len
                )
            )

        # Emmit record
        if cigar_op not in {4, 5}:
            if query_pos >= query_pos_start:
                lift_list.append((
                    last_query_pos, query_pos,
                    last_subject_pos, subject_pos,
                    cigar_op, cigar_len
                ))

        # Track CIGAR index for errors
        cigar_index += 1

    return lift_list

def get_lift_tree(lift_list, key_query=True):
    """
    Generate a tree of coordinates that can be queried.

    :param lift_list: List of lift tuples representing coordinates to lift:
        (query start, query end, subject start, subject end).
    :param query_ref: If `True`, generate a tree keyed by records by query (tig) coordinates, otherwise, key by
        subject (reference) coordinates.

    :return: Tree structure of lift coordinates.
    """

    lift_tree = intervaltree.IntervalTree()

    if key_query:
        for record in lift_list:
            if record[1] > record[0]:
                lift_tree[record[0]:record[1]] = record

    else:
        for record in lift_list:
            if record[3] > record[2]:
                lift_tree[record[2]:record[3]] = record

    return lift_tree

def tree_coords(query_pos, lift_tree, get_subject=True):
    """
    Get subject (reference) coordinates using a lift-tree (cigar_lift_to_subject converted to an
    intervaltree with get_lift_tree).

    :param query_pos: Positive on the query (contig).
    :param lift_tree: Lift tree.

    :return: Position on the subject (reference).
    """

    # Get intersecting record
    tree_record = lift_tree[query_pos]

    if len(tree_record) != 1:
        raise RuntimeError('Could not lift coordinate to reference "{}": Expected 1 matching record, found {}'.format(
            query_pos, len(tree_record)
        ))

    tree_record = list(tree_record)[0].data # From single-element set to the element (remove from set)

    # Iterpolate within record
    if get_subject:
        if tree_record[3] - tree_record[2] > 0:
            return tree_record[2] + (query_pos - tree_record[0])
        else:
            return tree_record[2]

    else:
        if tree_record[1] - tree_record[0] > 0:
            return tree_record[0] + (query_pos - tree_record[2])
        else:
            return tree_record[0]

def lift_coordinate(pos, subject_region, query_region, aln_file_name, ref_fa, to_subject=True):
    """
    Get subject (reference) coordinates using a lift-tree generated from a reference and contig region. If multiple
    coordinates must be generated from the region, then build the lift tree with `cigar_lift_to_subject` and
    `get_lift_tree', then query multiple times with `tree_coords`

    :param pos: Position to lift.
    :param subject_region: Subject (reference) region the query is aligned to.
    :param query_region: Query region to characterize. All coordinates within this region are returned.
    :param aln_file_name: Alignment file name.
    :param ref_fa: Reference file name.
    :param to_subject: `True` if lifting from query (aligned sequence) to subject (reference). If `False`, lift from
        subject to query.

    :return: Position on the subject (reference).
    """

    lift_list = cigar_lift_to_subject(
        subject_region,
        query_region,
        aln_file_name,
        ref_fa
    )

    lift_tree = get_lift_tree(lift_list, key_query=to_subject)

    return tree_coords(pos, lift_tree, get_subject=to_subject)


def get_matching_alignments(region_ref, region_tig, aln_file_name, ref_fa):
    """
    Get matching pysam alignment records for for two regions aligned.

    :param aln_file_name: Alignment file.
    :param ref_fa: Reference FASTA.
    :param region_ref: Reference region.
    :param region_tig: Contig region.

    :return: A list of pysam records for all contigs matching `region_tig` aligned to `region_ref`.
    """

    record_list = list()

    with pysam.AlignmentFile(aln_file_name, reference_filename=ref_fa) as aln_file:
        for record in aln_file.fetch(region_ref.chrom, region_ref.pos + 1, region_ref.end):

            if (
                    record.query_name == region_tig.chrom and
                    record.query_alignment_start <= region_tig.pos and
                    record.query_alignment_end >= region_tig.end
            ):
                record_list.append(record)

    return record_list
