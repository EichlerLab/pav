"""
Routines for calling inversions.
"""

import intervaltree
import numpy as np

import asmlib
import analib


#
# Constants
#

INITIAL_EXPAND = 4000      # Expand the flagged region by this much before starting.
EXPAND_FACTOR = 1.5        # Expand by this factor while searching

MIN_REGION_SIZE = 5000     # Extract a region at least this size
MAX_REGION_SIZE = 5000000  # Maximum region size

MAX_REF_KMER_COUNT = 100   # Skip low-complexity regions

MIN_INFORMATIVE_KMERS = 2000  # Minimum number of informative k-mers
MIN_KMER_STATE_COUNT = 20     # Remove states with fewer than this number of k-mers. Eliminates spikes in density

MIN_CONSECUTIVE_STATE_COUNT = 100    # When checking states, must have consecutive runs of the same state by this many k-mers. Used to remove noise for determining if region needs to be expanded or not.

DENSITY_SMOOTH_FACTOR = 1  # Estimate bandwith on Scott's rule then multiply by this factor to adjust bandwidth.


class InvCall:
    """
    Describes an inversion call with data supporting the call.

    :param region_ref_outer: Outer-flanks of inversion including inverted repeats, if present, that likely drove the
        inversion. The actual inversion breakpoints are likely between the outer and inner coordinates. This is the
        coordinate most inversion callers use.
    :param region_ref_inner: Inner coordinates of the inversion delineating the outermost boundary of strictly inverted
        sequence. This excludes reference and inverted repeats on the outer flanks of the inversion.
    :param region_tig_outer: Coordinates on the aligned contig of outer breakpoints corresponding to `region_ref_outer`.
    :param region_tig_inner: Coordinates on the aligned contig of inner breakpoints corresponding to `region_ref_inner`.

    :param region_ref_discovery: The reference region including the called inversion and surrounding unique
        sequence where the inversion was called.
    :param region_tig_discovery: The contig region matching `region_ref_discovery`.
    :param region_flag: The original flagged region, which was likely expanded to include the inversion and flanking
        unique sequence.

    :param df: A dataframe of k-mers and states for each k-mer with contig coordinate index relative to
        `region_tig_discovery`.
    """

    def __init__(
            self,
            region_ref_outer, region_ref_inner,
            region_tig_outer, region_tig_inner,
            region_ref_discovery, region_tig_discovery,
            region_flag,
            df
    ):

        self.region_ref_inner = region_ref_inner
        self.region_ref_outer = region_ref_outer

        self.region_tig_inner = region_tig_inner
        self.region_tig_outer = region_tig_outer

        self.region_ref_discovery = region_ref_discovery
        self.region_tig_discovery = region_tig_discovery
        self.region_flag = region_flag

        self.df = df

        self.id = '{}-{}-INV-{}'.format(region_ref_outer.chrom, region_ref_outer.pos + 1, len(region_ref_outer))

    def __repr__(self):
        return self.id


def scan_for_inv(region, ref_fa, aln_file_name, k_util, subseq_exe, log=None):
    """
    Scan region for inversions. Start with a flagged region (`region`) where variants indicated that an inversion
    might be. Scan that region for an inversion expanding as necessary.

    :param region: Flagged region to begin scanning for an inversion.
    :param ref_fa: Reference FASTA. Must also have a .fai file.
    :param aln_file_name: Contig alignment file.
    :param k_util: K-mer utility.
    :param subseq_exe: Subseq executable.
    :param log: Log file (open file handle).

    :return: A `InvCall` object describing the inversion found or `None` if no inversion was found.
    """

    # Init
    _write_log('Scanning for inversions in flagged region: {}'.format(region), log)

    region_flag = region.copy()  # Original flagged region

    df_fai = analib.ref.get_df_fai(ref_fa + '.fai')

    region.expand(INITIAL_EXPAND, min_pos=0, max_end=df_fai, shift=True)

    # Scan and expand
    while True:

        _write_log('Scanning region: {}'.format(region), log)

        ## Get reference k-mer counts ##
        ref_kmer_count = asmlib.seq.ref_kmers(region, ref_fa, k_util)

        if ref_kmer_count is None or len(ref_kmer_count) == 0:
            _write_log('No reference k-mers', log)
            return None

        # Skip low-complexity sites with repetitive k-mers
        max_mer_count = np.max(list(ref_kmer_count.values()))

        if max_mer_count > MAX_REF_KMER_COUNT:
            max_mer = [kmer for kmer, count in ref_kmer_count.items() if count == max_mer_count][0]

            _write_log('K-mer count exceeds max in {}: {} > {} ({})'.format(
                region_flag,
                max_mer_count,
                MAX_REF_KMER_COUNT,
                k_util.to_string(max_mer)
            ), log)

        ref_kmer_set = set(ref_kmer_count)


        ## Get contig k-mers ##

        tig_mer_region, tig_mer_stream = asmlib.seq.tig_mer_stream(region, aln_file_name, subseq_exe, k_util, index=True)

        if tig_mer_stream is None or len(tig_mer_stream) == 0:
            _write_log('No tig k-mers: May have expanded off contig boundaries', log)
            return None


        ## Density data frame ##

        df = asmlib.density.get_smoothed_density(
            tig_mer_stream,
            ref_kmer_set,
            k_util,
            min_informative_kmers=MIN_INFORMATIVE_KMERS,
            density_smooth_factor=DENSITY_SMOOTH_FACTOR,
            min_state_count=MIN_KMER_STATE_COUNT
        )

        # Note: States are 0 (fwd), 1 (fwd-rev), and 2 (rev) for k-mers found in forward orientation on the reference
        # region, in forward and reverse-complement, or reverese-complement, respectively.

        if df.shape[0] > 0:
            ## Check inversion ##

            # Get run-length encoded states (list of (state, count) tuples).
            state_rl = [record for record in asmlib.density.rl_encoder(df) if record[1] >= MIN_KMER_STATE_COUNT]
            condensed_states = [record[0] for record in state_rl]  # States only

            # Done if reference oriented k-mers (state == 0) found an both sides
            if len(condensed_states) > 2 and condensed_states[0] == 0 and condensed_states[-1] == 0:
                break

            # Expand
            last_len = len(region)
            expand_bp = np.int32(len(region) * EXPAND_FACTOR)

            if len(condensed_states) > 2:
                # More than one state. Expand disproprtionately if reference was found up or downstream.

                if condensed_states[0] == 0:
                    region.expand(
                        expand_bp, min_pos=0, max_end=df_fai, shift=True, balance=0.25
                    )  # Ref upstream: +25% upstream, +75% downstream

                elif condensed_states[-1] == 0:
                    region.expand(
                        expand_bp, min_pos=0, max_end=df_fai, shift=True, balance=0.75
                    )  # Ref downstream: +75% upstream, +25% downstream

                else:
                    region.expand(
                        expand_bp, min_pos=0, max_end=df_fai, shift=True, balance=0.5
                    )  # +50% upstream, +50% downstream

            else:
                region.expand(expand_bp, min_pos=0, max_end=df_fai, shift=True, balance=0.5)  # +50% upstream, +50% downstream

            if len(region) == last_len:
                # Stop if expansion had no effect

                _write_log(
                    'Reached reference limits, cannot expand',
                    log
                )

                return None

            # Continue with next iteration
        else:
            # No kmers, expand
            expand_bp = np.int32(len(region) * EXPAND_FACTOR)

            region.expand(
                expand_bp, min_pos=0, max_end=df_fai, shift=True, balance=0.5
            )  # +50% upstream, +50% downstream

    ## Characterize found region ##
    # Stop if no inverted sequence was found
    if not np.any([record[0] == 2 for record in state_rl]):
        _write_log('No inverted sequence found', log)
        return None

    # Call inversion
    if state_rl[0][0] != 0 or state_rl[-1][0] != 0:
        raise RuntimeError('Found INV region not flanked by reference sequence (program bug): {}'.format(region))

    # Find inverted repeat on left flank (upstream)
    coord_outer_l_pos = state_rl[1][2] + tig_mer_region.pos
    coord_inner_l_pos = coord_outer_l_pos

    index = 1

    while state_rl[index][0] != 2:
        coord_inner_l_pos = state_rl[index][3] + tig_mer_region.pos
        index += 1

    # Find inverted repeat on right flank
    coord_outer_r_pos = state_rl[-2][3] + tig_mer_region.pos
    coord_inner_r_pos = coord_outer_r_pos

    index = -2

    while state_rl[index][0] != 2:
        coord_inner_r_pos = state_rl[index][2] + tig_mer_region.pos
        index -= 1

    # Create contig coordinate records
    region_tig_outer = asmlib.seq.Region(tig_mer_region.chrom, coord_outer_l_pos, coord_outer_r_pos)
    region_tig_inner = asmlib.seq.Region(tig_mer_region.chrom, coord_inner_l_pos, coord_inner_r_pos)

    # Create an intervaltree for translating contig coordinates to reference coordinates
    lift_list = asmlib.seq.cigar_lift_to_subject(region, tig_mer_region, aln_file_name, ref_fa)

    lift_tree = intervaltree.IntervalTree()

    for record in lift_list:
        if record[1] > record[0]:
            lift_tree[record[0]:record[1]] = record

    # Get inner and outer reference coordinate breakpoints
    region_ref_outer = asmlib.seq.Region(
        region.chrom,
        tree_coords(region_tig_outer.pos, lift_tree),
        tree_coords(region_tig_outer.end, lift_tree)
    )

    region_ref_inner = asmlib.seq.Region(
        region.chrom,
        tree_coords(region_tig_inner.pos, lift_tree),
        tree_coords(region_tig_inner.end, lift_tree)
    )

    # Return inversion call
    return InvCall(
        region_ref_outer, region_ref_inner,
        region_tig_outer, region_tig_inner,
        region, tig_mer_region,
        region_flag, df
    )


def tree_coords(query_pos, lift_tree):
    """
    Get subject (reference) coordinates using a lift-tree (asmlib.seq.cigar_lift_to_subject converted to an
    intervaltree).

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
    if tree_record[3] - tree_record[2] > 0:
        return tree_record[2] + (query_pos - tree_record[0])
    else:
        return tree_record[2]




def _write_log(message, log):
    """
    Write message to log.

    :param message: Message.
    :param log: Log or `None`.
    """

    if log is None:
        return

    log.write(message)
    log.write('\n')

    log.flush()
