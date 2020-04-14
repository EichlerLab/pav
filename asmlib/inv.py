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
MAX_REGION_SIZE = 1000000  # Maximum region size

MAX_REF_KMER_COUNT = 100   # Skip low-complexity regions

MIN_INFORMATIVE_KMERS = 2000  # Minimum number of informative k-mers
MIN_KMER_STATE_COUNT = 20     # Remove states with fewer than this number of k-mers. Eliminates spikes in density

DENSITY_SMOOTH_FACTOR = 1  # Estimate bandwith on Scott's rule then multiply by this factor to adjust bandwidth.

MIN_INV_KMER_RUN = 100  # States must have a continuous run of this many strictly inverted k-mers

MIN_TIG_REF_PROP = 0.6  # The contig and reference region sizes must be within this factor (reciprocal) or the event
                        # is likely unbalanced (INS or DEL) and would already be in the callset

MIN_EXP_COUNT = 3  # The number of region expansions to try (including the initial expansion) and finding only fwd k-mer
                   # states after smoothing before giving up on the region.

# Matrix converting k-mer location to UP/DN match. Assign
# NA to missing or k-mers in both.
#
# KMER_LOC_STATE[in-upstream, in-dnstream]
KMER_LOC_STATE = np.asarray(
    [
        ['NA', 'OTHER'],
        ['SAME', 'NA']
    ]
)


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

        # Save INV regions and dataframe
        self.region_ref_outer = region_ref_outer
        self.region_ref_inner = region_ref_inner

        self.region_tig_outer = region_tig_outer
        self.region_tig_inner = region_tig_inner

        self.region_ref_discovery = region_ref_discovery
        self.region_tig_discovery = region_tig_discovery

        self.region_flag = region_flag

        self.df = df

        # Generate length and ID
        self.svlen = len(region_ref_outer)
        self.id = '{}-{}-INV-{}'.format(region_ref_outer.chrom, region_ref_outer.pos + 1, self.svlen)

        # Get max INV density height
        self.max_inv_den_diff = np.max(
            df.loc[
                df['STATE'] == 2, 'KERN_REV'
            ] - df.loc[
                df['STATE'] == 2, ['KERN_FWD', 'KERN_FWDREV']
            ].apply(
                np.max, axis=1
            )
        )

    def __repr__(self):
        return self.id


def scan_for_inv(region, ref_fa, aln_file_name, k_util, subseq_exe, log=None, flag_id=None):
    """
    Scan region for inversions. Start with a flagged region (`region`) where variants indicated that an inversion
    might be. Scan that region for an inversion expanding as necessary.

    :param region: Flagged region to begin scanning for an inversion.
    :param ref_fa: Reference FASTA. Must also have a .fai file.
    :param aln_file_name: Contig alignment file.
    :param k_util: K-mer utility.
    :param subseq_exe: Subseq executable.
    :param log: Log file (open file handle).
    :param flag_id: ID of the flagged region (printed so flagged regions can easily be tracked in logs).

    :return: A `InvCall` object describing the inversion found or `None` if no inversion was found.
    """

    # Init
    if flag_id is None:
        flag_id = '<No ID>'

    _write_log('Scanning for inversions in flagged region: {} (flagged region record id = {})'.format(region, flag_id), log)

    region_flag = region.copy()  # Original flagged region

    df_fai = analib.ref.get_df_fai(ref_fa + '.fai')

    region.expand(INITIAL_EXPAND, min_pos=0, max_end=df_fai, shift=True)

    expansion_count = 0

    # Scan and expand
    while True:

        expansion_count += 1

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

            return None

        ref_kmer_set = set(ref_kmer_count)


        ## Get contig k-mers ##
        stream_tuple = asmlib.seq.tig_mer_stream(region, aln_file_name, subseq_exe, k_util, index=True)

        if stream_tuple is not None:
            tig_mer_region, tig_mer_stream = stream_tuple
        else:
            _write_log('No tig k-mers: No contig sequence record covers region', log)
            return None

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
            state_rl = [record for record in asmlib.density.rl_encoder(df)]
            condensed_states = [record[0] for record in state_rl]  # States only

            if len(state_rl) == 1 and state_rl[0][0] == 0 and expansion_count >= MIN_EXP_COUNT:
                _write_log(
                    'Found no inverted k-mer states after {} expansion(s)'.format(expansion_count),
                    log
                )

                return None

            # Done if reference oriented k-mers (state == 0) found an both sides
            if len(condensed_states) > 2 and condensed_states[0] == 0 and condensed_states[-1] == 0:
                break

            # Expand
            last_len = len(region)
            expand_bp = np.int32(len(region) * EXPAND_FACTOR)

            if len(condensed_states) > 2:
                # More than one state. Expand disproportionately if reference was found up or downstream.

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
            _write_log(
                'No informative reference k-mers in forward or reverse orientation in region',
                log
            )

            return None

    ## Characterize found region ##
    # Stop if no inverted sequence was found
    if not np.any([record[0] == 2 for record in state_rl]):
        _write_log('No inverted sequence found', log)
        return None

    max_inv_run = np.max([record[1] for record in state_rl if record[0] == 2])

    if max_inv_run < MIN_INV_KMER_RUN:
        _write_log('Longest run of strictly inverted k-mers ({}) does not meet the minimum threshold ({})'.format(
            max_inv_run, MIN_INV_KMER_RUN
        ), log)

        return None


    # Code check - must be flanked by reference sequence
    if state_rl[0][0] != 0 or state_rl[-1][0] != 0:
        raise RuntimeError('Found INV region not flanked by reference sequence (program bug): {}'.format(region))

    # Subset to strictly inverted states (for calling inner breakpoints)
    state_rl_inv = [record for record in state_rl if record[0] == 2]

    # Find inverted repeat on left flank (upstream)
    tig_outer_up_stream = state_rl[1][2] + tig_mer_region.pos
    tig_inner_up_stream = state_rl_inv[0][2] + tig_mer_region.pos

    # Find inverted repeat on right flank
    tig_outer_dn_stream = state_rl[-2][3] + tig_mer_region.pos + k_util.k_size
    tig_inner_dn_stream = state_rl_inv[-1][3] + tig_mer_region.pos + k_util.k_size

    # Create contig coordinate records
    region_tig_outer = asmlib.seq.Region(tig_mer_region.chrom, tig_outer_up_stream, tig_outer_dn_stream)
    region_tig_inner = asmlib.seq.Region(tig_mer_region.chrom, tig_inner_up_stream, tig_inner_dn_stream)

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

    # Check size proportions
    if len(region_ref_outer) < len(region_tig_outer) * MIN_TIG_REF_PROP:
        _write_log(
            'Reference region too short: Reference region length ({:,d}) is not within {:.2f}% of the contig region length ({:,d})'.format(
                len(region_ref_outer),
                MIN_TIG_REF_PROP * 100,
                len(region_tig_outer)
            ),
            log
        )

        return None

    if len(region_tig_outer) < len(region_ref_outer) * MIN_TIG_REF_PROP:
        _write_log(
            'Contig region too short: Contig region length ({:,d}) is not within {:.2f}% of the reference region length ({:,d})'.format(
                len(region_tig_outer),
                MIN_TIG_REF_PROP * 100,
                len(region_ref_outer)
            ),
            log
        )

        return None

    # Get INV-DUP flanking annotation. Where there is an inverted duplication on the flanks, flag k-mers that belong
    # strictly to the upstream or downstream flanking duplication.
    df = annotate_inv_dup_mers(df, region_ref_outer, region_ref_inner, region_tig_outer, region_tig_inner, tig_mer_region, ref_fa, k_util)

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

def annotate_inv_dup_mers(
        df,
        region_ref_outer, region_ref_inner,
        region_tig_outer, region_tig_inner,
        region_tig_discovery,
        ref_fa, k_util
):

    # Get regions for duplications - ref
    region_dup_ref_up = asmlib.seq.Region(
        region_ref_outer.chrom,
        region_ref_outer.pos,
        region_ref_inner.pos
    )

    region_dup_ref_dn = asmlib.seq.Region(
        region_ref_outer.chrom,
        region_ref_inner.end,
        region_ref_outer.end
    )

    # Get regions for duplications - tig
    region_dup_tig_up = asmlib.seq.Region(
        region_tig_outer.chrom,
        region_tig_outer.pos,
        region_tig_inner.pos
    )

    region_dup_tig_dn = asmlib.seq.Region(
        region_tig_outer.chrom,
        region_tig_inner.end,
        region_tig_outer.end
    )

    # Get reference k-mer sets (canonical k-mers)
    ref_set_up = {
        k_util.canonical_complement(kmer) for kmer in asmlib.seq.ref_kmers(region_dup_ref_up, ref_fa, k_util).keys()
    }

    ref_set_dn = {
        k_util.canonical_complement(kmer) for kmer in asmlib.seq.ref_kmers(region_dup_ref_dn, ref_fa, k_util).keys()
    }

    # Set canonical k-mers in df
    df['KMER_CAN'] = df['KMER'].apply(lambda kmer: k_util.canonical_complement(kmer))

    # Add contig index
    df['TIG_INDEX'] = df['INDEX'] + region_tig_discovery.pos

    # Annotate upstream and downstream inverted duplications
    df['FLANK'] = np.nan

    df.loc[
        (df['TIG_INDEX'] >= region_dup_tig_up.pos) & (df['TIG_INDEX'] < region_dup_tig_up.end - k_util.k_size),
        'FLANK'
    ] = 'UP'

    df.loc[
        (df['TIG_INDEX'] >= region_dup_tig_dn.pos) & (df['TIG_INDEX'] < region_dup_tig_dn.end - k_util.k_size),
        'FLANK'
    ] = 'DN'

    # Annotate upstream/downstream k-mer matches
    df['MATCH'] = np.nan

    df.loc[df['FLANK'] == 'UP', 'MATCH'] = df.loc[
        df['FLANK'] == 'UP', 'KMER'
    ].apply(lambda kmer:
        KMER_LOC_STATE[
            int(kmer in ref_set_up),
            int(kmer in ref_set_dn)
        ]
    )

    df.loc[df['FLANK'] == 'DN', 'MATCH'] = df.loc[
        df['FLANK'] == 'DN', 'KMER'
    ].apply(lambda kmer:
        KMER_LOC_STATE[
            int(kmer in ref_set_dn),
            int(kmer in ref_set_up)
        ]
    )

    df['MATCH'] = df['MATCH'].apply(lambda val: val if val != 'NA' else np.nan)

    # Return updated DataFrame
    del(df['KMER_CAN'])
    del(df['TIG_INDEX'])

    return df


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
