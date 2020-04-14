"""
K-mer density used for calling inversions.
"""

import numpy as np
import pandas as pd

import scipy.stats

# K-mer orientation matrix: Tig k-mer against forward (vertical axis)
# and reverse (horizontal axis) k-mers in the reference
#
#       T         F
#   ---------------------------------
# T | No data (-1)  | Rev (2)       |
# F | Reference (0) | Fwd & Rev (1) |
#   ---------------------------------
KMER_ORIENTATION_STATE = np.asarray(
    [
        [-1, 2],
        [0, 1]
    ]
)


def get_smoothed_density(
        tig_mer_stream,
        ref_kmer_set,
        k_util,
        min_informative_kmers=2000,
        density_smooth_factor=1,
        min_state_count=20,
        state_run_smooth=20,
        state_run_smooth_delta=0.005
    ):

    """
    Transform a k-mer stream from a section of a contig (tig_mer_stream) and a set of reference k-mers it aligns to
    (ref_kmer_set) and generate a smoothed density plot showing forward, reverse, and forward-reverse (reference has
    both the forward and reverse complement) k-mers. These smoothed states can then be used to call inversion
    boundaries.

    There are two "spaces" indexed:
    * k-mer: The set of all k-mers indexed from the first k-mer in the sequence to the last. Skipped k-mers also skips
        indices in this space. Adding the genomic position of the first k-mer gives coordinates in the contig for each
        k-mer. This is the "INDEX" column in the output DataFrame.
    * density: All missing and uninformative k-mers (k-mers not in the reference region in either orientation) are
        removed and records are re-indexed starting from 0. This keeps the density from dropping off if there is
        disagreement between the sequence and reference (e.g. SVs or indels in the sequence). This index is "INDEX_DEN"
        in the output DataFrame. It is useful for following the density smoothing process, but is not used outside this
        function.
    
    Dataframe fields:
    * INDEX: K-mer index from the k-mer stream (skipped k-mers also skips indices).
    * INDEX_DEN: Index in density space. All uninformative k-mers are removed and remaining rows are re-indexed starting
        from 0. The kernel density function is computed over this range.
    * STATE_MER: State of the k-mer record.
    * STATE: Kernel density smoothed states.
    * KERN_FWD: Kernel density of forward-oriented k-mers.
    * KERN_FWDREV: Kernel density of forward- and reverse-oriented k-mers (found in both orientations in the reference).
    * KERN_REV: Kernel density of reverse-oriented k-mers.
    * INTERP: `True` if linear interpolation was used to calculated kernel values ("KERN_" columns). This is for
        troubleshooting inversion calls where slight changes in density might affect the call (not used by the
        pipeline). Interpolation for records with a `True` value were interpolated using values from the two nearest
        non-interpolated values (nearest upstream and downstream).

    :param tig_mer_stream: A list of (k-mer, count) tuples from the contig region.
    :param ref_kmer_set: A set of k-mers in the reference region where the contig region aligns.
    :param k_util: K-mer utility from kanapy package.
    :param min_informative_kmers: Do not attempt density if the number of informative k-mers does not reach this limit.
        Informative k-mers are defined as k-mers that are in forward and/or reverse-complement in the reference
        k-mer set.
    :param density_smooth_factor: Smooth density by this factor. Density bandwidth is estimated by Scott's rule then
        multiplied by this factor to further smooth (> 1) or restrict smoothing (< 1). See
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
    :param state_run_smooth: Smooth density calculations by this many records. Initially, calculate density only once
        in a run of this many k-mers. If there is a state change or a delta larger than `state_run_smooth_delta`, then
        density is filled in for all skipped elements. If there is no state change and the delta is small, use linear
        interpolation to fill in missing values.
    :param state_run_smooth_delta: Changes between state densities by this much or more will be filled in with actual
        density values instead of interpolated. See `state_run_smooth`.

    :return: A Pandas dataframe describing the density. The dataframe index is set to the "INDEX" column (index in
        k-mer/tig space).
    """

    # Make dataframe
    df = pd.DataFrame(tig_mer_stream, columns=['KMER', 'INDEX'])

    df['STATE'] = -1

    # Assign state
    df['STATE_MER'] = df['KMER'].apply(lambda kmer:
       KMER_ORIENTATION_STATE[
           int(kmer in ref_kmer_set),
           int(k_util.rev_complement(kmer) in ref_kmer_set)
       ]
    )

    # Subset to informative sites
    df = df.loc[df['STATE_MER'] != -1]

    # Remove low-count states. These cause spikes in density that rise above other states.
    state_count = df.groupby('STATE_MER')['STATE_MER'].count()

    for low_state in state_count[state_count < min_state_count].index:
        df_rm = df.loc[df['STATE_MER'] == low_state].copy()

        df_rm['KERN_FWD'] = np.nan
        df_rm['KERN_FWDREV'] = np.nan
        df_rm['KERN_REV'] = np.nan

        df = df.loc[df['STATE_MER'] != low_state]

    # Check for number of informative k-mers before computing density
    if df.shape[0] < min_informative_kmers:
        return pd.DataFrame([], columns=['INDEX', 'STATE_MER', 'STATE', 'KERN_FWD', 'KERN_FWDREV', 'KERN_REV', 'KMER'])

    # Setup bandwidth and index in informative k-mer space (ignore df['INDEX'] for density)
    density_bandwidth = df.shape[0] ** (-1.0 / 5.0) * density_smooth_factor

    # Index in condensed space
    df.insert(1, 'INDEX_DEN', np.arange(df.shape[0]))

    df.set_index('INDEX_DEN', inplace=True, drop=False)  # Switch index to density space

    # Get density functions
    if np.any(df['STATE_MER'] == 0):
        kernel_fwd = scipy.stats.gaussian_kde(
            df.loc[df['STATE_MER'] == 0, 'INDEX_DEN'],
            bw_method=density_bandwidth
        )
    else:
        kernel_fwd = lambda x: np.zeros(np.array(x).shape[0])

    if np.any(df['STATE_MER'] == 1):
        kernel_fwdrev = scipy.stats.gaussian_kde(
            df.loc[df['STATE_MER'] == 1, 'INDEX_DEN'],
            bw_method=density_bandwidth
        )
    else:
        kernel_fwdrev = lambda x: np.zeros(np.array(x).shape[0])

    if np.any(df['STATE_MER'] == 2):
        kernel_rev = scipy.stats.gaussian_kde(
            df.loc[df['STATE_MER'] == 2, 'INDEX_DEN'],
            bw_method=density_bandwidth
        )
    else:
        kernel_rev = lambda x: np.zeros(np.array(x).shape[0])

    # Assign density on sampled sites
    sum_state_fwd = np.sum(df['STATE_MER'] == 0)
    sum_state_fwdrev = np.sum(df['STATE_MER'] == 1)
    sum_state_rev = np.sum(df['STATE_MER'] == 2)

    sample_index_list = list(df.loc[df['INDEX_DEN'] % state_run_smooth == 0, 'INDEX_DEN'])

    if sample_index_list[-1] != np.int32(df.iloc[-1]['INDEX_DEN']):  # Add last table element if it does not exist
        sample_index_list.append(np.int32(df.iloc[-1]['INDEX_DEN']))

    df['STATE'] = -1
    df['KERN_FWD'] = np.nan
    df['KERN_FWDREV'] = np.nan
    df['KERN_REV'] = np.nan
    df['INTERP'] = False

    df.loc[sample_index_list, 'KERN_FWD'] = kernel_fwd(sample_index_list) * sum_state_fwd
    df.loc[sample_index_list, 'KERN_FWDREV'] = kernel_fwdrev(sample_index_list) * sum_state_fwdrev
    df.loc[sample_index_list, 'KERN_REV'] = kernel_rev(sample_index_list) * sum_state_rev

    # Get max state
    df.loc[sample_index_list, 'STATE'] = df.loc[
        sample_index_list, ['KERN_FWD', 'KERN_FWDREV', 'KERN_REV']
    ].apply(
        lambda vals: np.argmax(vals.array), axis=1
    )

    # Fill missing values
    for index in range(len(sample_index_list) - 1):

        range_start = sample_index_list[index]
        range_end = sample_index_list[index + 1]

        outer_range = np.arange(range_start, range_end + 1)  # All values including the first and last with densities
        inner_range = np.arange(range_start + 1, range_end)  # Inner values (without densities calculated)

        if range_end == range_start + 1:
            continue  # Indexes are next to each other, all records in range have density calculations

        state_change = \
            len(set(df.loc[outer_range, 'STATE_MER'])) > 1 or \
            df.loc[range_start, 'STATE'] != df.loc[range_end, 'STATE']

        density_change = np.max(np.abs(
            df.loc[range_start, ['KERN_FWD', 'KERN_FWDREV', 'KERN_REV']] - df.loc[range_end, ['KERN_FWD', 'KERN_FWDREV', 'KERN_REV']]
        )) > state_run_smooth_delta

        if state_change or density_change:
            # Calculate densities within range

            df.loc[inner_range, 'KERN_FWD'] = kernel_fwd(inner_range) * sum_state_fwd
            df.loc[inner_range, 'KERN_FWDREV'] = kernel_fwdrev(inner_range) * sum_state_fwdrev
            df.loc[inner_range, 'KERN_REV'] = kernel_rev(inner_range) * sum_state_rev

        else:
            # Interpolate densities within range

            df.loc[inner_range, 'KERN_FWD'] = np.interp(
                inner_range,
                [range_start, range_end],
                [df.loc[range_start, 'KERN_FWD'], df.loc[range_end, 'KERN_FWD']]
            )

            df.loc[inner_range, 'KERN_FWDREV'] = np.interp(
                inner_range,
                [range_start, range_end],
                [df.loc[range_start, 'KERN_FWDREV'], df.loc[range_end, 'KERN_FWDREV']]
            )

            df.loc[inner_range, 'KERN_REV'] = np.interp(
                inner_range,
                [range_start, range_end],
                [df.loc[range_start, 'KERN_REV'], df.loc[range_end, 'KERN_REV']]
            )

            df.loc[inner_range, 'INTERP'] = True

    # Penalize spikes above 1.0.
    df.loc[df['KERN_FWD'] > 1.0, 'KERN_FWD'] = 1 / df.loc[df['KERN_FWD'] > 1.0]
    df.loc[df['KERN_FWDREV'] > 1.0, 'KERN_FWDREV'] = 1 / df.loc[df['KERN_FWDREV'] > 1.0]
    df.loc[df['KERN_REV'] > 1.0, 'KERN_REV'] = 1 / df.loc[df['KERN_REV'] > 1.0]

    # Set max densities on all
    df['STATE'] = df[['KERN_FWD', 'KERN_FWDREV', 'KERN_REV']].apply(
        lambda vals: np.argmax(vals.array),
        axis=1
    )

    # Column order
    df = df[['INDEX', 'STATE_MER', 'STATE', 'KERN_FWD', 'KERN_FWDREV', 'KERN_REV', 'KMER']]

    # Return
    df.set_index(df['INDEX'], inplace=True, drop=False)

    return df

def rl_encoder(df, state_col='STATE'):
    """
    Take a density table containing INDEX and a state column (STATE_MER or STATE). Count consecutive states and track indices for each run of a state.

    :param df: Dataframe of states.
    :param state_col: State column to get state from. Should be "STATE_MER" for raw k-mer state or "STATE" for kernel-
        density max state at each locus. "STATE_KMER" is noisier, and "STATE" is smoothed.

    :return: Iterator of (state, count, pos, end) tuples.
    """

    state = None
    count = 0
    pos = None
    end = None

    for index, row in df[[state_col, 'INDEX']].iterrows():  # STATE/INDEX columns make row Series int (removing float columns so each row is integers)

        if row[state_col] == state:
            count += 1
            end = row['INDEX']

        else:
            if state is not None:
                yield (state, count, pos, end)

            state = row[state_col]
            count = 1
            pos = end = row['INDEX']

    if state is not None:
        yield (state, count, pos, end)
