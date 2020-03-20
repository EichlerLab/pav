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
        state_run_max=2000
    ):

    """
    Transform a k-mer stream from a section of a contig (tig_mer_stream) and a set of reference k-mers it aligns to
    (ref_kmer_set) and generate a smoothed density plot showing forward, reverse, and forward-reverse (reference has
    both the forward and reverse complement) k-mers. These smoothed states can then be used to call inversion
    boundaries.

    :param tig_mer_stream: A list of (k-mer, count) tuples from the contig region.
    :param ref_kmer_set: A set of k-mers in the reference region where the contig region aligns.
    :param k_util: K-mer utility from kanapy package.
    :param min_informative_kmers: Do not attempt density if the number of informative k-mers does not reach this limit.
        Informative k-mers are defined as k-mers that are in forward and/or reverse-complement in the reference
        k-mer set.
    :param density_smooth_factor: Smooth density by this factor. Density bandwidth is estimated by Scott's rule then
        multiplied by this factor to further smooth (> 1) or restrict smoothing (< 1). See
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
    :param state_run_max: Condense state consecutive runs of the same state longer than this length to one run of states
        with this length by removing the inner most states and leaving `state_run_max // 2` states on each side. If
        `None`, do not condense states.

    :return: A Pandas dataframe with columns "INDEX", "STATE", "
    """

    # Normalize state_run_max
    if state_run_max is not None:

        if state_run_max <= 1:
            raise RuntimeError(
                'Cannot condense runs of states: state_run_max must be 2 or greater: {}'.format(state_run_max)
            )

        state_run_max_flank = state_run_max // 2

        state_run_max = state_run_max_flank * 2  # If STATE_RUN_MAX was odd, ensure max run count is consistent with state_run_max_flank

    # Make dataframe
    df = pd.DataFrame(tig_mer_stream, columns=['KMER', 'INDEX'])

    # Assign state
    df['STATE'] = df['KMER'].apply(lambda kmer:
       KMER_ORIENTATION_STATE[
           int(kmer in ref_kmer_set),
           int(k_util.rev_complement(kmer) in ref_kmer_set)
       ]
   )

    # Clean k-mer column
    del (df['KMER'])

    # Subset to informative sites
    df = df.loc[df['STATE'] != -1]

    # Setup list to track removed records
    rm_record_list = list()

    # Remove low-count states. These cause spikes in density that rise above other states.
    state_count = df.groupby('STATE')['STATE'].count()

    for low_state in state_count[state_count < min_state_count].index:
        df_rm = df.loc[df['STATE'] == low_state].copy()

        df_rm['KERN_MAX'] = np.nan
        df_rm['KERN_FWD'] = np.nan
        df_rm['KERN_FWDREV'] = np.nan
        df_rm['KERN_REV'] = np.nan

        rm_record_list.append(df_rm)

        df = df.loc[df['STATE'] != low_state]

    # Check for number of informative k-mers before computing density
    if df.shape[0] < min_informative_kmers:
        return pd.DataFrame([], columns=['INDEX', 'STATE', 'KERN_MAX', 'KERN_FWD', 'KERN_FWDREV', 'KERN_REV'])

    # Subset runs of consecutive states
    if state_run_max is not None:

        # Get run-length encoded list of states
        state_rl = [element for element in rl_encoder(df, state_col='STATE') if element[1] > state_run_max]

        # Generate a table of states that should be removed
        df_long_run = pd.DataFrame(
            [
                (
                    element[2] + state_run_max_flank,
                    element[3] - state_run_max_flank,
                    element[0],
                    element[1] - state_run_max
                ) for element in state_rl
            ],
            columns=['INDEX_START', 'INDEX_END', 'STATE', 'RM_ELEMENTS']
        )

        for index, row in df_long_run.iterrows():

            # Save removed states
            df_rm = df.loc[(df['INDEX'] >= row['INDEX_START']) & (df['INDEX'] <= row['INDEX_END'])].copy()

            df_rm['KERN_MAX'] = df_rm['STATE']
            df_rm['KERN_FWD'] = np.nan
            df_rm['KERN_FWDREV'] = np.nan
            df_rm['KERN_REV'] = np.nan

            rm_record_list.append(df_rm)

            # Subset out removed states
            df = df.loc[(df['INDEX'] < row['INDEX_START']) | (df['INDEX'] > row['INDEX_END'])]

    # Setup bandwidth and index in informative k-mer space (ignore df['INDEX'] for density)
    density_bandwidth = df.shape[0] ** (-1.0 / 5.0) * density_smooth_factor

    df_pos_x  = np.arange(df.shape[0])  # Position for density. Index over informative k-mers only.

    # Get density
    if np.any(df['STATE'] == 0):
        kernel_fwd = scipy.stats.gaussian_kde(
            df_pos_x[df['STATE'] == 0],
            bw_method=density_bandwidth
        )
    else:
        kernel_fwd = lambda x: np.zeros(x.shape[0])

    if np.any(df['STATE'] == 1):
        kernel_fwdrev = scipy.stats.gaussian_kde(
            df_pos_x[df['STATE'] == 1],
            bw_method=density_bandwidth
        )
    else:
        kernel_fwdrev = lambda x: np.zeros(x.shape[0])

    if np.any(df['STATE'] == 2):
        kernel_rev = scipy.stats.gaussian_kde(
            df_pos_x[df['STATE'] == 2],
            bw_method=density_bandwidth
        )
    else:
        kernel_rev = lambda x: np.zeros(x.shape[0])

    # Scale kernel values by number of points (avoids penalizing states if they have fewer points in the whole region)
    kern_vals = np.array(
        [
            kernel_fwd(df_pos_x) * np.sum(df['STATE'] == 0),
            kernel_fwdrev(df_pos_x) * np.sum(df['STATE'] == 1),
            kernel_rev(df_pos_x) * np.sum(df['STATE'] == 2),
        ]
    )

    # Penalize spikes above 1.0.
    kern_vals[kern_vals > 1] = 1 / kern_vals[kern_vals > 1]

    # Assign to matrix
    df['KERN_MAX'] = np.argmax(kern_vals, 0)
    df['KERN_FWD'] = kern_vals[0, :]
    df['KERN_FWDREV'] = kern_vals[1, :]
    df['KERN_REV'] = kern_vals[2, :]

    # Check condensed states - Flanking max smoothed state must match the state of the rows that were removed
    if state_run_max is not None:
        for index, row in df_long_run.iterrows():

            state = row['STATE']
            state_up = df.loc[df.index < row['INDEX_START'], 'KERN_MAX'].iloc[0]
            state_dn = df.loc[df.index > row['INDEX_END'], 'KERN_MAX'].iloc[0]

            # Check smoothed state - must match original pre-smoothed state at each flank
            if state_up != state:
                raise RuntimeError((
                    'Condensed state reconstruction error: Index range {} - {} was condensed with state {}, but the '
                    'upstream flank in the condensed states does not match ({})'
                ).format(
                    row['INDEX_START'],
                    row['INDEX_END'],
                    state,
                    state_up
                ))

            if state_dn != state:
                raise RuntimeError((
                    'Condensed state reconstruction error: Index range {} - {} was condensed with state {}, but the '
                    'downstream flank in the condensed states does not match ({})'
                ).format(
                    row['INDEX_START'],
                    row['INDEX_END'],
                    state,
                    state_dn
                ))

    # Add dropped records
    if rm_record_list:
        df = pd.concat(
            [df] + rm_record_list
        ).sort_values('INDEX')

    # Return
    return df

def rl_encoder(df, state_col='KERN_MAX'):
    """
    Take a density table containing INDEX and a state column (STATE or KERN_MAX). Count consecutive STATEs and track indices for each run of a state.

    :param df: Dataframe of states.
    :param state_col: State column to get state from. Should be "STATE" for raw k-mer state or "KERN_MAX" for kernel-
        density max state at each locus. "STATE" is noisier, and "KERN_MAX" is smoothed.

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

# Overhauling method (2020-02-21) - rm commented code after
# def rl_encoder(val_list):
#     """
#     Take a list of values and return run-length encoded list of (value, length) where `value` is a value from the
#     list and `length` is the number of times it was repeated.
#
#     :param val_list: Values.
#
#     :return: A list of (value, length) tuples.
#     """
#
#     val = None
#     count = 0
#
#     for next_val in val_list:
#
#         if next_val == val:
#             count += 1
#         else:
#             if val is not None:
#                 yield (val, count)
#
#             val = next_val
#             count = 1
#
#     if val is not None:
#         yield (val, count)