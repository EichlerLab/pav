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
        min_state_count=20
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

    :return: A Pandas dataframe with columns "INDEX", "STATE", "
    """

    # Make dataframe
    df = pd.DataFrame(tig_mer_stream, columns=['KMER', 'INDEX'])

    # Assign state
    df['STATE'] = df['KMER'].apply(lambda kmer:
       KMER_ORIENTATION_STATE[
           int(kmer in ref_kmer_set),
           int(k_util.rev_complement(kmer) in ref_kmer_set)
       ]
   )

    # Subset to informative sites
    df = df.loc[df['STATE'] != -1]

    # Remove low-count states. These cause spikes in density that rise above other states.
    state_count = df.groupby('STATE')['STATE'].count()

    for low_state in state_count[state_count < min_state_count].index:
        df = df.loc[df['STATE'] != low_state]

    # Clean k-mer column
    del (df['KMER'])

    # Check for number of informative k-mers before computing density
    if df.shape[0] < min_informative_kmers:
        return pd.DataFrame([], columns=['INDEX', 'STATE', 'KERN_MAX', 'KERN_FWD', 'KERN_FWDREV', 'KERN_REV'])

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

    # Return
    return df

def rl_encoder(df):
    """
    Take a density table containing INDEX and STATE. Count consecutive STATEs and track indices for each run of a state.

    :param df: Dataframe of states.

    :return:
    """

    state = None
    count = 0
    pos = None
    end = None

    for index, row in df[['STATE', 'INDEX']].iterrows():  # STATE/INDEX columns make row Series int (removing float columns)

        if row['STATE'] == state:
            count += 1
            end = row['INDEX']

        else:
            if state is not None:
                yield (state, count, pos, end)

            state = row['STATE']
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