"""
K-mer density estimation.
"""

import importlib.util
import inspect
import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy.stats

# K-mer orientation matrix: Tig k-mer against forward (vertical axis)
# and reverse (horizontal axis) k-mers in the reference
#
#       F         T
#   ---------------------------------
# F | No data (-1)  | Rev (2)       |
# T | Reference (0) | Fwd & Rev (1) |
#   ---------------------------------
KMER_ORIENTATION_STATE = np.asarray(
    [
        [-1, 2],
        [0, 1]
    ]
)

SAMPLE_INDEX_CHUNK_SIZE = 400

class KdeTruncNorm(object):
    """
    Kernel Density Estimation (KDE) with a truncated-normal distribution. Uses FFT convolution to solve density.
    """

    def __init__(self, bandwidth=100.0, trunc_z=3.0, conv='auto'):

        # Check parameters
        if bandwidth <= 0:
            raise ValueError(f'Bandwidth must be > 0: {bandwidth}')

        if trunc_z <= 0:
            raise ValueError(f'Truncation SD must be > 0: {trunc_z}')

        self.bandwidth = float(bandwidth)
        self.band_sd = float(trunc_z)

        # Get convolution method
        if isinstance(conv, str):
            conv_lower = conv.strip().lower()

            is_auto = conv_lower == 'auto'

            self.conv_method = None

            if is_auto:
                conv_lower = 'fft'

            if conv_lower == 'fft' and self.conv_method is None:

                spec = importlib.util.find_spec('scipy.signal')

                if spec is None:
                    if not is_auto:
                        raise RuntimeError(f'Error initializing KdeTruncNorm: Missing package for KDE convolution method {conv}: scipy.signal')
                    else:
                        conv_upper = 'conv'  # Try next

                import scipy.signal
                self.conv_method = scipy.signal.fftconvolve

            if conv_lower == 'conv' and self.conv_method is None:
                self.conv_method = np.convolve

            if self.conv_method is None:
                if is_auto:
                    raise RuntimeError(f'Error initializing KdeTruncNorm: Could not automatically resolve convolution method')

                raise RuntimeError(f'Error initializing KdeTruncNorm: Unknown convolution method: {conv}')

        elif hasattr(conv, '__call__'):
            self.conv_method = conv

        # Inspect convolution method
        n_arg = len(inspect.getfullargspec(self.conv_method).args)
        n_default = len(inspect.getfullargspec(self.conv_method).defaults)

        if n_arg < 2:
            raise RuntimeError(f'Convolution method does not take at least 2 arguments (n = {n_arg}): {self.conv_method}')

        if n_arg - n_default > 2:
            raise RuntimeError(f'Convolution method requires more than 2 arguments (n = {n_arg - n_default} without default values): {self.conv_method}')

        # Set normal vector
        tnorm = scipy.stats.truncnorm(-trunc_z, trunc_z)  # Truncate at Z = band_sd

        self.band_bound = int(np.ceil(trunc_z * bandwidth))  # number of positions in the truncated normal vector

        self.v_tnorm = tnorm.pdf(
            np.arange(-self.band_bound, self.band_bound + 1) / self.bandwidth  # Range from -band_sd to +band_sd after accounting for bandwith
        )  # Pre-compute PDF at positions

        self.v_tnorm = self.v_tnorm / self.v_tnorm.sum()  # Vector sums to 1 (density 1 is a position with ones from -band_sd to band_sd)

    def __call__(self, x, n=None):
        """
        Get the density of an n-length vector with ones at x positions and zeros elsewhere.

        :param x: Location of 1's in the array (if `n` is defined) or an array to directly conolve with the density
            kernel and `n` becomes the length of this array.
        :param n: Length of the array. If `None`, `x` is the full-length array to convolve.

        :return: A density vector of length `n`.
        """

        if n is not None:
            v_state = np.zeros(n)

            for v in x:
                v_state[v] = 1
        else:
            if not isinstance(x, np.ndarray):
                v_state = np.array(x)
            else:
                v_state = x

        y = self.conv_method(v_state, self.v_tnorm)

        return y[self.band_bound:-self.band_bound]

def get_smoothed_density(
        tig_mer_stream,
        ref_kmer_set,
        k_util,
        threads=1,
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

    :param tig_mer_stream: A list of (k-mer, count) tuples from the contig region.
    :param ref_kmer_set: A set of k-mers in the reference region where the contig region aligns.
    :param k_util: K-mer utility from kanapy package.
    :param threads: Number of threads to use for computing densities.
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

    # Objects used by threads
    kernel_dict = None

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
    if df.shape[0] < min_informative_kmers or np.all(df['STATE_MER'] == 0):
        return pd.DataFrame([], columns=['INDEX', 'STATE_MER', 'STATE', 'KERN_FWD', 'KERN_FWDREV', 'KERN_REV', 'KMER'])

    # Setup bandwidth and index in informative k-mer space (ignore df['INDEX'] for density)
    density_bandwidth = df.shape[0] ** (-1.0 / 5.0) * density_smooth_factor

    # Index in condensed space
    df.insert(1, 'INDEX_DEN', np.arange(df.shape[0]))

    df.set_index('INDEX_DEN', inplace=True, drop=False)  # Switch index to density space

    # Count states for scaled-density
    sum_state_fwd = np.sum(df['STATE_MER'] == 0)
    sum_state_fwdrev = np.sum(df['STATE_MER'] == 1)
    sum_state_rev = np.sum(df['STATE_MER'] == 2)

    # Setup density models.
    # Initializes kernel_dict per-thread
    def init_density():
        global kernel_dict

        if kernel_dict is None:
            kernel_dict = dict()

            if np.any(df['STATE_MER'] == 0):
                kernel_dict['fwd'] = scipy.stats.gaussian_kde(
                    df.loc[df['STATE_MER'] == 0, 'INDEX_DEN'],
                    bw_method=density_bandwidth
                )
            else:
                kernel_dict['fwd'] = lambda x: np.zeros(np.array(x).shape[0])

            if np.any(df['STATE_MER'] == 1):
                kernel_dict['fwdrev'] = scipy.stats.gaussian_kde(
                    df.loc[df['STATE_MER'] == 1, 'INDEX_DEN'],
                    bw_method=density_bandwidth
                )
            else:
                kernel_dict['fwdrev'] = lambda x: np.zeros(np.array(x).shape[0])

            if np.any(df['STATE_MER'] == 2):
                kernel_dict['rev'] = scipy.stats.gaussian_kde(
                    df.loc[df['STATE_MER'] == 2, 'INDEX_DEN'],
                    bw_method=density_bandwidth
                )
            else:
                kernel_dict['rev'] = lambda x: np.zeros(np.array(x).shape[0])

    # Density calculation functions.
    def density_fwd(val):
        return kernel_dict['fwd'](val) * sum_state_fwd

    def density_fwdrev(val):
        return kernel_dict['fwdrev'](val) * sum_state_fwdrev

    def density_rev(val):
        return kernel_dict['rev'](val) * sum_state_rev

    # Interpolation functions.
    # Take a range of values where the flanks (values just before and after inner_range) ar filled in, then
    # interpolates density values between them.
    def interp_fwd(inner_range):
        range_start = inner_range[0] - 1
        range_end = inner_range[-1] + 1

        return \
        np.interp(
            inner_range,
            [range_start, range_end],
            [df.loc[range_start, 'KERN_FWD'], df.loc[range_end, 'KERN_FWD']]
        )

    def interp_fwdrev(inner_range):
        range_start = inner_range[0] - 1
        range_end = inner_range[-1] + 1

        return \
        np.interp(
            inner_range,
            [range_start, range_end],
            [df.loc[range_start, 'KERN_FWDREV'], df.loc[range_end, 'KERN_FWDREV']]
        )

    def interp_rev(inner_range):
        range_start = inner_range[0] - 1
        range_end = inner_range[-1] + 1

        return \
        np.interp(
            inner_range,
            [range_start, range_end],
            [df.loc[range_start, 'KERN_REV'], df.loc[range_end, 'KERN_REV']]
        )

    # Setup indexes (initial non-sampled sites density is initially calcuated over)
    sample_index_list = list(df.loc[df['INDEX_DEN'] % state_run_smooth == 0, 'INDEX_DEN'])

    if sample_index_list[-1] != np.int32(df.iloc[-1]['INDEX_DEN']):  # Add last table element if it does not exist
        sample_index_list.append(np.int32(df.iloc[-1]['INDEX_DEN']))

    # Break sample_index_list into equal-sized lists for parallel processing in chunks (last list may be shorter)
    sample_index_chunked = [
        sample_index_list[x:x + SAMPLE_INDEX_CHUNK_SIZE]
            for x in range(0, len(sample_index_list), SAMPLE_INDEX_CHUNK_SIZE)
    ]

    # Setup thread object
    pool = mp.Pool(threads, initializer=init_density)

    # Assign density on sampled sites
    df['STATE'] = -1
    df['KERN_FWD'] = np.nan
    df['KERN_FWDREV'] = np.nan
    df['KERN_REV'] = np.nan
    df['INTERP'] = False

    df.loc[sample_index_list, 'KERN_FWD'] = np.concatenate(pool.map(density_fwd, sample_index_chunked))
    df.loc[sample_index_list, 'KERN_FWDREV'] = np.concatenate(pool.map(density_fwdrev, sample_index_chunked))
    df.loc[sample_index_list, 'KERN_REV'] = np.concatenate(pool.map(density_rev, sample_index_chunked))

    # Get max state
    df.loc[sample_index_list, 'STATE'] = df.loc[
        sample_index_list, ['KERN_FWD', 'KERN_FWDREV', 'KERN_REV']
    ].apply(
        lambda vals: np.argmax(vals.array), axis=1
    )

    # Get indices where density must be computed and indices that can be interpolated
    density_range_list = list()  # Fully compute density
    interp_range_list = list()   # Interpolate from sampled loci

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
            df.loc[range_start, ['KERN_FWD', 'KERN_FWDREV', 'KERN_REV']] -
            df.loc[range_end, ['KERN_FWD', 'KERN_FWDREV', 'KERN_REV']]
        )) > state_run_smooth_delta

        if state_change or density_change:
            # Save indices where density must be calculated
            density_range_list.append(inner_range)

        else:
            # Interpolate densities within range
            interp_range_list.append(inner_range)

    # Interpolate values (where density does not need to be calculated)
    interp_chunk_size = SAMPLE_INDEX_CHUNK_SIZE // state_run_smooth
    interp_index = np.concatenate(interp_range_list)

    df.loc[interp_index, 'KERN_FWD'] = np.concatenate(pool.map(
        interp_fwd, interp_range_list, chunksize=interp_chunk_size
    ))

    df.loc[interp_index, 'KERN_FWDREV'] = np.concatenate(pool.map(
        interp_fwdrev, interp_range_list, chunksize=interp_chunk_size
    ))

    df.loc[interp_index, 'KERN_REV'] = np.concatenate(pool.map(
        interp_rev, interp_range_list, chunksize=interp_chunk_size
    ))

    df.loc[interp_index, 'INTERP'] = True

    # Calculate density values within chunks where values could not be interpolated
    # (too close to potential state changes)
    density_index = np.concatenate(density_range_list)

    density_index_chunked = [
        density_index[x:x + SAMPLE_INDEX_CHUNK_SIZE]
            for x in range(0, len(density_index), SAMPLE_INDEX_CHUNK_SIZE)
    ]

    df.loc[density_index, 'KERN_FWD'] = np.concatenate(pool.map(density_fwd, density_index_chunked))
    df.loc[density_index, 'KERN_FWDREV'] = np.concatenate(pool.map(density_fwdrev, density_index_chunked))
    df.loc[density_index, 'KERN_REV'] = np.concatenate(pool.map(density_rev, density_index_chunked))

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


def rl_encoder(df):
    """
    Take a full table of states and KDE estimates per site and generate a table of contiguous states.

    Recall that the state table has columns removed.

    The table returned has these fields:
    * STATE: State of a contiguous region.
    * POS_KDE: Position in the k-mer table.
    * LEN_KDE
    * POS_QRY
    * LEN_QRY
    * MAX_GAIN

    :param df: Dataframe of states.
    :param pos_qry: Query position at the start of the table (first query base is 'INDEX' + pos_qry).

    :return: DataFrame with STATE, POS_KDE, LEN_KDE, POS_QRY, LEN_QRY, and MAX_GAIN columns.
    """

    # Stop if dataframe is empty
    if df.shape[0] == 0:
        return pd.DataFrame(
            [], index=['STATE', 'POS_KDE', 'LEN_KDE', 'POS_QRY', 'LEN_QRY']
        )

    # Get a table of run-length encoded states (collapse consecutive states)
    df_list = list()

    state = int(df.loc[0, 'STATE'])
    pos = 0

    index = int(df.loc[0, 'INDEX'])

    for i in np.where(df['STATE'][:-1].array != df['STATE'][1:].array)[0] + 1:

        row = df.iloc[i]
        new_index = int(row['INDEX'])

        df_list.append(
            pd.Series(
                [state, pos, i - pos, index, new_index - index],
                index=['STATE', 'POS_KDE', 'LEN_KDE', 'POS_QRY', 'LEN_QRY']
            )
        )

        index = new_index
        pos = i
        state = int(row['STATE'])

    df_list.append(
        pd.Series(
            [state, pos, df.shape[0] - pos, index, df.iloc[-1]['INDEX'] - index],
            index=['STATE', 'POS_KDE', 'LEN_KDE', 'POS_QRY', 'LEN_QRY']
        )
    )

    df_rl = pd.concat(df_list, axis=1).T.astype(int)

    # Score states
    df_rl['MAX_GAIN'] = 0.0

    for index in range(df_rl.shape[0]):
        row = df_rl.loc[index]

        kde_matrix = np.asarray(
            df.loc[
                row['POS_KDE'] :
                row['POS_KDE'] + row['LEN_KDE'],
                ['KDE_FWD', 'KDE_FWDREV', 'KDE_REV']
            ]
        )

        df_rl.loc[index, 'MAX_GAIN'] = np.max(2 * kde_matrix[:, int(row['STATE'])] - kde_matrix.sum(axis=1))

    return df_rl
