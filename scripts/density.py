"""
K-mer density used for calling inversions.
"""

import argparse
import codecs
import multiprocessing as mp
import numpy as np
import os
import pickle
import pandas as pd
import sys

import scipy.stats

# Add PAV libraries and dependencies
PIPELINE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

sys.path.append(PIPELINE_DIR)  # pavlib
sys.path.append(os.path.join(PIPELINE_DIR, 'dep', 'svpop'))  # svpoplib
sys.path.append(os.path.join(PIPELINE_DIR, 'dep'))  # kanapy

import pavlib
import svpoplib
import kanapy



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

SAMPLE_INDEX_CHUNK_SIZE = 400  # Chunk density indices into groups of this many and compute density together.

MAX_REF_KMER_COUNT = 100   # Skip low-complexity regions


#
# Objects used by threads
#

kernel_dict = None

sum_state_fwd = None
sum_state_fwdrev = None
sum_state_rev = None

df = None

density_bandwidth = None


#
# Init Function
#

def init_process():

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


#
# Density calculation functions
#

# Direct density
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

    return np.interp(
        inner_range,
        [range_start, range_end],
        [df.loc[range_start, 'KERN_FWD'], df.loc[range_end, 'KERN_FWD']]
    )


def interp_fwdrev(inner_range):
    range_start = inner_range[0] - 1
    range_end = inner_range[-1] + 1

    return np.interp(
        inner_range,
        [range_start, range_end],
        [df.loc[range_start, 'KERN_FWDREV'], df.loc[range_end, 'KERN_FWDREV']]
    )


def interp_rev(inner_range):
    range_start = inner_range[0] - 1
    range_end = inner_range[-1] + 1

    return np.interp(
        inner_range,
        [range_start, range_end],
        [df.loc[range_start, 'KERN_REV'], df.loc[range_end, 'KERN_REV']]
    )


def get_smoothed_density():

    global df
    global density_bandwidth

    global sum_state_fwd
    global sum_state_fwdrev
    global sum_state_rev


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

        df = df.loc[df['STATE_MER'] != low_state].copy()

    # Check for number of informative k-mers before computing density
    if df.shape[0] < min_informative_kmers:
        return
        #return pd.DataFrame([], columns=['INDEX', 'STATE_MER', 'STATE', 'KERN_FWD', 'KERN_FWDREV', 'KERN_REV', 'KMER'])

    # Setup bandwidth and index in informative k-mer space (ignore df['INDEX'] for density)
    density_bandwidth = df.shape[0] ** (-1.0 / 5.0) * density_smooth_factor

    # Index in condensed space
    df.insert(1, 'INDEX_DEN', np.arange(df.shape[0]))

    df.set_index('INDEX_DEN', inplace=True, drop=False)  # Switch index to density space

    # Count states for scaled-density
    sum_state_fwd = np.sum(df['STATE_MER'] == 0)
    sum_state_fwdrev = np.sum(df['STATE_MER'] == 1)
    sum_state_rev = np.sum(df['STATE_MER'] == 2)

    # Setup indexes (initial non-sampled sites density is initially calcuated over)
    sample_index_list = list(df.loc[df['INDEX_DEN'] % state_run_smooth == 0, 'INDEX_DEN'])

    if sample_index_list[-1] != np.int32(df.iloc[-1]['INDEX_DEN']):  # Add last table element if it does not exist
        sample_index_list.append(np.int32(df.iloc[-1]['INDEX_DEN']))

    # Break sample_index_list into equal-sized lists for parallel processing in chunks (last list may be shorter)
    sample_index_chunked = [
        sample_index_list[x:x + SAMPLE_INDEX_CHUNK_SIZE]
            for x in range(0, len(sample_index_list), SAMPLE_INDEX_CHUNK_SIZE)
    ]

    # Init fields
    df['STATE'] = -1
    df['KERN_FWD'] = np.nan
    df['KERN_FWDREV'] = np.nan
    df['KERN_REV'] = np.nan
    df['INTERP'] = False

    # Setup thread object
    pool = mp.Pool(threads, initializer=init_process)

    # Assign density on samples sites
    df.loc[sample_index_list, 'KERN_FWD'] = np.concatenate(pool.map(density_fwd, sample_index_chunked))
    df.loc[sample_index_list, 'KERN_FWDREV'] = np.concatenate(pool.map(density_fwdrev, sample_index_chunked))
    df.loc[sample_index_list, 'KERN_REV'] = np.concatenate(pool.map(density_rev, sample_index_chunked))

    # Destroy thread object (must be re-created to see changes to df)
    pool.close()
    del(pool)

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

    # Create pool object (with updated df in their globals)
    pool = mp.Pool(threads, initializer=init_process)

    # Interpolate values (where density does not need to be calculated)
    if len(interp_range_list) > 0:

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
    if len(density_range_list) > 0:
        density_index = np.concatenate(density_range_list)

        density_index_chunked = [
            density_index[x:x + SAMPLE_INDEX_CHUNK_SIZE]
                for x in range(0, len(density_index), SAMPLE_INDEX_CHUNK_SIZE)
        ]

        df.loc[density_index, 'KERN_FWD'] = np.concatenate(pool.map(density_fwd, density_index_chunked))
        df.loc[density_index, 'KERN_FWDREV'] = np.concatenate(pool.map(density_fwdrev, density_index_chunked))
        df.loc[density_index, 'KERN_REV'] = np.concatenate(pool.map(density_rev, density_index_chunked))

    # Close pool object and free threads
    pool.close()
    del(pool)

    # Penalize spikes above 1.0.
    df.loc[df['KERN_FWD'] > 1.0, 'KERN_FWD'] = 1 / df.loc[df['KERN_FWD'] > 1.0]
    df.loc[df['KERN_FWDREV'] > 1.0, 'KERN_FWDREV'] = 1 / df.loc[df['KERN_FWDREV'] > 1.0]
    df.loc[df['KERN_REV'] > 1.0, 'KERN_REV'] = 1 / df.loc[df['KERN_REV'] > 1.0]

    # Set max densities on all
    df['STATE'] = df[['KERN_FWD', 'KERN_FWDREV', 'KERN_REV']].apply(
        lambda vals: np.argmax(vals.array),
        axis=1
    )

    # Column order and index
    df = df[['INDEX', 'STATE_MER', 'STATE', 'KERN_FWD', 'KERN_FWDREV', 'KERN_REV', 'KMER']]
    df.set_index(df['INDEX'], inplace=True, drop=False)


# def do_error(err_msg, ex_func=RuntimeError):
#     """
#     Handle an error. If stdout mode (no output files), then return exception to stdout as a PKL string. Otherwise,
#     raise the exception. The program is terminated either way.
#
#     :param err_msg: Error message.
#     :param ex_func: Exception constructor.
#     """
#
#     raise ex = ex_func(err_msg)
#
#     if do_stdout:
#         sys.stderr.write('Writing exception to stdout: \n')
#         sys.stderr.write(str(ex))
#         sys.stderr.write('\n')
#
#         sys.exit(1)
#
#     else:
#         raise ex


def write_table(out_file_name, test=False):
    """
    Write output dataframe (global df) to an output file. Recognized TSV (.tsv, .tsv.gz) or Excel (.xlsx). Raises an
    exception if the output file type cannot be determined.

    :param out_file_name: Output file name.
    :param test: If `True`, do not write, just check for known file names. Useful for testing output arguments before
        generating the dataframe.
    """

    out_lower = out_file_name.lower()

    if out_lower.endswith('.tsv'):
        if not test:
            df.to_csv(out_file_name, sep='\t', index=False)

    elif out_lower.endswith('.tsv.gz'):
        if not test:
            df.to_csv(out_file_name, sep='\t', index=False, compression='gzip')
    elif out_lower.endswith('.xlsx'):
        if not test:
            df.to_excel(out_file_name, index=False)

    else:
        raise RuntimeError(f'No recognized extension on output file: {out_file_name}: Expected tsv, tsv.gz, or xlsx')


def get_bool(bool_str):
    """
    Get boolean from string.

    * True: "true", "t", "1"
    * False: "false", "f", "0"

    Values are not case sensitive. All other values raise `RuntimeError`.

    :param bool_str: Boolean string.

    :return: Boolean `True` or `False`.
    """

    bs_lower = bool_str.lower()

    if bs_lower in {'true', 't', '1', 1, True}:
        return True

    if bs_lower in {'false', 'f', '0', 0, True}:
        return False

    raise RuntimeError(f'Unrecognized boolean string value: {bool_str}')


#
# Main
#

if __name__ == '__main__':

    # Get command-line arguments
    parser = argparse.ArgumentParser('Inversion density calculation')

    parser.add_argument('--tigregion', help='Contig region to extract.')

    parser.add_argument('--refregion', help='Reference region to extract.')

    parser.add_argument('--ref', help='Reference FASTA file')

    parser.add_argument('--tig', help='Contig FASTA file')

    parser.add_argument('-k', type=int, default=31, help='K-mer size')

    parser.add_argument('-t', '--threads', default=1, type=int, help='Number of concurrent threads')

    parser.add_argument('-r', '--revcompl',
                        help='Reverse-complement reference k-mers. Set this if the region aligned was '
                             'revese-complemented. If the inversion flanks are in reference orientation, then this '
                             'value should be "true". If the inversion flanks are not in reference orientation, then '
                             'this value should be "true". True values are "true", "t" and "1"; False values are '
                             '"false", "f", and 0; none are case sensitive. All other argument values raise an '
                             'execption.'
                        )

    parser.add_argument('--mininf', type=int, default=2000,
                        help='Do not attempt density if the number of informative k-mers does not reach this limit. '
                             'Informative k-mers are defined as k-mers that are in forward and/or reverse-complement '
                             'in the reference k-mer set.'
                        )

    parser.add_argument('--densmooth', type=int, default=1,
                        help='Smooth density by this factor. Density bandwidth is estimated by Scott\'s rule then '
                             'multiplied by this factor to further smooth (> 1) or restrict smoothing (< 1). See '
                             'https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html'
                        )

    parser.add_argument('--minstatecount', type=int, default=20,
                        help='States (fwd, fwd-rev, rev) that do not appear in at least this many k-mers are '
                             'discarded. This prevents extreme density spikes concentrated in a small number of '
                             'k-mers from states that are too sparse to be informative, which complicates downstream '
                             'analysis.')

    parser.add_argument('--staterunsmooth', type=int, default=20,
                        help='Smooth density calculations by this many records. Initially, calculate density only '
                             'once in a run of this many k-mers. If there is a state change or a delta larger than '
                             '"--staterundelta", then density is filled in for all skipped elements. If there '
                             'is no state change and the delta is small, use linear interpolation to fill in missing '
                             'values.'
                        )

    parser.add_argument('--staterundelta', type=float, default=0.005,
                        help='Changes between state densities by this much or more will be filled in with actual '
                             'density values instead of interpolated. See "--staterunsmooth".')

    parser.add_argument('outfile', nargs='*',
                        help='PKL (.pkl), TSV (.tsv, .tsv.gz), or excel (.xlsx) file containing the Pandas DataFrame '
                             'of fwd/fwd-rev/rev k-mer density information. File type is determined by extension. '
                             'Multiple files may be written if more than one output file is specified. Send PKL to '
                             'stdout as a base64 encoded string if no output files are specified.'
                        )

    args = parser.parse_args()

    do_stdout = (len(args.outfile) == 0)

    is_rev = get_bool(args.revcompl)


    # Test output file names before generating the dataframe.
    if not do_stdout:
        for out_file_name in args.outfile:
            write_table(out_file_name, test=True)


    ### Get regions and utilities ###

    region_ref = pavlib.seq.region_from_string(args.refregion)
    region_tig = pavlib.seq.region_from_string(args.tigregion)

    k_util = kanapy.util.kmer.KmerUtil(args.k)


    ### Get reference k-mer counts ###
    ref_kmer_count = pavlib.seq.ref_kmers(region_ref, args.ref, k_util)

    if ref_kmer_count is None or len(ref_kmer_count) == 0:
        raise RuntimeError(f'No reference k-mers for region {region_ref}')

    # Skip low-complexity sites with repetitive k-mers
    max_mer_count = np.max(list(ref_kmer_count.values()))

    if max_mer_count > MAX_REF_KMER_COUNT:
        max_mer = [kmer for kmer, count in ref_kmer_count.items() if count == max_mer_count][0]

        raise RuntimeError('K-mer count exceeds max: {} > {} ({}): {}'.format(
            max_mer_count,
            MAX_REF_KMER_COUNT,
            k_util.to_string(max_mer),
            region_ref
        ))

    ref_kmer_set = set(ref_kmer_count)

    if is_rev:
        ref_kmer_set = {k_util.rev_complement(kmer) for kmer in ref_kmer_set}

    ## Get contig k-mers as list ##

    seq_tig = pavlib.seq.region_seq_fasta(region_tig, args.tig)

    tig_mer_stream = list(kanapy.util.kmer.stream(seq_tig, k_util, index=True))

    ### Make args global variable (for convenience) ###
    threads = args.threads
    min_informative_kmers = args.mininf
    density_smooth_factor = args.densmooth
    min_state_count = args.minstatecount
    state_run_smooth = args.staterunsmooth
    state_run_smooth_delta = args.staterundelta


    ### Set global df with Density table ###
    get_smoothed_density()

    # Write
    if do_stdout:
        sys.stdout.write(
            codecs.encode(
                pickle.dumps(df),
                'base64'
            ).decode()
        )

    else:

        for out_file_name in args.outfile:
            write_table(out_file_name)
