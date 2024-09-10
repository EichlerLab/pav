"""
General utility functions.
"""

import numpy as np
import os
import pandas as pd


def as_bool(val, fail_to_none=False):
    """
    Translate value as a boolean. If `val` is boolean, return `val`. If val is a string, `True` if lower-case string
    is "true", "1", "yes", "t", or "y". `False` if lower-case string is "false", "0", "no", "f", or "n". All other
    values throw a `RuntimeError` unless `fail_to_none` is `True`. If `val` is not bool or string, it is converted to a
    string and the rules above are applied.

    :param val: Value to interpret as a boolean.
    :param fail_to_none: If `True`, return `None` if `val` is not a recognized boolean value.

    :return: `True` or `False` (see above).
    """

    if issubclass(val.__class__, bool):
        return val

    val = str(val).lower()

    if val in {'true', '1', 'yes', 't', 'y'}:

        return True

    if val in {'false', '0', 'no', 'f', 'n'}:
        return False

    if fail_to_none:
        return None

    raise RuntimeError('Cannot interpret as boolean value: {}'.format(val))


def region_merge(file_list, pad=500):
    """
    Merge regions from multiple BED files.

    :param file_list: List of files to merge. Each element may be a string (filename) or a DataFrame (in-memory table
        of alignments). Files with zero size are skipped.
    :param pad: Pad interval matches by this amount, but do not add it to the output intervals. Similar to bedtools
        merge "slop" parameter.
    """

    # Get BED list
    bed_list = list()

    for file_name in file_list:
        if isinstance(file_name, str):
            if os.stat(file_name).st_size > 0:  # Skip if file is empty
                bed_list.append(
                    pd.read_csv(
                        file_name, sep='\t', usecols=('#CHROM', 'POS', 'END')
                    )
                )

        elif isinstance(file_name, pd.DataFrame):
            bed_list.append(file_name)

        else:
            raise RuntimeError(f'File name is not a string or a DataFrame: {file_name} (type "{type(file_name)}")')

    df = pd.concat(bed_list, axis=0).sort_values(['#CHROM', 'POS', 'END'], ascending=[True, True, False]).reset_index(drop=True)

    # # Read regions
    # df = pd.concat(
    #     [
    #         pd.read_csv(
    #             file_name,
    #             sep='\t',
    #             usecols=('#CHROM', 'POS', 'END')
    #         ) for file_name in file_list if os.stat(file_name).st_size > 0
    #     ],
    #     axis=0
    # ).sort_values(['#CHROM', 'POS', 'END'], ascending=[True, True, False]).reset_index(drop=True)

    # Merge with intervals
    df_list = list()

    chrom = None
    pos = None
    end = None

    for index, row in df.iterrows():

        next_chrom = row['#CHROM']
        next_pos = row['POS'] - 500
        next_end = row['END'] + 500

        if row['#CHROM'] != chrom:

            # Write last record
            if chrom is not None:
                df_list.append(pd.Series([chrom, pos + pad, end - pad], index=['#CHROM', 'POS', 'END']))

            chrom = next_chrom
            pos = next_pos
            end = next_end

        else:

            if next_pos <= end:
                pos = np.min([pos, next_pos])
                end = np.max([end, next_end])

            else:
                df_list.append(pd.Series([chrom, pos + pad, end - pad], index=['#CHROM', 'POS', 'END']))

                pos = next_pos
                end = next_end

    # Write last record
    if chrom is not None:
        df_list.append(pd.Series([chrom, pos + pad, end - pad], index=['#CHROM', 'POS', 'END']))

    # Return
    if len(df_list) > 0:
        return pd.concat(df_list, axis=1).T
    else:
        return pd.DataFrame([], columns=['#CHROM', 'POS', 'END'])


def collapse_to_set(l, to_type=None):
    l = list(l)  # Copy so the original list is not modified
    s = set()

    while len(l) > 0:
        v = l.pop()

        if isinstance(v, (tuple, list)):
            l.extend(v)
        else:
            if to_type is not None:
                s.add(to_type(v))
            else:
                s.add(v)

    return s
