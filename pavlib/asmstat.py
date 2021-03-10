"""
Generate assembly stats.
"""

import numpy as np
import pandas as pd

from Bio import SeqIO

import svpoplib


def get_n_stat(len_list, n_stat=0.5, genome_size=None):
    """
    Get N statistic (e.g. N50 and NG50) from an array of contig lengths.

    :param vals: Numpy array of contig lengths. Does not need to be sorted.
    :param n_stat: N stat (e.g. "0.5" is N50, "0.1" is N10, etc).
    :param genome_size: Size of genome for NG50 calculation. If `None`, N50 using the assembly size is returned.

    :return: N50.
    """

    vals = pd.Series(len_list).sort_values(ascending=False)

    vals_csum = np.cumsum(vals)

    if genome_size is None:
        genome_size = np.int32(np.sum(vals) * n_stat)
    else:
        genome_size = np.int32(genome_size * 0.5)

    return vals.iloc[np.sum(vals_csum <= genome_size) + 1]


def get_stats(asm_name, hap, fa_file_name, genome_size=None, n_stat_list=[0.5]):
    """
    Get assembly stats.

    :param asm_name: Assembly name.
    :param hap: Haplotype.
    :param fa_file_name: Assembly FASTA file name.
    :param genome_size: If defined, add an NG50 field for a genome of this size.

    :return: A Pandas series with assembly stats.
    """

    # Get lengths

    with svpoplib.seq.PlainOrGzReader(fa_file_name) as in_file:
        len_list = [len(record.seq) for record in SeqIO.parse(in_file, 'fasta')]

    # Begin record
    record_list = [
        asm_name,
        hap,
        np.sum(len_list) / 1e9,
        len(len_list)
    ]

    index_list = ['ASM_NAME', 'HAP', 'GBP', 'CONTIGS']

    # Add N50
    for n_stat in n_stat_list:
        record_list.append(get_n_stat(len_list, n_stat=n_stat, genome_size=None) / 1e6)
        index_list.append('N{:d}_MBP'.format(np.int32(n_stat * 100)))

    if genome_size is not None:
        for n_stat in n_stat_list:
            record_list.append(get_n_stat(len_list, n_stat=n_stat, genome_size=genome_size) / 1e6)
            index_list.append('NG{:d}_MBP'.format(np.int32(n_stat * 100)))

    record_list.append(np.sum([contig_len for contig_len in len_list if contig_len >= 100000]))
    index_list.append('100k+')

    record_list.append(np.sum([contig_len for contig_len in len_list if contig_len >= 1000000]))
    index_list.append('1M+')

    # Return record
    return pd.Series(record_list, index_list)
