"""
Variant call utilities.
"""

import pandas as pd

import pavlib
import svpoplib
import kanapy


class CallerResources(object):
    """
    Container of resources needed by routines attempting to resolve variantns.

    Attributes:
        * df_align: Alignment table.
        * qry_fa: Query FASTA filename.
        * qry_fai: Query FASTA index filename.
        * ref_fa: Reference FASTA filename.
        * ref_fai: Reference FASTA index filename.
        * hap: Haplotype name.
        * score_model: Alignemnt score model.
        * tig_fai: Series of query sequence lengths keyed by query names.
        * align_lift: PAV alignment lift object (translating query/reference coordinates through the assembly).
        * k_util: K-mer utility object.
        * cache_qry_upper: Query sequence cache used for homology searches. Caches the last query sequence in
            upper-case.
        * cache_ref_upper: Reference sequence cache used for homology searches. Caches the last reference sequence in
            upper-case.
    """

    def __init__(self, df_align, df_align_tigref, qry_fa, ref_fa, hap, k_size=31, score_model=None):

        self.df_align = df_align
        self.df_align_tigref = df_align_tigref

        self.qry_fa = qry_fa
        self.qry_fai = self.qry_fa + '.fai'

        self.ref_fa = ref_fa
        self.ref_fai = self.ref_fa + '.fai'

        self.hap = hap

        if k_size is None:
            k_size = 31

        if score_model is None:
            score_model = pavlib.lgsv.align.AffineScoreModel()

        self.score_model = score_model

        self.tig_fai = svpoplib.ref.get_df_fai(self.qry_fai)
        #self.ref_fai = svpoplib.ref.get_df_fai(self.ref_fai)

        self.align_lift = pavlib.align.AlignLift(self.df_align, self.tig_fai)
        self.k_util = kanapy.util.kmer.KmerUtil(k_size)

        self.cache_qry_upper = pavlib.lgsv.SeqCache(self.qry_fa, uppercase=True)  # uppercase for homology searches
        self.cache_ref_upper = pavlib.lgsv.SeqCache(self.ref_fa, uppercase=True)


def record_to_paf(row_seg, ref_fai, qry_fai):
    """
    Convert the row of a segment table to PAF format.

    `row_seg` is a complex segment record with "MAPQ" and "CIGAR" fields added.

    :param row_seg: Segment table row.
    :param ref_fai: Reference FASTA index.
    :param qry_fai: Query FASTA index.

    :return: PAF record row.
    """

    match_n = 0
    align_len = 0
    cigar_index = -1

    cigar_list = list(pavlib.align.cigar_str_to_tuples(row_seg['CIGAR']))

    # Remove clipping and adjust coordinates
    if cigar_list[0][1] == 'H':
        cigar_list = cigar_list[1:]
        cigar_index += 1

    if cigar_list[0][1] == 'S':
        cigar_list = cigar_list[1:]
        cigar_index += 1

    if cigar_list[-1][1] == 'H':
        cigar_list = cigar_list[:-1]

    if cigar_list[-1][1] == 'S':
        cigar_list = cigar_list[:-1]

    cigar = ''.join([f'{op_len}{op_code}' for op_len, op_code in cigar_list])

    # Process CIGAR operations
    for op_len, op_code in cigar_list:
        cigar_index += 1

        if op_code == '=':
            match_n += op_len
            align_len += op_len

        elif op_code in {'X', 'I', 'D'}:
            align_len += op_len

        elif op_code in {'H', 'S'}:
            raise RuntimeError(f'Unhandled clipping in CIGAR string: {op_code} at CIGAR index {cigar_index}: Expected clipped bases at the beginning and end of the CIGAR string only.')

        else:
            raise RuntimeError(f'Unrecognized CIGAR op code: {op_code} at CIGAR index {cigar_index}')

    return pd.Series(
        [
            row_seg['QRY_ID'],
            qry_fai[row_seg['QRY_ID']],
            row_seg['QRY_POS'],
            row_seg['QRY_END'],
            row_seg['STRAND'],
            row_seg['#CHROM'],
            ref_fai[row_seg['#CHROM']],
            row_seg['POS'],
            row_seg['END'],
            match_n,
            align_len,
            row_seg['MAPQ'],
            cigar
        ],
        index=[
            'QRY_NAME',
            'QRY_LEN',
            'QRY_POS',
            'QRY_END',
            'STRAND',
            'CHROM',
            'CHROM_LEN',
            'CHROM_POS',
            'CHROM_END',
            'MISMATCH_N',
            'ALIGN_BLK_LEN',
            'MAPQ',
            'CIGAR'
        ]
    )
