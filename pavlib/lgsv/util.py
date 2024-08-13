"""
Variant call utilities.
"""

import pandas as pd

import Bio.Seq
import pavlib
import svpoplib


class CallerResources(object):
    """
    Container of resources needed by routines attempting to resolve variantns.

    Attributes:
        * df_align_qry: Alignment table (QRY trimmed)
        * df_align_qryref: Alignment table (QRY & REF trimmed)
        * qry_fa: Query FASTA filename.
        * qry_fai: Query FASTA index filename.
        * ref_fa: Reference FASTA filename.
        * ref_fai: Reference FASTA index filename.
        * hap: Haplotype name.
        * score_model: Alignemnt score model.
        * df_qry_fai: Query lengths keyed by query ID.
        * align_lift: Object for lifting alignment coordinates between query and reference through the alignment.
        * cache_qry_upper: Query sequence cache used for homology searches. Caches the last query sequence in
            upper-case.
        * cache_ref_upper: Reference sequence cache used for homology searches. Caches the last reference sequence in
            upper-case.
    """

    def __init__(self, df_align_qry, df_align_qryref, qry_fa, ref_fa, hap, score_model=None):

        self.df_align_qry = df_align_qry
        self.df_align_qryref = df_align_qryref

        self.qry_fa = qry_fa
        self.qry_fai = self.qry_fa + '.fai'

        self.ref_fa = ref_fa
        self.ref_fai = self.ref_fa + '.fai'

        self.hap = hap

        if score_model is None:
            score_model = pavlib.align.score.get_score_model()

        self.score_model = score_model

        self.df_qry_fai = svpoplib.ref.get_df_fai(self.qry_fai)
        #self.ref_fai = svpoplib.ref.get_df_fai(self.ref_fai)

        self.align_lift = pavlib.align.lift.AlignLift(self.df_align_qry, self.df_qry_fai)

        self.cache_qry_upper = pavlib.lgsv.util.SeqCache(self.qry_fa, uppercase=True)  # uppercase for homology searches
        self.cache_ref_upper = pavlib.lgsv.util.SeqCache(self.ref_fa, uppercase=True)


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

    cigar_list = list(pavlib.align.util.cigar_str_to_tuples(row_seg['CIGAR']))

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

class SeqCache:
    """
    Keep a cache of a sequence string in upper-case. Stores the last instance of the sequence and the ID. When a
    new ID is requested, the old sequnece is discarded and the new one is loaded.
    """

    def __init__(self, fa_filename, uppercase=True):
        """
        Create a cache object to read from indexed FASTA file `fa_filename`.

        :param fa_filename: Indexed FASTA file name.
        :param uppercase: `True` if sequences should be made upper-case, otherwise, preserve case.
        """

        self.fa_filename = fa_filename
        self.uppercase = uppercase

        self.id = None
        self.seq = None
        self.is_rev = None

    def get(self, sequence_id, is_rev):
        """
        Get a sequence. Returns the cached version if ID matches the current ID, otherwise, the correct sequence is
        retrieved, cached, and returned.

        :param sequence_id: Sequence ID string or Region.
        :param is_rev: `True` if the sequence is reverse-complemented. Retrieving the same sequence ID as the cached
            sequence with `is_rev` mismatch will reload the sequence in the requested orientation.

        :return: Sequence.
        """

        if self.id != sequence_id or self.is_rev != is_rev:
            new_seq = pavlib.seq.region_seq_fasta(
                sequence_id, self.fa_filename
            )

            if self.uppercase:
                new_seq = new_seq.upper()

            if is_rev:
                new_seq = str(Bio.Seq.Seq(new_seq).reverse_complement())

            self.seq = new_seq

            self.id = sequence_id
            self.is_rev = is_rev

        return self.seq
