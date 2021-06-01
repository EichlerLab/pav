"""
Routines for handling alignments.
"""

import collections
import intervaltree
import numpy as np
import pandas as pd
import pysam

import pavlib.seq

import svpoplib.ref

_INT_STR_SET = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}
_CIGAR_OP_SET = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}

# Indices for tuples returned by trace_cigar_to_zero()
TC_INDEX = 0
TC_OP_LEN = 1
TC_OP_CODE = 2
TC_DIFF_CUM = 3
TC_DIFF = 4
TC_EVENT_CUM = 5
TC_EVENT = 6
TC_SUB_BP = 7
TC_QRY_BP = 8
TC_CLIPS_BP = 9
TC_CLIPH_BP = 10


def trim_alignments(df, min_trim_tig_len, tig_fai):
    """
    Do alignment trimming from prepared alignment BED file. This BED contains information about reference and
    contig coordinates mapped, CIGAR string and flags.

    :param df: Alignment dataframe.
    :param min_trim_tig_len: Minimum alignment record length. Alignment records smaller will be discarded.
    :param tig_fai: FAI file from the contig FASTA that was mapped. Used for an alignment sanity check after trimming.

    :return: Trimmed alignments as an alignment DataFrame. Same format as `df` with columns added describing the
        number of reference and contig bases that were trimmed. Dropped records (mapped inside another or too shart) are
        removed.
    """

    # Add fields for the number of bases that are removed from each end
    df['CUT_REF_L'] = 0
    df['CUT_REF_R'] = 0
    df['CUT_TIG_L'] = 0
    df['CUT_TIG_R'] = 0

    # Remove short alignments
    for index in df.index:
        if df.loc[index, 'QUERY_TIG_END'] - df.loc[index, 'QUERY_TIG_POS'] < min_trim_tig_len:
            df.loc[index, 'INDEX'] = -1


    ###                                          ###
    ### Trim overlapping contigs in contig space ###
    ###                                          ###

    # Discard fully trimmed records
    df = df.loc[df['INDEX'] >= 0].copy()

    # Sort by alignment lengths in contig space
    df['QUERY_LEN'] = df['QUERY_END'] - df['QUERY_POS']
    df['SUB_LEN'] = df['END'] - df['POS']

    df.sort_values(['QUERY_ID', 'QUERY_LEN'], ascending=(True, False), inplace=True)

    df.reset_index(inplace=True, drop=True)

    # Do trim in contig space
    iter_index_l = 0
    index_max = df.shape[0]

    while iter_index_l < index_max:

        iter_index_r = iter_index_l + 1

        while iter_index_r < index_max and df.loc[iter_index_l, 'QUERY_ID'] == df.loc[iter_index_r, 'QUERY_ID']:

            # Get index in order of contig placement
            if df.loc[iter_index_l, 'QUERY_TIG_POS'] <= df.loc[iter_index_r, 'QUERY_TIG_POS']:
                index_l = iter_index_l
                index_r = iter_index_r
            else:
                index_l = iter_index_r
                index_r = iter_index_l

            # Skip if one record was already removed
            if df.loc[index_l, 'INDEX'] < 0 or df.loc[index_r, 'INDEX'] < 0:
                iter_index_r += 1
                continue

            # Skip if there is no overlap
            if df.loc[index_r, 'QUERY_TIG_POS'] >= df.loc[index_l, 'QUERY_TIG_END']:
                iter_index_r += 1
                continue

            # Check for record fully contained within another
            if df.loc[index_r, 'QUERY_TIG_END'] <= df.loc[index_l, 'QUERY_TIG_END']:
                df.loc[index_r, 'INDEX'] = -1
                iter_index_r += 1
                continue

            # Determine trim orientation (right side of index_l is to trimmed, so must be reversed so
            # trimmed CIGAR records are at the beginning; left side of index_r is to be trimmed, which is
            # already at the start of the CIGAR string).
            rev_l = not df.loc[index_l, 'REV']  # Trim right end of index_l
            rev_r = df.loc[index_r, 'REV']      # Trim left end of index_r

            # Detect overlap in contig space
            if rev_l == rev_r or df.loc[index_l, '#CHROM'] != df.loc[index_r, '#CHROM']:
                # Contigs were aligned in different reference chromosomes or orientations, no overlap
                ref_overlap = False

            else:
                if df.loc[index_l, 'POS'] < df.loc[index_r, 'POS']:
                    ref_overlap = df.loc[index_r, 'POS'] < df.loc[index_r, 'END']

                elif df.loc[index_r, 'POS'] < df.loc[index_l, 'POS']:
                    ref_overlap = df.loc[index_l, 'POS'] < df.loc[index_r, 'END']

                else:
                    # POS placement was the same,
                    ref_overlap = False

            # If contig alignments overlap in reference space and are in the same orientation, preferentially
            # trim the downstream end more to left-align alignment-truncating SVs (first argument to
            # trim_alignment_record is preferentially trimmed)
            if ref_overlap:
                # Try both trim orientations

                # a: Try with record l as first and r as second
                record_l_a, record_r_a = trim_alignment_record(
                    df.loc[index_l], df.loc[index_r], 'query',
                    rev_l=rev_l,
                    rev_r=rev_r
                )

                # b: Try with record r as first and l as second
                record_l_b, record_r_b = trim_alignment_record(
                    df.loc[index_r], df.loc[index_l], 'query',
                    rev_l=rev_r,
                    rev_r=rev_l
                )

                ### Determine which left-aligns best ###
                keep = None

                # Case: Alignment trimming completely removes one of the records
                rm_l_a = record_l_a['QUERY_TIG_END'] - record_l_a['QUERY_TIG_POS'] < min_trim_tig_len
                rm_l_b = record_l_b['QUERY_TIG_END'] - record_l_b['QUERY_TIG_POS'] < min_trim_tig_len

                rm_r_a = record_r_a['QUERY_TIG_END'] - record_r_a['QUERY_TIG_POS'] < min_trim_tig_len
                rm_r_b = record_r_b['QUERY_TIG_END'] - record_r_b['QUERY_TIG_POS'] < min_trim_tig_len

                rm_any_a = rm_l_a or rm_r_a
                rm_any_b = rm_l_b or rm_r_b

                # Break tie if one way removes a record and the other does not
                if rm_any_a and not rm_any_b:
                    if not rm_l_a and rm_r_a:
                        keep = 'a'

                elif rm_any_b and not rm_any_a:
                    if not rm_l_b and rm_r_b:
                        keep = 'b'

                # Break tie if both are removed in one (do not leave short alignments).
                if keep is None and rm_any_a:  # Both l and r are None, case where one is None was checked
                    keep = 'a'

                if keep is None and rm_any_b:  # Both l and r are None, case where one is None was checked
                    keep = 'b'

                # Break tie on most left-aligned base
                if keep is None:

                    # Get position at end of trim
                    trim_pos_l_a = record_l_a['END'] if not record_l_a['REV'] else record_l_a['POS']
                    trim_pos_l_b = record_l_b['END'] if not record_l_b['REV'] else record_l_b['POS']

                    if trim_pos_l_a <= trim_pos_l_b:
                        keep = 'a'
                    else:
                        keep = 'b'

                # Set record_l and record_r to the kept record
                if keep == 'a':
                    record_l = record_l_a
                    record_r = record_r_a

                else:
                    # Note: record at index_l became record_r_b (index_r become record_l_b)
                    # Swap back to match the index
                    record_l = record_r_b
                    record_r = record_l_b

            else:

                # Switch record order if they are on the same contig and same orientation (note: rev_l and rev_r are
                # opposite if contigs are mapped in the same orientation, one was swapped to trim the dowstream-
                # aligned end).
                if df.loc[index_l, '#CHROM'] == df.loc[index_r, '#CHROM'] and rev_l != rev_r:

                    # Get position of end to be trimmed
                    trim_pos_l = df.loc[index_l, 'END'] if not df.loc[index_l, 'REV'] else df.loc[index_l, 'POS']
                    trim_pos_r = df.loc[index_r, 'POS'] if not df.loc[index_r, 'REV'] else df.loc[index_r, 'END']

                    # Swap positions so the upstream-aligned end of the contig is index_l. The left end is
                    # preferentially trimmed shorter where there are equal breakpoints effectively left-aligning
                    # around large SVs (e.g. large DELs).
                    if trim_pos_r < trim_pos_l:
                        # Swap
                        rev_tmp = rev_l
                        rev_l = rev_r
                        rev_r = rev_tmp

                        index_tmp = index_l
                        index_l = index_r
                        index_r = index_tmp

                # Trim record
                record_l, record_r = trim_alignment_record(
                    df.loc[index_l], df.loc[index_r], 'query',
                    rev_l=rev_l,
                    rev_r=rev_r
                )

            # Modify if new aligned size is at least min_trim_tig_len, remove if shorter
            if record_l['QUERY_TIG_END'] - record_l['QUERY_TIG_POS'] >= min_trim_tig_len:
                df.loc[index_l] = record_l
            else:
                df.loc[index_l, 'INDEX'] = -1

            if (record_r['QUERY_TIG_END'] - record_r['QUERY_TIG_POS']) >= min_trim_tig_len:
                df.loc[index_r] = record_r
            else:
                df.loc[index_r, 'INDEX'] = -1

            # Next r record
            iter_index_r += 1

        # Next l record
        iter_index_l += 1

    ###                                             ###
    ### Trim overlapping contigs in reference space ###
    ###                                             ###

    # Discard fully trimmed records
    df = df.loc[df['INDEX'] >= 0].copy()

    # Sort by contig alignment length in reference space
    df['QUERY_LEN'] = df['QUERY_END'] - df['QUERY_POS']
    df['SUB_LEN'] = df['END'] - df['POS']

    df.sort_values(['#CHROM', 'SUB_LEN'], ascending=(True, False), inplace=True)

    df.reset_index(inplace=True, drop=True)

    # Do trim in reference space
    iter_index_l = 0
    index_max = df.shape[0]

    while iter_index_l < index_max:
        iter_index_r = iter_index_l + 1

        while (
                iter_index_r < index_max and
                df.loc[iter_index_l, '#CHROM'] == df.loc[iter_index_r, '#CHROM']
        ):

            # Skip if one record was already removed
            if df.loc[iter_index_l, 'INDEX'] < 0 or df.loc[iter_index_r, 'INDEX'] < 0:
                iter_index_r += 1
                continue

            # Get indices ordered by contig placement
            if df.loc[iter_index_l, 'POS'] <= df.loc[iter_index_r, 'POS']:
                index_l = iter_index_l
                index_r = iter_index_r
            else:
                index_l = iter_index_r
                index_r = iter_index_l

            # Check for overlaps
            if df.loc[index_r, 'POS'] < df.loc[index_l, 'END']:
                # Found overlapping records
                # print('Ref Overlap: {}-{} ({}:{}-{},{} vs {}:{}-{},{}) [iter {}, {}]'.format(
                #     df.loc[index_l, 'INDEX'], df.loc[index_r, 'INDEX'],
                #     df.loc[index_l, 'QUERY_ID'], df.loc[index_l, 'QUERY_TIG_POS'], df.loc[index_l, 'QUERY_TIG_END'], ('-' if df.loc[index_l, 'REV'] else '+'),
                #     df.loc[index_r, 'QUERY_ID'], df.loc[index_r, 'QUERY_TIG_POS'], df.loc[index_r, 'QUERY_TIG_END'], ('-' if df.loc[index_r, 'REV'] else '+'),
                #     iter_index_l, iter_index_r
                # ))

                # Check for record fully contained within another
                if df.loc[index_r, 'END'] <= df.loc[index_l, 'END']:
                    # print('\t* Fully contained')

                    df.loc[index_r, 'INDEX'] = -1

                else:

                    record_l, record_r = pavlib.align.trim_alignment_record(df.loc[index_l], df.loc[index_r], 'subject')

                    if record_l is not None and record_r is not None:

                        # Modify if new aligned size is at least min_trim_tig_len, remove if shorter
                        if record_l['QUERY_TIG_END'] - record_l['QUERY_TIG_POS'] >= min_trim_tig_len:
                            df.loc[index_l] = record_l
                        else:
                            df.loc[index_l, 'INDEX'] = -1

                        if (record_r['QUERY_TIG_END'] - record_r['QUERY_TIG_POS']) >= min_trim_tig_len:
                            df.loc[index_r] = record_r
                        else:
                            df.loc[index_r, 'INDEX'] = -1

            # Next r record
            iter_index_r += 1

        # Next l record
        iter_index_l += 1

    ###                      ###
    ### Post trim formatting ###
    ###                      ###

    # Clean and re-sort
    df = df.loc[df['INDEX'] >= 0].copy()

    df = df.loc[(df['END'] - df['POS']) > 0]  # Should never occur, but don't allow 0-length records
    df = df.loc[(df['QUERY_END'] - df['QUERY_POS']) > 0]

    df.sort_values(['#CHROM', 'POS', 'END', 'QUERY_ID'], ascending=[True, True, False, True], inplace=True)

    del(df['QUERY_LEN'])
    del(df['SUB_LEN'])

    # Check sanity
    df_tig_fai = svpoplib.ref.get_df_fai(tig_fai)

    df.apply(pavlib.align.check_record, df_tig_fai=df_tig_fai, axis=1)

    # Return trimmed alignments
    return df


def trim_alignment_record(record_l, record_r, match_coord, rev_l=True, rev_r=False):
    """
    Trim ends of overlapping alignments until ends no longer overlap. In repeat-mediated events, aligners (e.g.
    minimap2) will align the same parts of a contig to both reference copies (e.g. large DEL) or two parts of a contig
    to the same region (e.g. tandem duplication). This function trims back the alignments using the CIGAR string until
    the overlap is resolved using a simple greedy algorithm that minimizes variation.

    For example, a large repeat-mediated deletion will have two reference copies, but one copy in the contig, and the
    single contig copy is aligned to both by breaking the alignment record into two (one up to the deletion, and one
    following it). If the contig coordinates were ignored, the alignment gap is smaller than the actual deletion event
    and one or both sides of the deletion are filled with false variants. In this example, the alignment is walked-
    out from both ends of the deletion until there is no duplication of aligned contig (e.g. the alignment stops at
    one contig base and picks up at the next contig base). In this case, this function would be asked to resolve the
    contig coordinates (match_coord = "query").

    A similar situation occurs for large tandem duplications, except there is one copy in the reference and two
    (or more) in the contig. Aligners may align through the reference copy, break the alignment, and start a new
    alignment through the second copy in the contig. In this case, this function would be asked to resolve reference
    coordinates (match_coord = "subject").

    :param record_l: Pandas Series alignment record (generated by align_get_read_bed).
    :param record_r: Pandas Series alignment record (generated by align_get_read_bed).
    :param match_coord: "query" to trim contig alignments, or "subject" to match reference alignments.
    :param rev_l: Trim `record_l` from the downstream end (alignment end) if `True`, otherwise, trim from the upstream
        end (alignment start).
    :param rev_r: Trim `record_r` from the downstream end (alignment end) if `True`, otherwise, trim from the upstream
        end (alignment start).

    :return: A tuple of modified `record_l` and 'record_r`.
    """

    record_l = record_l.copy()
    record_r = record_r.copy()

    # Check arguments
    if match_coord not in {'query', 'subject'}:
        raise RuntimeError('Unknown match_coord parameter: {}: Expected "query" or "subject"'.format(match_coord))

    # Get cigar operations (list of (op_len, op_code) tuples)
    cigar_l = list(cigar_str_to_tuples(record_l))
    cigar_r = list(cigar_str_to_tuples(record_r))

    # Orient CIGAR operations so regions to be trimmed are at the head of the list
    if rev_l:
        cigar_l = cigar_l[::-1]

    if rev_r:
        cigar_r = cigar_r[::-1]

    # Get number of bases to trim. Assumes records overlap and are oriented in
    if match_coord == 'query':

        if record_l['QUERY_TIG_POS'] < record_r['QUERY_TIG_POS']:
            diff_bp = record_l['QUERY_TIG_END'] - record_r['QUERY_TIG_POS']

        else:
            #raise RuntimeError('Contigs are incorrectly ordered in query space: {} ({}:{}) vs {} ({}:{}), match_coord={}'.format(
            #    record_l['QUERY_ID'], record_l['#CHROM'], record_l['POS'],
            #    record_r['QUERY_ID'], record_r['#CHROM'], record_r['POS'],
            #    match_coord
            #))

            diff_bp = record_r['QUERY_TIG_END'] - record_l['QUERY_TIG_POS']

        if diff_bp <= 0:
            raise RuntimeError('Cannot trim to negative distance {}: {} ({}:{}) vs {} ({}:{}), match_coord={}'.format(
                diff_bp,
                record_l['QUERY_ID'], record_l['#CHROM'], record_l['POS'],
                record_r['QUERY_ID'], record_r['#CHROM'], record_r['POS'],
                match_coord
            ))
    else:
        if record_l['POS'] > record_r['POS']:
            raise RuntimeError('Contigs are incorrectly ordered in subject space: {} ({}:{}) vs {} ({}:{}), match_coord={}'.format(
                record_l['QUERY_ID'], record_l['#CHROM'], record_l['POS'],
                record_r['QUERY_ID'], record_r['#CHROM'], record_r['POS'],
                match_coord
            ))

        diff_bp = record_l['END'] - record_r['POS']

        if diff_bp <= 0:
            raise RuntimeError('Cannot trim to negative distance {}: {} ({}:{}) vs {} ({}:{}), match_coord={}'.format(
                diff_bp,
                record_l['QUERY_ID'], record_l['#CHROM'], record_l['POS'],
                record_r['QUERY_ID'], record_r['#CHROM'], record_r['POS'],
                match_coord
            ))

    # Find the number of upstream (l) bases to trim to get to 0 (or contig start)
    trace_l = trace_cigar_to_zero(cigar_l, diff_bp, record_l, match_coord == 'query')

    # Find the number of downstream (r) bases to trim to get to 0 (or contig start)
    trace_r = trace_cigar_to_zero(cigar_r, diff_bp, record_r, match_coord == 'query')


    # For each upstream alignment cut-site, find the best matching downstream alignment cut-site. Not all cut-site
    # combinations need to be tested since trimmed bases and event count is non-decreasing as it moves away from the
    # best cut-site (residual overlapping bases 0 and maximum events consumed)

    # Find optimal cut sites.
    # cut_idx_l and cut_idx_r are indices to trace_l and trace_r. These trace records point to the last CIGAR operation to
    # survive the cut, although they may be truncated. The whole record will not be removed.
    cut_idx_l, cut_idx_r = find_cut_sites(trace_l, trace_r, diff_bp)

    # Check for no cut-sites. Should not occur at this stage
    if cut_idx_l is None or cut_idx_r is None:
        raise RuntimeError('Program bug: Found no cut-sites: {} (INDEX={}) vs {} (INDEX={}), match_coord={}'.format(
            record_l['QUERY_ID'], record_l['INDEX'],
            record_r['QUERY_ID'], record_r['INDEX'],
            match_coord
        ))

    # Get cut records
    cut_l = trace_l[cut_idx_l]
    cut_r = trace_r[cut_idx_r]

    # Set mid-record cuts (Left-align cuts, mismatch first)
    residual_bp = diff_bp - (cut_l[TC_DIFF_CUM] + cut_r[TC_DIFF_CUM])
    trim_l = 0
    trim_r = 0

    if residual_bp > 0 and cut_r[TC_OP_CODE] == 'X':  # Right mismatch
        trim_r += np.min([residual_bp, cut_r[TC_OP_LEN] - 1])
        residual_bp -= trim_r

    if residual_bp > 0 and cut_l[TC_OP_CODE] == 'X':  # Left mismatch
        trim_l += np.min([residual_bp, cut_l[TC_OP_LEN] - 1])
        residual_bp -= trim_l

    if residual_bp > 0 and cut_l[TC_OP_CODE] == '=':  # Left match
        trim_l += np.min([residual_bp, cut_l[TC_OP_LEN] - 1])
        residual_bp -= trim_l

    if residual_bp > 0 and cut_r[TC_OP_CODE] == '=':  # Right match
        trim_r += np.min([residual_bp, cut_r[TC_OP_LEN] - 1])
        residual_bp -= trim_r

    # Get cut CIGAR String
    cigar_l_mod = cigar_l[cut_l[TC_INDEX]:]
    cigar_r_mod = cigar_r[cut_r[TC_INDEX]:]

    # Shorten last alignment record if set.
    cigar_l_mod[0] = (cigar_l_mod[0][0] - trim_l, cigar_l_mod[0][1])
    cigar_r_mod[0] = (cigar_r_mod[0][0] - trim_r, cigar_r_mod[0][1])

    # Modify alignment records
    record_l_mod = record_l.copy()
    record_r_mod = record_r.copy()

    cut_sub_l = cut_l[TC_SUB_BP] + trim_l
    cut_qry_l = cut_l[TC_QRY_BP] + trim_l

    cut_sub_r = cut_r[TC_SUB_BP] + trim_r
    cut_qry_r = cut_r[TC_QRY_BP] + trim_r

    if rev_l:
        record_l_mod['END'] -= cut_sub_l
        record_l_mod['QUERY_END'] -= cut_qry_l

        # Adjust positions in contig space
        if record_l_mod['REV']:
            record_l_mod['QUERY_TIG_POS'] += cut_qry_l
        else:
            record_l_mod['QUERY_TIG_END'] -= cut_qry_l

        # Track cut bases
        record_l_mod['CUT_REF_R'] += cut_sub_l
        record_l_mod['CUT_TIG_R'] += cut_qry_l

    else:
        record_l_mod['POS'] += cut_sub_l
        record_l_mod['QUERY_POS'] += cut_qry_l

        # Adjust positions in contig space
        if record_l_mod['REV']:
            record_l_mod['QUERY_TIG_END'] -= cut_qry_l
        else:
            record_l_mod['QUERY_TIG_POS'] += cut_qry_l

        # Track cut bases
        record_l_mod['CUT_REF_L'] += cut_sub_l
        record_l_mod['CUT_TIG_L'] += cut_qry_l

    if rev_r:
        record_r_mod['END'] -= cut_sub_r
        record_r_mod['QUERY_END'] -= cut_qry_r

        # Adjust positions in contig space
        if record_r_mod['REV']:
            record_r_mod['QUERY_TIG_POS'] += cut_qry_r
        else:
            record_r_mod['QUERY_TIG_END'] -= cut_qry_r

        # Track cut bases
        record_r_mod['CUT_REF_R'] += cut_sub_r
        record_r_mod['CUT_TIG_R'] += cut_qry_r

    else:
        record_r_mod['POS'] += cut_sub_r
        record_r_mod['QUERY_POS'] += cut_qry_r

        # Adjust positions in contig space
        if record_r_mod['REV']:
            record_r_mod['QUERY_TIG_END'] -= cut_qry_r
        else:
            record_r_mod['QUERY_TIG_POS'] += cut_qry_r

        # Track cut bases
        record_r_mod['CUT_REF_L'] += cut_sub_r
        record_r_mod['CUT_TIG_L'] += cut_qry_r

    # Add clipped bases to CIGAR
    if cut_l[TC_CLIPH_BP] > 0:
        cigar_l_pre = [(cut_l[TC_CLIPH_BP], 'H')]
    else:
        cigar_l_pre = []

    if cut_r[TC_CLIPH_BP] > 0:
        cigar_r_pre = [(cut_r[TC_CLIPH_BP], 'H')]
    else:
        cigar_r_pre = []

    clip_s_l = cut_l[TC_CLIPS_BP] + cut_l[TC_QRY_BP] + trim_l
    clip_s_r = cut_r[TC_CLIPS_BP] + cut_r[TC_QRY_BP] + trim_r

    if clip_s_l > 0:
        cigar_l_pre.append((clip_s_l, 'S'))

    if clip_s_r > 0:
        cigar_r_pre.append((clip_s_r, 'S'))

    # Append remaining CIGAR
    cigar_l_mod = cigar_l_pre + cigar_l_mod
    cigar_r_mod = cigar_r_pre + cigar_r_mod

    if rev_l:
        cigar_l_mod = cigar_l_mod[::-1]

    if rev_r:
        cigar_r_mod = cigar_r_mod[::-1]

    record_l_mod['CIGAR'] = ''.join([str(cigar_len) + cigar_op for cigar_len, cigar_op in cigar_l_mod])
    record_r_mod['CIGAR'] = ''.join([str(cigar_len) + cigar_op for cigar_len, cigar_op in cigar_r_mod])

    # Adjust reference and contig span after trimming
    record_l_mod['QUERY_LEN'] = record_l_mod['QUERY_END'] - record_l_mod['QUERY_POS']
    record_l_mod['SUB_LEN'] = record_l_mod['END'] - record_l_mod['POS']

    record_r_mod['QUERY_LEN'] = record_r_mod['QUERY_END'] - record_r_mod['QUERY_POS']
    record_r_mod['SUB_LEN'] = record_r_mod['END'] - record_r_mod['POS']

    # Return trimmed records
    return record_l_mod, record_r_mod


def find_cut_sites(trace_l, trace_r, diff_bp):
    """
    Find best cut-sites for left and right alignments to consume `diff_bp` bases.

    Optimize by:
    1) `diff_bp` or more bases removed.
    2) Maximize events (I, D, X)
    3) Tie-break by:
      a) Total removed bases closest to `diff_bp`.
      b) Left-align break (trace_l is preferentially trimmed when there is a tie).

    :param trace_l: List of tuples for the left alignment generated by `trace_cigar_to_zero()`. This list
        should be reversed after it was generated (traversed start to end).
    :param trace_r: List of tuples for the right alignment generated by `trace_cigar_to_zero()`.
    :param diff_bp: Target removing this many bases. Could be subject (ref) or query (tig) depending on how the
        traces were constructed.

    :return: Tuple of (cut_idx_l, cut_idx_r). cut_idx_l and cut_idx_r are the left contig and right contig cigar list index (argument to
        trace_cigar_to_zero()), index element of `trace_l` and `trace_r`) where the alignment cuts should occur.
    """

    # Right-index traversal
    tc_idx_r = 0        # Current right-index in trace record list (tc)
    #tc_idx_r_last = -1  # Last right-index. Used for early-exit.

    len_r = len(trace_r)  # End of r-indexes

    # Optimal cut-site for this pair of alignments
    cut_idx_l = None  # Record where cut occurs in left trace
    cut_idx_r = None  # Record where cut occurs in right trace
    
    max_event = 0            # Maximum number of events that may be cut
    max_diff_optimal = None  # Optimal difference in the number of bases cut over diff_bp. closest to 0 means cut-site
                             # can be placed exactly and does not force over-cutting to remove overlap)

    # Debugging: List of operations
    # diff_list_all = list()  # DBGTMP

    # Traverse l cut-sites
    for tc_idx_l in range(len(trace_l) - 1, -1, -1):

        # Optimal cutsite for this pair of alignments at a given this left index.
        cut_idx_part_l = None
        cut_idx_part_r = None

        max_event_part = 0
        max_diff_optimal_part = None

        # Note: "=" and "X" consume both subj and qry, so diff calculation is the same if subject or query is cut.
        # The following code assumes an = or X record is being processed.

        # Get min and max base differences achievable by cutting at the end or beginning of this l-record.
        min_bp_l = trace_l[tc_idx_l][TC_DIFF_CUM]
        max_bp_l = trace_l[tc_idx_l][TC_DIFF_CUM] + trace_l[tc_idx_l][TC_DIFF] - 1  # Cut all but one left base

        # Traverse r cut-sites until max-left + max-right base difference diff_bp or greater.
        while (
                tc_idx_r + 1 < len_r and
                max_bp_l + trace_r[tc_idx_r][TC_DIFF_CUM] + trace_r[tc_idx_r][TC_DIFF] - 1 < diff_bp  # Cut all but one right base
        ):
            tc_idx_r += 1

        # Traverse all cases where max-cutting the left event crosses 0 residual bases (or the single case resulting in
        # over-cutting). After this loop, the range of acceptable right indices spans tc_idx_r_start to tc_idx_r (exclusive on right side).
        tc_idx_r_start = tc_idx_r

        while (
                tc_idx_r < len_r and (
                    min_bp_l + trace_r[tc_idx_r][TC_DIFF_CUM] <= diff_bp or  # Acceptable cut site not found
                    tc_idx_r == tc_idx_r_start  # Find at least one cut-site on the right side, even if it over-cuts.
                )
        ):

            # Collect cut-site stats
            min_bp = min_bp_l + trace_r[tc_idx_r][TC_DIFF_CUM]
            max_bp = max_bp_l + trace_r[tc_idx_r][TC_DIFF_CUM] + trace_r[tc_idx_r][TC_DIFF] - 1

            diff_min = diff_bp - max_bp

            # Count number of events if the minimal cut at these sites are made.
            event_count = trace_l[tc_idx_l][TC_EVENT_CUM] + trace_r[tc_idx_r][TC_EVENT_CUM]

            if diff_min <= 0:
                # Target cut length is within the minimum and maximum bp by cutting at this site

                # Add up to bases_mismatch (number of X sites) to cut at target length (diff_bp)
                event_count += np.min([
                    # Difference on min-bases to target
                    diff_bp - diff_min,

                    # Number of events that could be added if either or both ends X and are fully cut. Since records
                    # Cannot be fully truncated, remove one event for X records.
                    (
                        trace_l[tc_idx_l][TC_EVENT] +  # Number of left events by dropping whole left record
                        trace_r[tc_idx_r][TC_EVENT] -  # Number of right events by dropping whole right record
                        (1 if trace_l[tc_idx_l][TC_EVENT] > 0 else 0) -  # Cannot cut whole left record
                        (1 if trace_r[tc_idx_r][TC_EVENT] > 0 else 0)    # Cannot cut whole right record
                    )
                ])

                # Cannot cut whole record, so remove one from event_count if count for left or right is greater than 0
                # if trace_l[tc_idx_l][TC_EVENT] > 0 and event_count > 0:
                #     event_count -= 1
                #
                # if trace_r[tc_idx_r][TC_EVENT] > 0 and event_count > 0:
                #     event_count -= 1

                diff_optimal = 0  # diff_bp is exactly achievable
            else:
                # Must over-cut to use these sites.
                diff_optimal = diff_min

            # else: Cut will remove more than the target number of bases; do not cut into these (accept full records)

            # DBGTMP: Turn on to get full list
            # diff_list_all.append(pd.Series(
            #     [
            #         tc_idx_l, tc_idx_r,
            #         trace_l[tc_idx_l][TC_INDEX],
            #         trace_r[tc_idx_r][TC_INDEX],
            #         diff_min,
            #         diff_bp - min_bp,
            #         diff_optimal,
            #         event_count
            #     ],
            #     index=['TC_IDX_L', 'TC_IDX_R',
            #            'INDEX_L', 'INDEX_R',
            #            'DIFF_MIN', 'DIFF_MAX', 'DIFF_BP',
            #            'EVENT_COUNT'
            #    ]
            # ))

            # Save max
            if (
                event_count > max_event_part or (  # Better event count, or
                    event_count == max_event_part and (  # Same event count, and
                        max_diff_optimal_part is None or diff_optimal < max_diff_optimal_part  # Optimal difference is closer to 0 (less over-cut)
                    )
                )
            ):
                cut_idx_part_l = tc_idx_l
                cut_idx_part_r = tc_idx_r
                max_event_part = event_count
                max_diff_optimal_part = diff_optimal

            tc_idx_r += 1

        # Save max
        if (
            max_event_part > max_event or (  # Better event count, or
                max_event_part == max_event and (  # Same event count, and
                    max_diff_optimal is None or max_diff_optimal_part < max_diff_optimal  # Optimal difference is closer to 0 (less over-cut)
                )
            )
        ):
            cut_idx_l = cut_idx_part_l
            cut_idx_r = cut_idx_part_r
            max_event = max_event_part
            max_diff_optimal = max_diff_optimal_part

        # Reset right index
        tc_idx_r = tc_idx_r_start

    # df_diff_all = pd.concat(diff_list_all, axis=1).T  # DBGTMP

    return cut_idx_l, cut_idx_r


def trace_cigar_to_zero(cigar_list, diff_bp, aln_record, diff_query):
    """
    Trace CIGAR operations back until diff_bp query bases are discarded from the alignment. CIGAR operations must only
    contain operators "IDSH=X" (no "M"). The list returned is only alignment match ("=" or "X" records) for the
    optimal-cut algrothm (can only cut at aligned bases).

    Returns a list of tuples for each CIGAR operation traversed:
        * TC_INDEX = 0: Index in cigar_list.
        * TC_OP_LEN = 1: CIGAR operation length.
        * TC_OP_CODE = 2: CIGAR operation code (character, e.g. "I", "=").
        * TC_DIFF_CUM = 3: Cumulative base difference up this event, but not including it.
        * TC_DIFF = 4: Base difference for this event. Will be oplen depending on the operation code.
        * TC_EVENT_CUM = 5: Cumulative event difference (number of insertions, deletions, and SNVs) up to this event,
            but not including it.
        * TC_EVENT = 6: Event differences for this event. "1" for insertions or deletions, OP_LEN for mismatches
            "X", SNV).
        * TC_SUB_BP = 7: Cumulative number of subject (ref) bases consumed up to this event, but not including it.
        * TC_QRY_BP = 8: Cumulative number of query (tig) bases consumed up to this event, but not including it.
        * TC_CLIPS_BP = 9: Cumulative number of soft-clipped bases up to AND INCLUDING this event. Alignments are not
            cut on clipped records, so cumulative and including does not affect the algorithm.
        * TC_CLIPH_BP = 10: Cumulative number of hard-clipped bases up to AND INCLUDING this event.

    :param cigar_list: List of cigar operation tuples (cigar_len, cigar_op) with cigar_op as characters (e.g. "X", "=").
    :param diff_bp: Number of query bases to trace back. Final record will traverse past this value.
    :param aln_record: Alignment record for error reporting.
    :param diff_query: Compute base differences for query (tig) sequence if `True`. If `False`, compute for subject
        (reference).

    :return: A list of tuples tracing the effects of truncating an alignment at a given CIGAR operation.
    """

    index = 0
    index_end = len(cigar_list)

    cigar_count = 0

    diff_cumulative = 0
    event_cumulative = 0

    sub_bp_sum = 0
    qry_bp_sum = 0
    clip_s_sum = 0
    clip_h_sum = 0

    trace_list = list()

    while index < index_end and (diff_cumulative <= diff_bp or cigar_list[index][1] not in {'=', 'X'} or len(trace_list) == 0):
        cigar_count += 1
        cigar_len, cigar_op = cigar_list[index]

        if cigar_op == '=':
            event_count = 0

            sub_bp = cigar_len
            qry_bp = cigar_len

        elif cigar_op == 'X':
            event_count = cigar_len

            sub_bp = cigar_len
            qry_bp = cigar_len

        elif cigar_op == 'I':
            event_count = 1

            sub_bp = 0
            qry_bp = cigar_len

        elif cigar_op == 'D':
            event_count = 1

            sub_bp = cigar_len
            qry_bp = 0

        elif cigar_op == 'S':
            event_count = 0

            sub_bp = 0
            qry_bp = 0

            clip_s_sum += cigar_len

        elif cigar_op == 'H':
            event_count = 0

            sub_bp = 0
            qry_bp = 0

            clip_h_sum += cigar_len

        else:
            raise RuntimeError((
                'Illegal operation in contig alignment while trimming alignment: {}{} '
                '(start={}:{}): CIGAR operation #{}: Expected CIGAR op in "IDSH=X"'
            ).format(cigar_len, cigar_op, aln_record['#CHROM'], aln_record['POS'], index))

        # Get number of bases affected by this event
        if diff_query:
            diff_change = qry_bp
        else:
            diff_change = sub_bp

        # Add to trace list
        if cigar_op in {'=', 'X'}:
            trace_list.append(
                (
                    index,
                    cigar_len, cigar_op,
                    diff_cumulative, diff_change,
                    event_cumulative, event_count,
                    sub_bp_sum, qry_bp_sum,
                    clip_s_sum, clip_h_sum
                )
            )

        # Increment cumulative counts
        diff_cumulative += diff_change
        event_cumulative += event_count

        sub_bp_sum += sub_bp
        qry_bp_sum += qry_bp

        index += 1

    return trace_list


def cigar_str_to_tuples(record):
    """
    Get an iterator for cigar operation tuples. Each tuple is (cigar-len, cigar-op).

    :param record: Alignment record.

    :return: Iterator of CIGAR operation tuples.
    """

    cigar = record['CIGAR']

    pos = 0
    max_pos = len(cigar)

    while pos < max_pos:

        len_pos = pos

        while cigar[len_pos] in _INT_STR_SET:
            len_pos += 1

        if len_pos == pos:
            raise RuntimeError('Missing length in CIGAR string for contig {} alignment starting at {}:{}: CIGAR index {}'.format(
                record['QUERY_ID'], record['#CHROM'], record['POS'], pos
            ))

        if cigar[len_pos] not in _CIGAR_OP_SET:
            raise RuntimeError('Unknown CIGAR operation for contig {} alignment starting at {}:{}: CIGAR operation {}'.format(
                record['QUERY_ID'], record['#CHROM'], record['POS'], cigar[pos]
            ))

        yield((int(cigar[pos:len_pos]), cigar[len_pos]))

        pos = len_pos + 1


def match_bp(record, right_end):
    """
    Get the number of matching bases at the end of an alignment. Used by variant callers to left-align SVs through
    alignment-truncating events.

    :param record: Alignment record (from alignment BED) with CIGAR string.
    :param right_end: `True` if matching alignments from the right end of `record`, or `False` to match from
        the left end.

    :return: Minimum of the number of matched bases at the end of two alignment records.
    """

    cigar = list(cigar_str_to_tuples(record))

    if right_end:
        cigar = cigar[::-1]

    # Get match base count (CIGAR op "=") on a
    match_count = 0

    for cigar_len, cigar_op in cigar:
        if cigar_op in {4, 5}:  # Skip clipped bases: S, H
            continue

        if cigar_op == 7:  # Matched bases: =
            match_count += cigar_len

        elif cigar_op == 0:
            raise RuntimeError(
                'Detected "M" opcodes in CIGAR string for record INDEX={}: Sequence match/mismatch opcodes are required ("=", "X")'.format(
                    record['INDEX'] if 'INDEX' in record.index else '<UNKNOWN>'
                )
            )
        else:
            break  # No more bases to traverse

    return match_count


def get_max_cluster(df, chrom, min_prop=0.85, min_aln_len=1e6):
    """
    For a chromosome, get the name of the cluster with the most aligned bases. When contigs are assigned to a chromosome cluster
    as part of the assembly process (Strand-seq phased assembly pipeline), this can be used to filter erroneously mapped contigs
    belonging to other chromosomes (SD driven mapping errors).

    :param df: DataFrame of alignments.
    :param chrom: Chromosome to check.
    :param min_prop: Only return if the cluster accounts for at least this proportion of aligned bases on the chromosome.
    """

    subdf = df.loc[df['#CHROM'] == chrom]

    if np.sum(subdf['SUB_LEN']) < min_aln_len or subdf.shape[0] == 1:
        return None

    cluster_count = df.loc[df['#CHROM'] == chrom].groupby('CLUSTER')['SUB_LEN'].sum().sort_values(ascending=False)

    if cluster_count.shape[0] == 0 or cluster_count.iloc[0] / np.sum(cluster_count) < min_prop:
        return None

    return cluster_count.index[0]


class AlignLift:
    """
    Create an alignment liftover object for translating between reference and contig alignments (both directions). Build
    liftover from alignment data in a DataFrame (requires #CHROM, POS, END, QUERY_ID, QUERY_POS, and QUERY_END). The
    DataFrame index must not have repeated values. The data-frame is also expected not to change, so feed this object
    a copy if it might be altered while this object is in use.
    """

    def __init__(self, df, df_fai, cache_align=10):
        """
        Create AlignLift instance.

        :param df: Alignment BED file (post-cut) as DataFrame.
        :param df_fai: Contig FAI file as Series.
        :param cache_align: Number of alignment records to cache.
        """

        self.df = df
        self.df_fai = df_fai
        self.cache_align = cache_align

        # Check df
        if len(set(df.index)) != df.shape[0]:
            raise RuntimeError('Cannot create AlignLift object with duplicate index values')

        # Build a reference-coordinate and a tig-coordinate tree
        self.ref_tree = collections.defaultdict(intervaltree.IntervalTree)
        self.tig_tree = collections.defaultdict(intervaltree.IntervalTree)

        for index, row in df.iterrows():
            self.ref_tree[row['#CHROM']][row['POS']:row['END']] = index
            self.tig_tree[row['QUERY_ID']][row['QUERY_TIG_POS']:row['QUERY_TIG_END']] = index

        # Build alignment caching structures
        self.cache_queue = collections.deque()

        self.ref_cache = dict()
        self.tig_cache = dict()

    def lift_to_sub(self, query_id, coord, gap=False):
        """
        Lift coordinates from query (tig) to subject (reference).

        Returns tuple(s) of:
            0: Query ID.
            1: Position.
            2: Is reverse complemented. `True` if the alignment was lifted through a reverse-complemented
                alignment record.
            3: Minimum position.
            4: Maximum position.
            5: Align record index.

        If the lift hits an alignment gap, then the function can try to resolve the lift to the best subject position
        based on the alignment using a strategy outlined by the `gap` parameter.

        :param query_id: Query record ID.
        :param coord: Query coordinates. May be a single value, list, or tuple.
        :param gap: Interpolate into an alignment gap between two contigs if `True`.

        :return: Single coordinate tuple or list of coordinate tuples if `coord` is a list or tuple. Returns `None`
            for coordinates that cannot be lifted.
        """

        # Determine type
        if issubclass(coord.__class__, list) or issubclass(coord.__class__, tuple):
            ret_list = True
        else:
            ret_list = False
            coord = (coord,)

        # Do lift
        lift_coord_list = list()

        for pos in coord:
            pos_org = pos  # Pre-reverse position

            # Find matching records
            match_set = self.tig_tree[query_id][pos:(pos + 1)]

            if len(match_set) == 1:
                index = list(match_set)[0].data

            elif len(match_set) == 0 and gap:
                lift_coord_list.append(self._get_subject_gap(query_id, pos))
                continue

            else:
                lift_coord_list.append(None)
                continue

            # Get lift tree
            if index not in self.tig_cache.keys():
                self._add_align(index)

            lift_tree = self.tig_cache[index]

            # Get row
            row = self.df.loc[index]

            # Reverse coordinates of pos if the alignment is reverse-complemented. Translates from QUERY_TIG_POS-space
            # to QUERY_POS-space.
            if row['REV']:
                pos = self.df_fai[query_id] - pos

            # Get match record
            match_set = lift_tree[pos]

            if len(match_set) == 1:
                match_interval = list(match_set)[0]

            if len(match_set) == 0:
                # Allow queries to match if they end exactly at the alignment end
                match_set = lift_tree[pos - 1]

                match_interval = list(match_set)[0] if len(match_set) == 1 else None

                if not match_interval or match_interval.end != pos:
                    lift_coord_list.append(None)

                    raise RuntimeError(
                        (
                            'Found no matches in a lift-tree for a record within a '
                            'global to-subject tree: {}:{} (index={}, gap={})'
                        ).format(query_id, pos_org, index, gap)
                    )

            elif len(match_set) > 1:
                lift_coord_list.append(None)

                raise RuntimeError(
                    (
                        'Found multiple matches in a lift-tree for a record within a '
                        'global to-subject tree: {}:{} (index={}, gap={})'
                    ).format(query_id, pos_org, index,  gap)
                )

            # Interpolate coordinates
            if match_interval.data[1] - match_interval.data[0] > 1:
                lift_pos = match_interval.data[0] + (pos - match_interval.begin)

                lift_coord_list.append((
                    row['#CHROM'],
                    lift_pos,
                    row['REV'],
                    lift_pos,
                    lift_pos,
                    (row['INDEX'],)
                ))

            else:  # Lift from missing bases on the target (insertion or deletion)
                lift_coord_list.append((
                    row['#CHROM'],
                    match_interval.data[1],
                    row['REV'],
                    match_interval.data[1],
                    match_interval.data[1],
                    (row['INDEX'],)
                ))

        # Return coordinates
        if ret_list:
            return lift_coord_list
        else:
            return lift_coord_list[0]

    def lift_to_qry(self, subject_id, coord):
        """
        Lift coordinates from subject (reference) to query (tig).

        Returns tuple(s) of:
            0: Query ID.
            1: Position.
            2: Is reverse complemented. `True` if the alignment was lifted through a reverse-complemented
                alignment record.
            3: Minimum position.
            4: Maximum position.
            5: Align record index.

        :param subject_id: Subject ID.
        :param coord: Subject coordinates. May be a single value, list, or tuple.

        :return: Single coordinate tuple or list of coordinate tuples if `coord` is a list or tuple. Returns `None`
            for coordinates that cannot be lifted.
        """

        # Determine type
        if issubclass(coord.__class__, list) or issubclass(coord.__class__, tuple):
            ret_list = True
        else:
            ret_list = False
            coord = (coord,)

        # Do lift
        lift_coord_list = list()

        for pos in coord:
            match_set = self.ref_tree[subject_id][pos:(pos + 1)]

            # Check coordinates
            if len(match_set) == 0:
                lift_coord_list.append(None)
                continue

                # raise ValueError('Subject region {}:{} has no to-query lift records'.format(subject_id, pos))

            if len(match_set) > 1:
                lift_coord_list.append(None)
                continue

                # raise ValueError(
                #     'Subject region {}:{} has {} to-query lift records'.format(subject_id, pos, len(match_set))
                # )

            # Get lift tree
            index = list(match_set)[0].data

            if index not in self.ref_cache.keys():
                self._add_align(index)

            lift_tree = self.ref_cache[index]

            # Save row
            row = self.df.loc[index]

            # Get match record
            match_set = lift_tree[pos:(pos + 1)]

            if len(match_set) != 1:
                raise RuntimeError(
                    (
                        'Program bug: Found no matches in a lift-tree for a record withing a '
                        'global to-query tree: {}:{} (index={})'
                    ).format(subject_id, pos, index)
                )

            match_interval = list(match_set)[0]

            # Interpolate coordinates
            if match_interval.data[1] - match_interval.data[0] > 1:
                qry_pos = match_interval.data[0] + (pos - match_interval.begin)

            else:  # Lift from missing bases on the target (insertion or deletion)
                qry_pos = match_interval.data[1]

            if row['REV']:
                qry_pos = self.df_fai[row['QUERY_ID']] - qry_pos

            lift_coord_list.append((
                row['QUERY_ID'],
                qry_pos,
                row['REV'],
                qry_pos,
                qry_pos,
                (row['INDEX'],)
            ))

        # Return coordinates
        if ret_list:
            return lift_coord_list
        else:
            return lift_coord_list[0]

    def lift_region_to_sub(self, region, gap=False):
        """
        Lift region to subject.

        :param region: Query region.
        :param gap: Interpolate within gap if `True`.

        :return: Subject region or `None` if it could not be lifted.
        """

        # Lift
        sub_pos, sub_end = self.lift_to_sub(region.chrom, (region.pos, region.end), gap)

        # Check lift: Must lift both ends to the same subject ID
        if sub_pos is None or sub_end is None:
            return None

        if sub_pos[0] != sub_end[0] or (sub_pos[2] is not None and sub_end[2] is not None and sub_pos[2] != sub_end[2]):
            return None

        # Return
        return pavlib.seq.Region(
            sub_pos[0], sub_pos[1], sub_end[1],
            is_rev=False,
            pos_min=sub_pos[3], pos_max=sub_pos[4],
            end_min=sub_end[3], end_max=sub_end[4],
            pos_aln_index=(sub_pos[5],),
            end_aln_index=(sub_end[5],)
        )

    def lift_region_to_qry(self, region):
        """
        Lift region to query.

        :param region: Subject region.

        :return: Query region or `None` if it could not be lifted.
        """

        # Lift
        query_pos, query_end = self.lift_to_qry(region.chrom, (region.pos, region.end))

        # Check lift: Must lift both ends to the same query ID
        if query_pos is None or query_end is None:
            return None

        if query_pos[0] != query_end[0] or query_pos[2] != query_end[2]:
            return None

        # Return
        return pavlib.seq.Region(
            query_pos[0], query_pos[1], query_end[1],
            is_rev=query_pos[2],
            pos_min=query_pos[3], pos_max=query_pos[4],
            end_min=query_end[3], end_max=query_end[4],
            pos_aln_index=(query_pos[5],),
            end_aln_index=(query_end[5],)
        )

    def _get_subject_gap(self, query_id, pos):
        """
        Interpolate lift coordinates to an alignment gap.

        :param pos: Position on the contig.

        :return: A tuple of (subject, pos, rev, min, max).
        """

        # Check arguments
        if pos is None:
            return None

        # Get alignment records for this contig
        subdf = self.df.loc[self.df['QUERY_ID'] == query_id]

        # Must be flanked by two contigs on either side
        if not np.any(subdf['QUERY_TIG_END'] < pos):
            return None

        if not np.any(subdf['QUERY_TIG_POS'] > pos):
            return None

        # Get left and right rows for the alignment record flanking this position
        row_l = subdf.loc[
            subdf.loc[subdf['QUERY_TIG_END'] < pos, 'QUERY_TIG_END'].sort_values().index[-1]
        ]

        row_r = subdf.loc[
            subdf.loc[subdf['QUERY_TIG_POS'] > pos, 'QUERY_TIG_POS'].sort_values().index[0]
        ]

        # Rows must be mapped to the same subject
        if row_l['#CHROM'] != row_r['#CHROM']:
            return None

        return (
            (
                row_l['#CHROM'],
                np.int64((row_l['QUERY_TIG_END'] + row_r['QUERY_TIG_POS']) / 2),
                row_l['REV'] if row_l['REV'] == row_r['REV'] else None,
                row_l['QUERY_TIG_END'],
                row_r['QUERY_TIG_POS'],
                (row_l['INDEX'], row_r['INDEX'])
            )
        )

    def _add_align(self, index):
        """
        Add an alignment from DataFrame index `index`.

        :param index: DataFrame index.
        """

        # No alignment to add if it's already cached.
        if index in self.ref_cache.keys():

            # Append to end (last used) and skip adding alignment
            while index in self.cache_queue:
                self.cache_queue.remove(index)

            self.cache_queue.appendleft(index)

            return

        # Make space for this alignment
        self._check_and_clear()

        # Get row
        row = self.df.loc[index]

        # Build lift trees
        sub_bp = row['POS']
        qry_bp = 0

        itree_ref = intervaltree.IntervalTree()
        itree_tig = intervaltree.IntervalTree()

        # Get CIGAR and check query start position
        cigar_op_list = list(cigar_str_to_tuples(self.df.loc[index]))

        clipped_bp = 0

        cigar_index = 0

        while cigar_index < len(cigar_op_list) and cigar_op_list[cigar_index][1] in {'S', 'H'}:
            clipped_bp += cigar_op_list[cigar_index][0]
            cigar_index += 1

        if row['QUERY_POS'] != clipped_bp:
            raise RuntimeError(
                'Number of clipped bases ({}) does not match query start position {}: Alignment {}:{} ({}:{})'.format(
                    clipped_bp, row['QUERY_POS'], row['#CHROM'], row['POS'], row['QUERY_ID'], row['QUERY_POS']
                )
            )

        # Build trees
        for cigar_len, cigar_op in cigar_op_list:

            if cigar_op in {'=', 'X', 'M'}:

                itree_ref[sub_bp:(sub_bp + cigar_len)] = (qry_bp, qry_bp + cigar_len)
                itree_tig[qry_bp:(qry_bp + cigar_len)] = (sub_bp, sub_bp + cigar_len)

                sub_bp += cigar_len
                qry_bp += cigar_len

            elif cigar_op == 'I':

                itree_tig[qry_bp:(qry_bp + cigar_len)] = (sub_bp, sub_bp + 1)

                qry_bp += cigar_len

            elif cigar_op == 'D':

                itree_ref[sub_bp:(sub_bp + cigar_len)] = (qry_bp, qry_bp + 1)

                sub_bp += cigar_len

            elif cigar_op in {'S', 'H'}:

                qry_bp += cigar_len

                clipped_bp += cigar_len

            else:
                row = self.df.loc[index]

                raise RuntimeError('Unhandled CIGAR operation: {}: Alignment {}:{} ({}:{})'.format(
                    cigar_op, row['#CHROM'], row['POS'], row['QUERY_ID'], row['QUERY_POS']
                ))

        # Cache trees
        self.ref_cache[index] = itree_ref
        self.tig_cache[index] = itree_tig

        # Add index to end of queue
        self.cache_queue.appendleft(index)

    def _check_and_clear(self):
        """
        Check alignment cache and clear if necessary to make space.
        """

        while len(self.cache_queue) >= self.cache_align:
            index = self.cache_queue.pop()

            del(self.ref_cache[index])
            del(self.tig_cache[index])


def check_record(row, df_tig_fai):
    """
    Check alignment DatFrame record for sanity. Throws exceptions if there are problems. Returns nothing if everything
    passes.

    :param record: Alignment table record (Pandas Series).
    :param df_tig_fai: Panadas Series with contig names as keys and contig lengths as values.
    """

    try:
        ref_bp, tig_bp, clip_h_l, clip_s_l, clip_h_r, clip_s_r = count_cigar(row)

    except Exception as ex:
        raise RuntimeError('CIGAR parsing error: {} (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(ex, row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']))

    tig_len = df_tig_fai[row['QUERY_ID']]

    # POS and END agree with length
    if row['POS'] + ref_bp != row['END']:

        raise RuntimeError(
            'END mismatch: POS + len != END ({} != {}) (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
                row['POS'] + ref_bp, row['END'], row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
            )
        )

    # Query POS and END  agree with length
    if row['QUERY_POS'] + tig_bp != row['QUERY_END']:
        raise RuntimeError(
            'QUERY_END mismatch: {} != {} (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
                row['QUERY_POS'] + tig_bp, row['QUERY_END'], row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
            )
        )

    # Query contig (non-rev-compl if alignment was reversed) POS and END agree with length
    if row['QUERY_TIG_POS'] + tig_bp != row['QUERY_TIG_END']:
        raise RuntimeError(
            'QUERY_TIG_END mismatch: QUERY_TIG_POS + tig_bp != QUERY_TIG_END ({} != {}) (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
                row['QUERY_TIG_POS'] + tig_bp, row['QUERY_TIG_END'], row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
            )
        )

    # Check QUERY_ and QUERY_TIG_ concordance against contig length
    if row['REV']:
        if row['QUERY_TIG_POS'] != tig_len - row['QUERY_END']:
            raise RuntimeError(
                'Rev and QUERY_END does not translate to QUERY_TIG_POS: QUERY_TIG_POS != tig_len - QUERY_END ({} != {}) (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
                    row['QUERY_TIG_POS'], tig_len - row['QUERY_END'], row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
                )
            )

        if row['QUERY_TIG_END'] != tig_len - row['QUERY_POS']:
            raise RuntimeError(
                'Rev and QUERY_POS does not translate to QUERY_TIG_END: QUERY_TIG_END != tig_len - QUERY_POS ({} != {}) (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
                    row['QUERY_TIG_END'], tig_len - row['QUERY_POS'], row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
                )
            )
    else:
        if row['QUERY_TIG_POS'] != row['QUERY_POS']:
            raise RuntimeError(
                'Fwd and QUERY_POS does not match QUERY_TIG_POS: QUERY_TIG_POS != QUERY_POS ({} != {}) (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
                    row['QUERY_TIG_POS'], row['QUERY_POS'], row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
                )
            )

        if row['QUERY_TIG_END'] != row['QUERY_END']:
            raise RuntimeError(
                'Fwd and QUERY_END does not match QUERY_TIG_END: QUERY_TIG_END != QUERY_END ({} != {}) (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
                    row['QUERY_TIG_END'], row['QUERY_END'], row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
                )
            )

    # Query and reference positions are in the right order
    if row['QUERY_POS'] >= row['QUERY_END']:
        raise RuntimeError('QUERY_POS >= QUERY_END ({} >= {}) (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
            row['QUERY_POS'], row['QUERY_END'], row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
        ))

    if row['QUERY_TIG_POS'] >= row['QUERY_TIG_END']:
        raise RuntimeError('QUERY_TIG_POS >= QUERY_TIG_END ({} >= {}) (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
            row['QUERY_TIG_POS'], row['QUERY_TIG_END'], row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
        ))

    if row['POS'] >= row['END']:
        raise RuntimeError('POS >= END ({} >= {}) (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
            row['POS'], row['END'], row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
        ))

    # Contig ends are not longer than contig length
    if row['QUERY_TIG_END'] > tig_len:
        raise RuntimeError('QUERY_TIG_END > tig_len ({} > {}) (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
            row['QUERY_TIG_END'], tig_len, row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
        ))

    if row['QUERY_END'] > tig_len:
        raise RuntimeError('QUERY_END > tig_len ({} > {}) (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
            row['QUERY_END'], tig_len, row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
        ))

    # Clipping agrees with contig starts
    if row['QUERY_POS'] != clip_h_l + clip_s_l:
        raise RuntimeError(
            'QUERY_POS != clip_h_l + clip_s_l ({} != {}) (INDEX={}, QUERY={}, LOC={}:{}-{})'.format(
                row['QUERY_POS'], clip_h_l + clip_s_l, row['INDEX'], row['QUERY_ID'], row['#CHROM'], row['POS'], row['END']
            )
        )


def check_record_err_string(df, df_tig_fai):
    """
    Runs check_record on each row of `df`, captures exceptions, and returns a Series of error message strings instead
    of failing on the first error. The Series can be added as a column to `df`. For each record where there was no
    error, the field for that record in the returned series is NA (`np.nan`). This function may not be used by the
    pipeline, but is here for troubleshooting alignments.

    :param df: Dataframe of alignment records.
    :param df_tig_fai: Panadas Series with contig names as keys and contig lengths as values.

    :return: A Series of error messages (or NA) for each record in `df`.
    """

    def _check_record_err_string_check_row(row):
        try:
            check_record(row, df_tig_fai)
            return None
        except Exception as ex:
            return str(ex)

    return df.apply(_check_record_err_string_check_row, axis=1)


def count_cigar(row):
    """
    Count bases affected by CIGAR operations in an alignment record (row is a Pandas Series from an ailgnment BED).

    Returns a tuple of:
    * ref_bp: Reference (subject) bases traversed by CIGAR operations.
    * tig_bp: Contig (query) bases traversed by CIGAR operations. Does not include clipped bases.
    * clip_h_l: Hard-clipped bases on the left (upstream) side.
    * clip_s_l: Soft-clipped bases on the left (upstream) side.
    * clip_h_r: Hard-clipped bases on the right (downstream) side.
    * clip_s_r: Soft-clipped bases on the right (downstream) side.

    :param row: Row with CIGAR records as a CIGAR string.

    :return: A tuple of (ref_bp, tig_bp, clip_h_l, clip_s_l, clip_h_r, clip_s_r).
    """

    ref_bp = 0
    tig_bp = 0

    clip_s_l = 0
    clip_h_l = 0

    clip_s_r = 0
    clip_h_r = 0

    cigar_list = list(pavlib.align.cigar_str_to_tuples(row))

    cigar_n = len(cigar_list)

    index = 0

    while index < cigar_n and cigar_list[index][1] in {'S', 'H'}:
        cigar_len, cigar_op = cigar_list[index]

        if cigar_op == 'S':
            if clip_s_l > 0:
                raise RuntimeError('Duplicate S records (left) at index {}'.format(index))
            clip_s_l = cigar_len

        if cigar_op == 'H':
            if clip_h_l > 0:
                raise RuntimeError('Duplicate H records (left) at index {}'.format(index))

            if clip_s_l > 0:
                raise RuntimeError('S record before H (left) at index {}'.format(index))

            clip_h_l = cigar_len

        index += 1

    while index < cigar_n:
        cigar_len, cigar_op = cigar_list[index]

        # print('{}: {} {}'.format(index, cigar_len, cigar_op))

        if cigar_op in {'=', 'X'}:

            if clip_s_r > 0 or clip_h_r > 0:
                raise RuntimeError(
                    'Found clipped bases before last non-clipped CIGAR operation at operation {} ({}{})'.format(
                        index, cigar_len, cigar_op
                    )
                )

            ref_bp += cigar_len
            tig_bp += cigar_len

        elif cigar_op == 'I':

            if clip_s_r > 0 or clip_h_r > 0:
                raise RuntimeError(
                    'Found clipped bases before last non-clipped CIGAR operation at operation {} ({}{})'.format(
                        index, cigar_len, cigar_op
                    )
                )

            tig_bp += cigar_len

        elif cigar_op == 'D':

            if clip_s_r > 0 or clip_h_r > 0:
                raise RuntimeError(
                    'Found clipped bases before last non-clipped CIGAR operation at operation {} ({}{})'.format(
                        index, cigar_len, cigar_op
                    )
                )

            ref_bp += cigar_len

        elif cigar_op == 'S':

            if clip_s_r > 0:
                raise RuntimeError('Duplicate S records (right) at operation {}'.format(index))

            if clip_h_r > 0:
                raise RuntimeError('H record before S record (right) at operation {}'.format(index))

            clip_s_r = cigar_len

        elif cigar_op == 'H':

            if clip_h_r > 0:
                raise RuntimeError('Duplicate H records (right) at operation {}'.format(index))

            clip_h_r = cigar_len

        else:
            raise RuntimeError('Bad CIGAR op: ' + cigar_op)

        index += 1

    return ref_bp, tig_bp, clip_h_l, clip_s_l, clip_h_r, clip_s_r


def get_align_bed(align_file, df_tig_fai, hap, chrom_cluster=False, min_mapq=0):
    """
    Read alignment file as a BED file that PAV can process. Drops any records marked as unaligned by the SAM flag.

    :param align_file: SAM, CRAM, BAM, anything `pysam.AlignmentFile` can read.
    :param df_tig_fai: Pandas Series with contig names as keys and contig lengths as values. Index should be cast as
        type str if contig names are numeric.
    :param hap: Haplotype assinment for this alignment file (h1 or h2).
    :param chrom_cluster: If contigs were grouped by chromosome before aligning (PGAS with Strand-seq clustering does
        this), then assign `CLUSTER_MATCH` to alignment records if the clusters can be assigned to chromosomes and
        this record's cluster matches the aligned chromosome. PAV assumes cluster names are everything before the
        first underscore ("_") in the contig name. If `False`, `CLUSTER_MATCH` is set to `np.nan`. If a cluster cannot
        be assigned to a chromosome, `CLUSTER_MATCH` is `np.nan` for those records.
    :param min_mapq: Minimum MAPQ. If 0, then all alignments are accepted as long as the unmapped flag is not set.

    :return: BED file of alignment records.
    """

    # Get records from SAM
    record_list = list()

    align_index = 0

    with pysam.AlignmentFile(align_file, 'rb') as in_file:
        for record in in_file:

            # Increment align_index
            align_index += 1

            # Skipped unmapped reads
            if record.is_unmapped or record.mapping_quality < min_mapq or len(record.cigar) == 0:
                continue

            # Get length for computing real tig positions for rev-complemented records
            tig_len = df_tig_fai[record.query_name]

            # Read tags
            tags = dict(record.get_tags())

            # Determine left hard-clipped bases.
            # pysam query alignment functions are relative to the sequence in the alignment record, not the original
            # sequence. The left-most hard-clipped bases must be added to the query positions to translate to the
            # correct contig coordinates (https://github.com/pysam-developers/pysam/issues/1017).
            cigar_tuples = record.cigartuples

            clip_h_l = 0
            cigar_index = 0

            while cigar_tuples[cigar_index][0] == 5 and cigar_index < len(record.cigar):
                clip_h_l += cigar_tuples[cigar_index][1]
                cigar_index += 1

            # Disallow alignment match (M) in CIGAR (requires =X for base match/mismatch)
            if 'M' in record.cigarstring:
                raise RuntimeError((
                    'Found alignment match CIGAR operation (M) for record {} (Start = {}:{}): '
                    'Alignment requires CIGAR base-level match/mismatch (=X)'
                ).format(record.query_name, record.reference_name, record.reference_start))

            # Save record
            record_list.append(pd.Series(
                [
                    record.reference_name,
                    record.reference_start,
                    record.reference_end,

                    align_index,

                    record.query_name,
                    record.query_alignment_start + clip_h_l,
                    record.query_alignment_end + clip_h_l,

                    tig_len - (record.query_alignment_end + clip_h_l) if record.is_reverse else record.query_alignment_start + clip_h_l,
                    tig_len - (record.query_alignment_start + clip_h_l) if record.is_reverse else record.query_alignment_end + clip_h_l,

                    tags['RG'] if 'RG' in tags else 'NA',
                    tags['AO'] if 'AO' in tags else 'NA',

                    record.mapping_quality,

                    record.is_reverse,
                    '0x{:04x}'.format(record.flag),

                    hap,
                    record.cigarstring
                ],
                index=[
                    '#CHROM', 'POS', 'END',
                    'INDEX',
                    'QUERY_ID', 'QUERY_POS', 'QUERY_END',
                    'QUERY_TIG_POS', 'QUERY_TIG_END',
                    'RG', 'AO',
                    'MAPQ',
                    'REV', 'FLAGS', 'HAP',
                    'CIGAR'
                ]
            ))

    # Merge records
    if len(record_list) > 0:
        df = pd.concat(record_list, axis=1).T
    else:
        df = pd.DataFrame(
            [],
            columns=[
                '#CHROM', 'POS', 'END',
                'INDEX',
                'QUERY_ID', 'QUERY_POS', 'QUERY_END',
                'QUERY_TIG_POS', 'QUERY_TIG_END',
                'RG', 'AO',
                'MAPQ',
                'REV', 'FLAGS', 'HAP',
                'CIGAR'
            ]
        )

    df.sort_values(['#CHROM', 'POS', 'END', 'QUERY_ID'], ascending=[True, True, False, True], inplace=True)

    # Check sanity
    df.apply(pavlib.align.check_record, df_tig_fai=df_tig_fai, axis=1)

    # Find max cluster match for each chromosome
    if chrom_cluster and df.shape[0] > 0:
        df['SUB_LEN'] = df['END'] - df['POS']

        df['CLUSTER'] = df['QUERY_ID'].apply(lambda val: val.split('_')[0])
        max_cluster = {chrom: pavlib.align.get_max_cluster(df, chrom) for chrom in set(df['#CHROM'])}

        df['CLUSTER_MATCH'] = df.apply(lambda row: row['CLUSTER'] == max_cluster[row['#CHROM']], axis=1)
        df['CLUSTER_MATCH'] = df.apply(lambda row: row['CLUSTER_MATCH'] if max_cluster[row['#CHROM']] is not None else np.nan, axis=1)

        del(df['SUB_LEN'])
    else:
        df['CLUSTER_MATCH'] = np.nan

    # Return BED
    return df