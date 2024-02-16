# Alignment trimming functions

import numpy as np

import pavlib
import svpoplib

from .align import *


def trim_alignments(df, min_trim_tig_len, tig_fai, match_tig=False, mode='both'):
    """
    Do alignment trimming from prepared alignment BED file. This BED contains information about reference and
    contig coordinates mapped, CIGAR string and flags.

    :param df: Alignment dataframe.
    :param min_trim_tig_len: Minimum alignment record length. Alignment records smaller will be discarded.
    :param tig_fai: FAI file from the contig FASTA that was mapped. Used for an alignment sanity check after trimming.
    :param match_tig: When trimming in reference space (removing redundantly mapped reference bases), only trim
        alignment records if the contig (query ID) matches. This allows multiple contigs to map to the same location,
        which may produce a redundant set of variant calls from mis-aligned contigs (e.g. incorrectly mapped paralogs).
        This mode is useful for generating a maximal callset where multiple contigs may cover the same region (e.g.
        alleles in squashed assembly or mixed haplotypes (e.g. subclonal events). Additional QC should be applied to
        callsets using this option.
    :param mode: Trim alignments to remove redundantly mapped tig bases (mode = "tig"), reference bases with more than
        one alignment (mode = "ref"), or both (mode = "both"). If None, assume "both".

    :return: Trimmed alignments as an alignment DataFrame. Same format as `df` with columns added describing the
        number of reference and contig bases that were trimmed. Dropped records (mapped inside another or too shart) are
        removed.
    """

    # Check mode
    if mode is None:
        mode = 'both'

    mode = mode.lower()

    if mode == 'tig':
        do_trim_tig = True
        do_trim_ref = False
    elif mode == 'ref':
        do_trim_tig = False
        do_trim_ref = True
    elif mode == 'both':
        do_trim_tig = True
        do_trim_ref = True
    else:
        raise RuntimeError(f'Unrecognized trimming mode "{mode}": Expected "tig", "ref", or "both"')

    # Remove short alignments
    for index in df.index:
        if df.loc[index, 'QRY_END'] - df.loc[index, 'QRY_POS'] < min_trim_tig_len:
            df.loc[index, 'INDEX'] = -1

    # Discard fully trimmed records (shouldn't occur at this stage)
    df = df.loc[df['INDEX'] >= 0].copy()


    ###                                          ###
    ### Trim overlapping contigs in contig space ###
    ###                                          ###

    if do_trim_tig:

        # Sort by alignment lengths in contig space
        df.sort_values(['QRY_ID', 'QRY_LEN'], ascending=(True, False), inplace=True)

        df.reset_index(inplace=True, drop=True)

        # Do trim in contig space #
        iter_index_l = 0
        index_max = df.shape[0]

        while iter_index_l < index_max:

            iter_index_r = iter_index_l + 1

            while iter_index_r < index_max and df.loc[iter_index_l, 'QRY_ID'] == df.loc[iter_index_r, 'QRY_ID']:

                # Get index in order of contig placement
                if df.loc[iter_index_l, 'QRY_POS'] <= df.loc[iter_index_r, 'QRY_POS']:
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
                if df.loc[index_r, 'QRY_POS'] >= df.loc[index_l, 'QRY_END']:
                    iter_index_r += 1
                    continue

                # Check for record fully contained within another
                if df.loc[index_r, 'QRY_END'] <= df.loc[index_l, 'QRY_END']:
                    df.loc[index_r, 'INDEX'] = -1
                    iter_index_r += 1
                    continue

                # Determine trim orientation (right side of index_l is to trimmed, so must be reversed so
                # trimmed CIGAR records are at the beginning; left side of index_r is to be trimmed, which is
                # already at the start of the CIGAR string).
                rev_l = not df.loc[index_l, 'REV']  # Trim right end of index_l
                rev_r = df.loc[index_r, 'REV']      # Trim left end of index_r

                # Detect overlap in reference
                if rev_l == rev_r or df.loc[index_l, '#CHROM'] != df.loc[index_r, '#CHROM']:
                    # Contigs were aligned in different reference chromosomes or orientations, no overlap
                    ref_overlap = False

                else:
                    if df.loc[index_l, 'POS'] < df.loc[index_r, 'POS']:
                        ref_overlap = df.loc[index_r, 'POS'] < df.loc[index_l, 'END']

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
                    rm_l_a = record_l_a['QRY_END'] - record_l_a['QRY_POS'] < min_trim_tig_len
                    rm_l_b = record_l_b['QRY_END'] - record_l_b['QRY_POS'] < min_trim_tig_len

                    rm_r_a = record_r_a['QRY_END'] - record_r_a['QRY_POS'] < min_trim_tig_len
                    rm_r_b = record_r_b['QRY_END'] - record_r_b['QRY_POS'] < min_trim_tig_len

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
                if record_l['QRY_END'] - record_l['QRY_POS'] >= min_trim_tig_len:
                    df.loc[index_l] = record_l
                else:
                    df.loc[index_l, 'INDEX'] = -1

                if (record_r['QRY_END'] - record_r['QRY_POS']) >= min_trim_tig_len:
                    df.loc[index_r] = record_r
                else:
                    df.loc[index_r, 'INDEX'] = -1

                # Next r record
                iter_index_r += 1

            # Next l record
            iter_index_l += 1

        # Discard fully trimmed records
        df = df.loc[df['INDEX'] >= 0].copy()

    # end: if do_trim_tig

    ###                                             ###
    ### Trim overlapping contigs in reference space ###
    ###                                             ###

    if do_trim_ref:

        # Sort by contig alignment length in reference space
        df = df.loc[
            pd.concat(
                [df['#CHROM'], df['END'] - df['POS']], axis=1
            ).sort_values(
                ['#CHROM', 0],
                ascending=(True, False)
            ).index
        ].reset_index(drop=True)

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

                # Skip if match_tig and query names differ
                if match_tig and df.loc[iter_index_l, 'QRY_ID'] != df.loc[iter_index_r, 'QRY_ID']:
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

                    # Check for record fully contained within another
                    if df.loc[index_r, 'END'] <= df.loc[index_l, 'END']:
                        # print('\t* Fully contained')

                        df.loc[index_r, 'INDEX'] = -1

                    else:

                        record_l, record_r = pavlib.align.trim_alignment_record(df.loc[index_l], df.loc[index_r], 'subject')

                        if record_l is not None and record_r is not None:

                            # Modify if new aligned size is at least min_trim_tig_len, remove if shorter
                            if record_l['QRY_END'] - record_l['QRY_POS'] >= min_trim_tig_len:
                                df.loc[index_l] = record_l
                            else:
                                df.loc[index_l, 'INDEX'] = -1

                            if (record_r['QRY_END'] - record_r['QRY_POS']) >= min_trim_tig_len:
                                df.loc[index_r] = record_r
                            else:
                                df.loc[index_r, 'INDEX'] = -1

                # Next r record
                iter_index_r += 1

            # Next l record
            iter_index_l += 1

        # Discard fully trimmed records
        df = df.loc[df['INDEX'] >= 0].copy()

    # end: if do_trim_ref


    ###                      ###
    ### Post trim formatting ###
    ###                      ###

    # Clean and re-sort
    df = df.loc[df['INDEX'] >= 0].copy()

    df = df.loc[(df['END'] - df['POS']) > 0]  # Should never occur, but don't allow 0-length records
    df = df.loc[(df['QRY_END'] - df['QRY_POS']) > 0]

    df.sort_values(['#CHROM', 'POS', 'END', 'QRY_ID'], ascending=[True, True, False, True], inplace=True)

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
    the overlap is resolved using a simple greedy algorithm that maximizes the number of variants removed from the
    alignment during trimming (each variant is an insertion, deletion, or SNVs; currently, no bonus is given to removing
    larger insertions or deletions vs smaller ones).

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

        if record_l['QRY_POS'] < record_r['QRY_POS']:
            diff_bp = record_l['QRY_END'] - record_r['QRY_POS']

        else:
            #raise RuntimeError('Contigs are incorrectly ordered in query space: {} ({}:{}) vs {} ({}:{}), match_coord={}'.format(
            #    record_l['QRY_ID'], record_l['#CHROM'], record_l['POS'],
            #    record_r['QRY_ID'], record_r['#CHROM'], record_r['POS'],
            #    match_coord
            #))

            diff_bp = record_r['QRY_END'] - record_l['QRY_POS']

        if diff_bp <= 0:
            raise RuntimeError('Cannot trim to negative distance {}: {} ({}:{}) vs {} ({}:{}), match_coord={}'.format(
                diff_bp,
                record_l['QRY_ID'], record_l['#CHROM'], record_l['POS'],
                record_r['QRY_ID'], record_r['#CHROM'], record_r['POS'],
                match_coord
            ))
    else:
        if record_l['POS'] > record_r['POS']:
            raise RuntimeError('Contigs are incorrectly ordered in subject space: {} ({}:{}) vs {} ({}:{}), match_coord={}'.format(
                record_l['QRY_ID'], record_l['#CHROM'], record_l['POS'],
                record_r['QRY_ID'], record_r['#CHROM'], record_r['POS'],
                match_coord
            ))

        diff_bp = record_l['END'] - record_r['POS']

        if diff_bp <= 0:
            raise RuntimeError('Cannot trim to negative distance {}: {} ({}:{}) vs {} ({}:{}), match_coord={}'.format(
                diff_bp,
                record_l['QRY_ID'], record_l['#CHROM'], record_l['POS'],
                record_r['QRY_ID'], record_r['#CHROM'], record_r['POS'],
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
            record_l['QRY_ID'], record_l['INDEX'],
            record_r['QRY_ID'], record_r['INDEX'],
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

        # Adjust positions in contig space
        if record_l_mod['REV']:
            record_l_mod['QRY_POS'] += cut_qry_l
        else:
            record_l_mod['QRY_END'] -= cut_qry_l

        # Track cut bases
        record_l_mod['TRIM_REF_R'] += cut_sub_l
        record_l_mod['TRIM_QRY_R'] += cut_qry_l

    else:
        record_l_mod['POS'] += cut_sub_l

        # Adjust positions in contig space
        if record_l_mod['REV']:
            record_l_mod['QRY_END'] -= cut_qry_l
        else:
            record_l_mod['QRY_POS'] += cut_qry_l

        # Track cut bases
        record_l_mod['TRIM_REF_L'] += cut_sub_l
        record_l_mod['TRIM_QRY_L'] += cut_qry_l

    if rev_r:
        record_r_mod['END'] -= cut_sub_r

        # Adjust positions in contig space
        if record_r_mod['REV']:
            record_r_mod['QRY_POS'] += cut_qry_r
        else:
            record_r_mod['QRY_END'] -= cut_qry_r

        # Track cut bases
        record_r_mod['TRIM_REF_R'] += cut_sub_r
        record_r_mod['TRIM_QRY_R'] += cut_qry_r

    else:
        record_r_mod['POS'] += cut_sub_r

        # Adjust positions in contig space
        if record_r_mod['REV']:
            record_r_mod['QRY_END'] -= cut_qry_r
        else:
            record_r_mod['QRY_POS'] += cut_qry_r

        # Track cut bases
        record_r_mod['TRIM_REF_L'] += cut_sub_r
        record_r_mod['TRIM_QRY_L'] += cut_qry_r

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

    # # Adjust reference and contig span after trimming
    # record_l_mod['QRY_MAP_LEN'] = record_l_mod['QRY_END'] - record_l_mod['QRY_POS']
    # record_l_mod['SUB_MAP_LEN'] = record_l_mod['END'] - record_l_mod['POS']
    #
    # record_r_mod['QRY_MAP_LEN'] = record_r_mod['QRY_END'] - record_r_mod['QRY_POS']
    # record_r_mod['SUB_MAP_LEN'] = record_r_mod['END'] - record_r_mod['POS']

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

    last_no_match = False  # Continue until the last element is a match

    while index < index_end and (diff_cumulative <= diff_bp or last_no_match or len(trace_list) == 0):
        cigar_count += 1
        cigar_len, cigar_op = cigar_list[index]

        if cigar_op == '=':
            event_count = 0

            sub_bp = cigar_len
            qry_bp = cigar_len

            last_no_match = False

        elif cigar_op == 'X':
            event_count = cigar_len

            sub_bp = cigar_len
            qry_bp = cigar_len

            last_no_match = True

        elif cigar_op == 'I':
            event_count = 1

            sub_bp = 0
            qry_bp = cigar_len

            last_no_match = True

        elif cigar_op == 'D':
            event_count = 1

            sub_bp = cigar_len
            qry_bp = 0

            last_no_match = True

        elif cigar_op == 'S':
            event_count = 0

            sub_bp = 0
            qry_bp = 0

            clip_s_sum += cigar_len

            last_no_match = True

        elif cigar_op == 'H':
            event_count = 0

            sub_bp = 0
            qry_bp = 0

            clip_h_sum += cigar_len

            last_no_match = True

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
