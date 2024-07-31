"""
Large variant discovery call methods
"""

import numpy as np
import pandas as pd

import pavlib

#
# Variant call objects
#

CALL_SOURCE = 'ALNTRUNC'

class Variant(object):
    """
    Base class for variant call objects.

    Attributes:
        chrom: Variant chromosome
        pos: Variant position
        end: Variant end
        variant_id: Variant ID
        vartype: Variant class
        svtype: Variant type
        svsubtype: Variant subtype (e.g. TANDEM)
        svlen: Variant length
        qry_region: Query region (pavlib.seq.Region or str)
        call_source: Type of evidence supporting the variant call
        seq: Variant sequence or "*" if the sequence is not set.
        is_rev: True if the query is on the reverse strand of the reference
        score_variant: Variant score

        hom_ref: Breakpoint homology in the reference (set by complete_anno())
        hom_qry: Breakpoint homology in the query (set by complete_anno())

        interval: Variant interval
        caller_resources: CallerResources object
        df_ref_trace: Variant reference trace

        is_complete_anno: Set when variant annotations are complete (only on accepted variants)
    """

    def __init__(self,
                 interval, caller_resources,
                 score_variant=-np.inf,
                 df_ref_trace=None
                 ):

        # Check and set
        if interval is None:
            raise RuntimeError('Cannot create variant with interval: None')

        self.interval = interval

        if caller_resources is None:
            raise RuntimeError('Cannot create variant with caller_resources: None')

        self.caller_resources = caller_resources

        try:
            self.score_variant = float(score_variant)
        except ValueError:
            raise RuntimeError(f'Variant parameter "score_variant" is not a number: {score_variant}')

        self.df_ref_trace = df_ref_trace

        self.strand = '-' if interval.is_rev else '+'

        # Variant call fields (most must be set by the variant type implementation)
        self.chrom = interval.chrom
        self.pos = None
        self.end = None
        self.vartype = None
        self.svtype = None
        self.svsubtype = None
        self.svlen = 0
        self.is_rev = interval.is_rev

        self.qry_region = None

        self.hom_ref = None
        self.hom_qry = None

        self.call_source = CALL_SOURCE
        self.seq = '*'

        # Set to True when annotations are complete
        self.is_complete_anno = False

    def __getattr__(self, name):
        """
        Get additional attributes by name.

        :param name: Attribute name.

        :return: Attribute value.
        """

        if name == 'variant_id':
            if self.is_null():
                return '*-*-NULLVAR-0'

            # Check for unset required fields
            unset_fields = list()

            for field in ['chrom', 'pos', 'svtype', 'svlen']:
                if getattr(self, field, None) is None:
                    unset_fields.append(field)

            if unset_fields:
                raise RuntimeError(f'Cannot complete annotations: Missing required variant fields: {", ".join(unset_fields)}')

            return f'{self.chrom}-{self.pos}-{str(self.svtype).upper()}-{self.svlen}'

        raise AttributeError(f'Variant has no attribute: {name}')

    def __repr__(self):
        """
        String representation of a variant object.

        :return: String representation.
        """

        return f'Variant({self.variant_id}, ref={self.interval.ref_region}, qry={self.interval.qry_region}, strand={self.strand})'

    def complete_anno(self):
        """
        Complete annotations on this variant call.

        Initial variant calls are not completed so they can be prioritized
        before spending the CPU cycles and IO time pulling sequences to complete annotations.

        Implementations of this method should complete annotations and set `self.is_complete_anno` to `True`.

        :param caller_resources: Caller resources including qry/ref sequences, k-mer utilities, and scoring models.
        """

        if self.is_complete_anno or self.is_null():
            return

        # Variant type implementation
        self.complete_anno_type()

        # Check for unset required fields
        unset_fields = list()

        for field in ['chrom', 'pos', 'end', 'vartype', 'svtype', 'svlen', 'is_rev', 'qry_region']:
            if getattr(self, field) is None:
                unset_fields.append(field)

        if unset_fields:
            raise RuntimeError(f'Cannot complete annotations: Missing required variant fields: {", ".join(unset_fields)}')

        self.is_complete_anno = True

    def is_null(self):
        """
        Determine if variant is a null call.

        :return: `True` if the variant call is Null (does not represent a real variant), and 'False' if the variant call
            exists.
        """

        return self.score_variant == -np.inf

    def complete_anno_type(self):
        pass

    def row(self):
        """
        Get a Pandas Series object representing a variant call table row for this variant call.

        Must be called after `complete_anno()`.

        :return: Pandas Series object of this variant call.
        """

        if self.is_null():
            raise RuntimeError('Cannot get row for a null variant')

        self.row_type()

    def row_impl(self):
        """
        Variant type implementation of row().

        :return: A Pandas Series object representing a variant call table row for this variant call.
        """
        raise NotImplementedError('Must be implemented by derived class')


# TODO: correct untemplated insertion breakpoints for homology at breakpoints
#
# Find untemplated insertion
#
# Trim reference if ends overlap around an untemplated insertion. This is likely homology at the breakpoint
# where the homologous sequence at each end of the insertion was aligned to the reference.
#
#                             [    HOM REGION    ]
#     Ref: --------->--------->--------->--------->--------->--------->--------->
# Align 1: --------->--------->--------->--------->
# Align 2:                     --------->--------->--------->--------->--------->
#
# Alignments 1 & 2 skip query sequence between them (untemplated insertion). This does not occur for strict templated
# insertions because query trimming will remove redundancy.

class InsertionVariant(Variant):
    """
    Insertion variant call.
    """

    # Simple insertion (INS = unaligned or aligned to another site)
    # INS:                           -------------------
    #                                        ||
    # Qry: --->--------->---------->---------> --------->-------->--------->--------->-----
    # Ref: --------------------------------------------------------------------------------

    def __init__(self, interval, caller_resources, df_ref_trace=None):
        Variant.__init__(self, interval, caller_resources, df_ref_trace=df_ref_trace)

        # Return immediately to leave the variant call a Null type
        if interval.seg_n != 1 or interval.len_qry == 0:
            return

        # Reference gap penalty: Penalize insertions for overlaps or deletions in the reference at the INS breakpoint
        len_ref = np.abs(self.interval.len_ref)

        self.ref_overlap = len_ref if self.interval.len_ref < 0 else 0
        self.svlen = len(interval.qry_region) + self.ref_overlap

        self.score_variant = \
            self.caller_resources.score_model.gap(self.svlen) + \
            (self.caller_resources.score_model.gap(len_ref) if len_ref > 0 else 0.0)

        self.pos = self.interval.qry_region.pos
        self.end = self.pos + 1
        self.vartype = 'SV' if self.svlen > 50 else 'INDEL'
        self.svtype = 'INS'


class DeletionVariant(Variant):
    """
    Deletion variant call.
    """

    # Simple deletion
    # Qry: ---->---------->--------->                  --------->-------->--------->-------
    # Ref: --------------------------------------------------------------------------------

    def __init__(self, interval, caller_resources, df_ref_trace=None):
        Variant.__init__(self, interval, caller_resources, df_ref_trace=df_ref_trace)

        # Return immediately to leave the variant call a Null type
        if interval.len_ref <= 0 or interval.seg_n != 0:
            return

        self.pos = self.interval.qry_region.pos
        self.end = self.interval.qry_region.end
        self.svlen = len(self.interval.qry_region)
        self.score_variant = self.caller_resources.score_model.gap(self.svlen)

        self.vartype = 'SV' if self.svlen > 50 else 'INDEL'
        self.svtype = 'DEL'


class InversionVariant(Variant):
    """
    Balanced inversion variant call.
    """

    def __init__(self, interval, caller_resources, df_ref_trace=None):
        Variant.__init__(self, interval, caller_resources, df_ref_trace=df_ref_trace)

        # Check segment - only allow inversion if there is a single inverted aligned segment
        if interval.seg_n != 1 or interval.df_segment.iloc[1]['IS_REV'] == interval.is_rev or interval.len_ref <= 0:
            return

        # Note: Breakpoints may not be placed consistently in inverted repeats on each end, the aligner makes an
        # alignment decision independently for each breakpoint. Therefore, the inverted segment may not align to the
        # reference gap between the repeats. An alignment penalty should be applied to discourage inversion calls
        # that do not fit the classic model (an inverted alignment in a reference gap), but not penalize alignment
        # this common alignment artifact.
        #
        #            [ > >  Repeat  > > ]                                        [ < <  Repeat  < < ]
        # Flank: --->--------->------->                                                          --->----
        #   INV:       <---------<---------<---------<---------<---------<---------
        # Trace:  NML    |   DUP      |                  INV                       |    DEL      |   NML
        #
        # The pattern above is apparent in trimmed alignments (by contig trimming). The un-trimmed alignment will
        # typically align through both repeats:
        #
        #            [ > >  Repeat  > > ]                                        [ < <  Repeat  < < ]
        # Flank: --->--------->------->--                                        --------->------->------
        #   INV:     --<---------<---------<---------<---------<---------<---------<---------<-------
        # Trace:  NML    |   DUP      |                  INV                       |    DEL      |   NML
        #
        # These two sources of information are combined. The tig-trimmed alignment is used to score the inversion after
        # redundant alignments are removed so complex patterns can be penalized appropriately. The untrimmed alignment
        # is used for setting inner- and outer-breakpoint locations.

        self.size_gap = np.abs(interval.len_ref - interval.df_segment.iloc[1]['LEN_REF'])
        self.svlen = interval.len_ref

        self.score_variant = (
            self.caller_resources.score_model.gap(self.size_gap) +  # Penalize differences between reference gap and inverted alignment
            self.caller_resources.score_model.template_switch * 2
        )

        # Set inner/outer reference regions
        self.region_ref_outer = pavlib.seq.Region(
            interval.chrom,
            interval.df_segment.iloc[1]['POS'] - interval.df_segment.iloc[1]['TRIM_REF_L'],
            interval.df_segment.iloc[1]['END'] + interval.df_segment.iloc[1]['TRIM_REF_R']
        )

        self.region_ref_inner = pavlib.seq.Region(
            interval.chrom,
            interval.df_segment.iloc[0]['END'] + interval.df_segment.iloc[0]['TRIM_REF_R'],
            interval.df_segment.iloc[-1]['POS'] - interval.df_segment.iloc[-1]['TRIM_REF_L']
        )

        # Set inner/outer query regions
        self.region_qry_outer = pavlib.seq.Region(
            interval.qry_id,
            interval.df_segment.iloc[1]['QRY_POS'] - interval.df_segment.iloc[1]['TRIM_QRY_L'],
            interval.df_segment.iloc[1]['QRY_END'] + interval.df_segment.iloc[1]['TRIM_QRY_R']
        )

        if interval.is_rev:
            row_l = interval.df_segment.iloc[-1]
            row_r = interval.df_segment.iloc[0]

            trim_qry_r = row_l['TRIM_QRY_L']
            trim_qry_l = row_r['TRIM_QRY_R']

        else:
            row_l = interval.df_segment.iloc[0]
            row_r = interval.df_segment.iloc[-1]

            trim_qry_r = row_r['TRIM_QRY_R']
            trim_qry_l = row_l['TRIM_QRY_L']

        self.region_qry_inner = pavlib.seq.Region(
            interval.qry_id,
            row_l['QRY_END'] + trim_qry_r,
            row_r['QRY_POS'] - trim_qry_l
        )

        # TODO: Verify query coordinates on an inversion with flanks in forward orientation

        self.pos = self.region_ref_outer.pos
        self.end = self.region_ref_outer.end
        self.svtype = 'INV'
        self.svlen = len(self.region_ref_outer)

        return


class TandemDuplicationVariant(Variant):
    """
    Tandem duplication variant call.
    """

    # Tandem duplication (TD)
    # Repeats:                                [> REP >]            [> REP >]
    # Qry 1:    --------->--------->--------->--------->--------->--------->
    # Qry 2:                                  --------->--------->--------->--------->--------->--------->
    #
    # Directly-oriented repeats may mediated tandem repeats. Look at alignment-trimming in three ways:
    # * Qry trimming: Identifies TD if redundant query bases are removed, but queries still overlap
    # * Qry & Ref trimming: Find a breakpoint for an insertion call.
    # * No trimming: Identify the repeats at the ends of the TD.
    #
    # The resulting call is an INS call with a breakpoint placed in a similar location if the TD was called as an
    # insertion event in the CIGAR string. The variant is annotated with the DUP locus, and the sequence of both copies
    # is

    def __init__(self, interval, caller_resources, df_ref_trace=None):
        Variant.__init__(self, interval, caller_resources, df_ref_trace=df_ref_trace)

        if interval.seg_n != 0 or interval.len_ref >= 0:
            return

        if interval.df_segment.iloc[0]['INDEX'] not in caller_resources.df_align_tigref:
            return

        if interval.df_segment.iloc[-1]['INDEX'] not in caller_resources.df_align_tigref:
            return

        #raise RuntimeError('Tandem duplication not implemented - Test and complete')

        # Determine left-most breakpoint using alignment trimming for homology
        self.pos = caller_resources.df_align_tigref.loc[interval.df_segment.iloc[0]['INDEX'], 'END']
        #sub_end = caller_resources.df_align_tigref.loc[interval.df_segment.iloc[-1]['INDEX'], 'POS']

        self.end = self.pos + 1

        qry_pos = caller_resources.df_align_tigref.loc[interval.df_segment.iloc[0]['INDEX'], 'QRY_END']
        qry_end = caller_resources.df_align_tigref.loc[interval.df_segment.iloc[-1]['INDEX'], 'QRY_POS']

        self.svlen = qry_end - qry_pos

        self.score_variant = \
            self.caller_resources.score_model.gap(self.svlen)

        self.svtype = 'INS'
        self.svsubtype = 'TD'


class ComplexVariant(Variant):
    """
    Complex variant call.
    """

    def __init__(self, interval, caller_resources, df_ref_trace=None):
        Variant.__init__(self, interval, caller_resources, df_ref_trace=df_ref_trace)

        # Compute variant score
        self.score_variant = caller_resources.score_model.template_switch * (interval.df_segment.shape[0] - 1)  # Template switches between segments

        for i in range(1, interval.df_segment.shape[0] - 1):  # Gap per segment
            self.score_variant += caller_resources.score_model.gap(
                interval.df_segment.loc[i, 'LEN_QRY']
            )

        #self.score_variant += caller_resources.score_model.gap(interval.len_ref) / 4  # Penalize reference gap between anchors

        self.pos = self.interval.ref_region.pos
        self.end = self.interval.ref_region.end

        self.svtype = 'CPX'
        self.svlen = len(interval.qry_region) + np.abs(interval.len_ref)

    def complete_anno_type(self):

        # Get reference trace
        if self.df_ref_trace is None:
            self.df_ref_trace = get_reference_trace(self.interval)

        self.struct_qry = get_qry_struct_str(self.df_segment)
        self.struct_ref = get_ref_struct_str(self.df_ref_trace)

        self.row = pd.Series(
            [
                self.chrom, self.ref_region.pos, self.ref_region.end,
                'CPX', self.svlen,
                str(self.qry_region),
                self.struct_ref, self.struct_qry,
                '-' if self.is_rev else '+',
            ],
            index=[
                '#CHROM', 'POS', 'END',
                'SVTYPE', 'SVLEN',
                'QRY_REGION',
                'CPX_REF', 'CPX_QRY',
                'STRAND',
            ]
        )


def get_reference_trace(interval):

    # Set local template switch boundaries
    # Distance from reference positions (pos & end) must be no more half of:
    #   * CPX event in contig bases.
    #   * Distance between reference anchoring breakpoints.
    # May not expand beyond archoring alignments.
    local_dist = max([len(interval.ref_region), len(interval.qry_region)]) // 2

    local_pos = max([
        interval.ref_region.pos - local_dist,
        min([
            interval.df_segment.iloc[0]['POS'],
            interval.df_segment.iloc[0]['END'],
            interval.df_segment.iloc[-1]['POS'],
            interval.df_segment.iloc[-1]['END']
        ])
    ])

    local_end = min([
        interval.ref_region.end + local_dist,
        max([
            interval.df_segment.iloc[0]['POS'],
            interval.df_segment.iloc[0]['END'],
            interval.df_segment.iloc[-1]['POS'],
            interval.df_segment.iloc[-1]['END']
        ])
    ])

    # Make depth table
    df_depth = pavlib.align.align_bed_to_depth_bed(interval.df_segment.loc[interval.df_segment['IS_ALIGNED']], df_fai=None)

    # Refine depth to local regions (remove distal, trim segments crossing local region pos & end)
    df_depth['POS'] = df_depth['POS'].apply(lambda val: max(val, local_pos))
    df_depth['END'] = df_depth['END'].apply(lambda val: min(val, local_end))

    df_depth = df_depth.loc[
        (df_depth['#CHROM'] == interval.ref_region.chrom) &
        (df_depth['END'] - df_depth['POS'] > 0)
    ]

    df_depth['INDEX'] = df_depth['INDEX'].fillna('').astype(str).apply(lambda val: val.strip())

    # Remove depth records that cover only anchor alignments at the flanks
    anchor_index_set = {
        str(interval.df_segment.iloc[0]['INDEX']),
        str(interval.df_segment.iloc[-1]['INDEX'])
    }

    while df_depth.iloc[-1]['INDEX'] in anchor_index_set and df_depth.iloc[-1]['DEPTH'] == 1:
        df_depth = df_depth.iloc[:-1]

    while df_depth.iloc[0]['INDEX'] in anchor_index_set and df_depth.iloc[0]['DEPTH'] == 1:
        df_depth = df_depth.iloc[1:]

    # Make reference context
    df_ref_trace_list = list()  # Type, length, depth fwd, depth rev, alignment index

    del_len = 0

    for index, row in df_depth.iterrows():

        depth_len = row['END'] - row['POS']
        depth = row['DEPTH']

        if depth > 0:
            index_set = {int(val) for val in row['INDEX'].split(',')}

            if len(index_set) != depth:
                raise RuntimeError(f'Depth record index list length {len(index_set)} does not match depth {depth}')

            fwd_count = np.sum(interval.df_segment.loc[interval.df_segment['INDEX'].isin(index_set), 'STRAND'] == '+')
            rev_count = depth - fwd_count

        else:
            fwd_count = 0
            rev_count = 0

        if interval.is_rev:
            fwd_count, rev_count = rev_count, fwd_count

        # Get segment type
        if depth == 0:
            seg_type = 'DEL'
            del_len += depth_len

        elif depth == 1:
            if rev_count == 1:
                seg_type = 'INV'
            else:
                seg_type = 'NML'

        else:
            dup_type = {
                0: 'DEL',
                1: 'NML',
                2: 'DUP',
                3: 'TRP',
                4: 'QUAD',
            }.get(depth, 'HDUP')  # Default to high-copy duplication

            if fwd_count == depth:
                seg_type = dup_type
            elif rev_count >= (depth - 1):
                seg_type = f'INV{dup_type}'
            else:
                seg_type = f'MIX{dup_type}'

        # Add trace record
        df_ref_trace_list.append(
            pd.concat([
                row,
                pd.Series(
                    [seg_type, depth_len, fwd_count, rev_count],
                    index=['TYPE', 'LEN', 'FWD_COUNT', 'REV_COUNT']
                )
            ], axis=0)
        )

    df_ref_trace = pd.concat(df_ref_trace_list, axis=1).T

    return df_ref_trace


def get_qry_struct_str(df_segment):
    """
    Get a comuplex SV structure following the reference through template switches, templated insertions, and
    untemplated insertions.

    :param df_segment: Segment table.

    :return: String describing the complex SV structure.
    """

    last_chrom = df_segment.iloc[0]['#CHROM']
    last_pos = df_segment.iloc[0]['END']

    struct_list = list()

    for i in range(1, df_segment.shape[0] - 1):
        row = df_segment.iloc[i]

        if row['IS_ALIGNED']:
            if row['#CHROM'] == last_chrom:
                dist = row['POS'] - last_pos
                struct_list.append(f'TS({dist})')
            else:
                struct_list.append(f'TS({row["#CHROM"]})')
                last_chrom = row['#CHROM']

            struct_list.append(f'TINS({row["STRAND"]}{row["LEN_QRY"]})')

            last_pos = row['END']
        else:
            struct_list.append(f'INS({row["LEN_REF"]})')

    row = df_segment.iloc[-1]
    dist = row['POS'] - last_pos
    struct_list.append(f'TS({dist})')

    return ':'.join(struct_list)


def get_ref_struct_str(df_ref_trace):
    """
    Get reference structure string describing a complex SV from the reference perspective.

    :param df_ref_trace: Reference trace table.

    :return: A string describing the reference structure.
    """

    return ':'.join(df_ref_trace.apply(lambda row: f'{row["TYPE"]}({row["LEN"]})', axis=1))
