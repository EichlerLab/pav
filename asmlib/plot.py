"""
Make plots of inversion calls.
"""

import numpy as np

import asmlib
import kanapy


def dotplot_inv_call(inv_call, ref_fa, aln_file_name, coords='discovery'):
    """
    Make a dotplot of an inversion call.

    :param inv_call: asmlib.inv.InvCall object describing the inversion.
    :param ref_fa: Reference FASTA.
    :param aln_file_name: Alignment file name.
    :param coords: Coordinates to plot. Defaults to discovery region, which includes unique sequence on each side
        of the inversion and shows the full context of the inversion call. Other valid values are "outer" and "inner"
        for outer and inner breakpoints.


    :return: A plot object. Before it is discarded, this object should be closed with `matplotlib.pyplot.close()` to
        free memory.
    """

    # Get regions
    if coords == 'discovery':
        region_ref = inv_call.region_ref_discovery
        region_tig = inv_call.region_tig_discovery

    elif coords == 'outer':
        region_ref = inv_call.region_ref_outer,
        region_tig = inv_call.region_tig_outer

    elif coords == 'inner':
        region_ref = inv_call.region_ref_inner
        region_tig = inv_call.region_tig_inner

    else:
        raise RuntimeError('Unrecognized coords arguments: {}'.format(coords))

    # Get reference sequence
    seq_ref = asmlib.seq.fa_region(region_ref, ref_fa)

    # Get contig sequence
    contig_record_list = asmlib.seq.get_matching_alignments(region_ref, region_tig, aln_file_name, ref_fa)

    if len(contig_record_list) != 1:
        raise RuntimeError('Expected 1 overlapping region for {} aligned to {}: Found {}'.format(
            region_tig, region_ref, len(contig_record_list)
        ))

    contig_record = contig_record_list[0]
    del contig_record_list

    # Count hard-clipped bases
    query_start = 0

    index = 0
    while contig_record.cigar[index][0] == 5:  # 5 = H cigar op
        query_start += contig_record.cigar[index][1]
        index += 1

    seq_tig = contig_record.query_sequence[(region_tig.pos - query_start):(region_tig.end - query_start)]

    # Create plot config
    plot_config = {
        'label_x': '{} ({:,d} - {:,d})'.format(region_tig.chrom, region_tig.pos + 1, region_tig.end),
        'label_y': '{} ({:,d} - {:,d})'.format(region_ref.chrom, region_ref.pos + 1, region_ref.end),
        'start_x': region_tig.pos,
        'start_y': region_ref.pos,
    }

    # Add annotations
    anno_list = [
        {
            'type': 'hshade',
            'y1': inv_call.region_ref_outer.pos + 1,  # From BED to 1-based closed coordinates
            'y2': inv_call.region_ref_outer.end,
            'color': 'lightseagreen'
        },
        {
            'type': 'vline',
            'x': np.array([inv_call.region_tig_outer.pos, inv_call.region_tig_outer.end]),
            'ymin': inv_call.region_ref_discovery.pos + 1,
            'ymax': inv_call.region_ref_discovery.end,
            'color': 'lightseagreen',
            'style': 'solid',
            'alpha': 0.4
        }
    ]


    if inv_call.region_ref_outer != inv_call.region_ref_inner:
        anno_list.append(
            {
                'type': 'hshade',
                'y1': inv_call.region_ref_inner.pos + 1,  # From BED to 1-based closed coordinates
                'y2': inv_call.region_ref_inner.end,
                'color': 'seagreen'
            }
        )

        anno_list.append(
            {
                'type': 'vline',
                'x': np.array([inv_call.region_tig_inner.pos, inv_call.region_tig_inner.end]),
                'ymin': inv_call.region_ref_discovery.pos + 1,
                'ymax': inv_call.region_ref_discovery.end,
                'color': 'seagreen',
                'style': 'dashed',
                'alpha': 0.4
            }
        )

    # Make plot
    fig = kanapy.plot.dotplot.dotplot(
        seq_x=seq_tig,
        seq_y=seq_ref,
        config=plot_config,
        title=inv_call.id,
        anno_list=anno_list
    )

    # Return figure
    return fig