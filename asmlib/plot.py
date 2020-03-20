"""
Make plots of inversion calls.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

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

def kmer_density_plot(inv_call, hap=None, width=7, height=4, dpi=300):
    """
    Get plot of k-mer density from the called region.

    :param inv_call: Inversion call
    :param width: Plot width (in).
    :param height: Plot height (in).
    :param dpi: DPI

    :return: A plot object. Before it is discarded, this object should be closed with `matplotlib.pyplot.close()` to
        free memory.
    """
    
    fig, ax1, ax2 = kmer_density_plot_base(inv_call.df, width=width, height=height, dpi=dpi)
    

    # Add lines
    if inv_call.region_tig_outer is not None:
        ax1.vlines(
            x=[inv_call.region_tig_outer.pos, inv_call.region_tig_outer.end],
            ymin=-1, ymax=1,
            color='lightseagreen',
            linestyles='solid'
        )

        ax2.vlines(
            x=[inv_call.region_tig_outer.pos, inv_call.region_tig_outer.end],
            ymin=0, ymax=1,
            color='lightseagreen',
            linestyles='solid'
        )

    if inv_call.region_tig_inner is not None and (inv_call.region_tig_outer is None or inv_call.region_tig_outer != inv_call.region_tig_inner):
        ax1.vlines(
            x=[inv_call.region_tig_inner.pos, inv_call.region_tig_inner.end],
            ymin=-1, ymax=1,
            color='seagreen',
            linestyles='dashed'
        )

        ax2.vlines(
            x=[inv_call.region_tig_inner.pos, inv_call.region_tig_inner.end],
            ymin=0, ymax=1,
            color='seagreen',
            linestyles='dashed'
        )

    # Plot aestetics

    plot_title = inv_call.id

    if hap is not None:
        plot_title += ' ({})'.format(hap)

    fig.suptitle(plot_title)

    fig.subplots_adjust(bottom=0.15)  # Test/tune - Add space to labels don't run off bottom end

    # Return figure
    return fig


def kmer_density_plot_base(df, region_tig, width=7, height=4, dpi=300):
    
    # Make figure
    fig = plt.figure(figsize=(width, height), dpi=dpi)

    ax1, ax2 = fig.subplots(2, 1)


    ## Smoothed state (top pane) ##

    # Points
    ax1.scatter(
        x=np.asarray(df.loc[df['STATE'] == 0, 'INDEX']) + region_tig.pos,
        y=np.repeat(1, np.sum(df['STATE'] == 0)),
        color='blue',
        alpha=0.2
    )

    ax1.scatter(
        x=np.asarray(df.loc[df['STATE'] == 1, 'INDEX']) + region_tig.pos,
        y=np.repeat(0, np.sum(df['STATE'] == 1)),
        color='purple',
        alpha=0.2
    )

    ax1.scatter(
        x=np.asarray(df.loc[df['STATE'] == 2, 'INDEX']) + region_tig.pos,
        y=np.repeat(-1, np.sum(df['STATE'] == 2)),
        color='red',
        alpha=0.2
    )

    # Max density line (smoothed state call)
    ax1.plot(
        df['INDEX'] + region_tig.pos,
        (df['KERN_MAX'] - 1) * -1,
        color='black'
    )

    # Plot aestetics
    ax1.get_xaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    )

    ax1.set_yticks(np.asarray([-1, 0, 1]))
    ax1.set_yticklabels(np.array(['Rev', 'Fwd+Rev', 'Fwd']))



    ## Density (bottom pane) ##

    ax2.plot(
        df['INDEX'] + region_tig.pos,
        df['KERN_FWD'],
        color='blue'
    )

    ax2.plot(
        df['INDEX'] + region_tig.pos,
        df['KERN_FWDREV'],
        color='purple'
    )

    ax2.plot(
        df['INDEX'] + region_tig.pos,
        df['KERN_REV'],
        color='red'
    )

    # Plot aestetics

    ax2.get_xaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    )

    ax2.set_xlabel('{} ({:,d} - {:,d})'.format(
        region_tig.chrom,
        region_tig.pos + 1,
        region_tig.end
    ))

    ax1.tick_params(labelbottom=False)

    for label in ax2.get_xticklabels():
        label.set_rotation(30)
        label.set_ha('right')


    ## Return plot and axes ##
    return fig, ax1, ax2
