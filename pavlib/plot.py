"""
Make plots of inversion calls.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np

import pavlib
import kanapy


def dotplot_inv_call(
    inv_call, ref_fa, tig_fa=None, seq_tig=None,
    region_ref=None, region_tig=None,
    k=32
):
    """
    Make a dotplot of an inversion call.

    :param inv_call: pavlib.inv.InvCall object describing the inversion.
    :param ref_fa: Reference FASTA.
    :param tig_fa: Contig FASTA file name.
    :param seq_tig: Contig alignment sequence or `None`. If `None`, then `tig_fa` must be set.
    :param region_ref: Reference region. If None, uses inv_call.region_ref_discovery.
    :param region_tig: Contig region. If none, uses inv_call.region_tig_discovery.
    :param k: K-mer size.

    :return: A plot object. Before it is discarded, this object should be closed with `matplotlib.pyplot.close()` to
        free memory.
    """

    if region_ref is None:
        region_ref = inv_call.region_ref_discovery

    if region_tig is None:
        region_tig = inv_call.region_tig_discovery

    # Get reference sequence
    seq_ref = pavlib.seq.region_seq_fasta(region_ref, ref_fa, False)

    # Get contig sequence
    if seq_tig is None:

        if tig_fa is None:
            raise RuntimeError('Cannot get contig sequence: tig_fa is None')

        seq_tig = pavlib.seq.region_seq_fasta(region_tig, tig_fa)

    # Create plot config
    plot_config = {
        'label_x': '{} ({:,d} - {:,d})'.format(region_tig.chrom, region_tig.pos + 1, region_tig.end),
        'label_y': '{} ({:,d} - {:,d})'.format(region_ref.chrom, region_ref.pos + 1, region_ref.end),
        'start_x': region_tig.pos,
        'start_y': region_ref.pos,
        'k': k
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
            'ymin': region_ref.pos + 1,
            'ymax': region_ref.end,
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
                'ymin': region_ref.pos + 1,
                'ymax': region_ref.end,
                'color': 'seagreen',
                'style': 'dashed',
                'alpha': 0.4
            }
        )

    # Add short dotted lines for discovery region - Tig
    if inv_call.region_tig_discovery is not None:
        if region_tig.pos <= inv_call.region_tig_discovery.pos <= region_tig.end:
            anno_list.append(
                {
                    'type': 'vline',
                    'x': inv_call.region_tig_discovery.pos,
                    'ymin': region_ref.pos + 1,
                    'ymax': region_ref.pos + 1 + int(len(region_ref) * 0.1),
                    'color': '0.3',
                    'style': 'dotted',
                    'alpha': 0.4
                }
            )

            anno_list.append(
                {
                    'type': 'vline',
                    'x': inv_call.region_tig_discovery.pos,
                    'ymin': region_ref.end - int(len(region_ref) * 0.1),
                    'ymax': region_ref.end,
                    'color': '0.3',
                    'style': 'dotted',
                    'alpha': 0.4
                }
            )

        if region_tig.pos <= inv_call.region_tig_discovery.end <= region_tig.end:
            anno_list.append(
                {
                    'type': 'vline',
                    'x': inv_call.region_tig_discovery.end,
                    'ymin': region_ref.pos + 1,
                    'ymax': region_ref.pos + 1 + int(len(region_ref) * 0.1),
                    'color': '0.3',
                    'style': 'dotted',
                    'alpha': 0.4
                }
            )

            anno_list.append(
                {
                    'type': 'vline',
                    'x': inv_call.region_tig_discovery.end,
                    'ymin': region_ref.end - int(len(region_ref) * 0.1),
                    'ymax': region_ref.end,
                    'color': '0.3',
                    'style': 'dotted',
                    'alpha': 0.4
                }
            )
    
    if inv_call.region_ref_discovery is not None:
        if region_ref.pos <= inv_call.region_ref_discovery.pos <= region_ref.end:
            anno_list.append(
                {
                    'type': 'hline',
                    'y': inv_call.region_ref_discovery.pos,
                    'xmin': region_tig.pos + 1,
                    'xmax': region_tig.pos + 1 + int(len(region_tig) * 0.1),
                    'color': '0.3',
                    'style': 'dotted',
                    'alpha': 0.4
                }
            )

            anno_list.append(
                {
                    'type': 'hline',
                    'y': inv_call.region_ref_discovery.pos,
                    'xmin': region_tig.end - int(len(region_tig) * 0.1),
                    'xmax': region_tig.end,
                    'color': '0.3',
                    'style': 'dotted',
                    'alpha': 0.4
                }
            )

        if region_ref.pos <= inv_call.region_ref_discovery.end <= region_ref.end:
            anno_list.append(
                {
                    'type': 'hline',
                    'y': inv_call.region_ref_discovery.end,
                    'xmin': region_tig.pos + 1,
                    'xmax': region_tig.pos + 1 + int(len(region_tig) * 0.1),
                    'color': '0.3',
                    'style': 'dotted',
                    'alpha': 0.4
                }
            )

            anno_list.append(
                {
                    'type': 'hline',
                    'y': inv_call.region_ref_discovery.end,
                    'xmin': region_tig.end - int(len(region_tig) * 0.1),
                    'xmax': region_tig.end,
                    'color': '0.3',
                    'style': 'dotted',
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


def kmer_density_plot(inv_call, hap=None, width=7, height=4, dpi=300, flank_whiskers=True):
    """
    Get plot of k-mer density from the called region.

    :param inv_call: Inversion call
    :param width: Plot width (in).
    :param height: Plot height (in).
    :param dpi: DPI
    :param flank_whiskers: If `True`, show whiskers above or below points indicating if they match the upstream or
        downstream flanking inverted duplication.

    :return: A plot object. Before it is discarded, this object should be closed with `matplotlib.pyplot.close()` to
        free memory.
    """
    
    fig = kmer_density_plot_base(
        inv_call.df, inv_call.region_tig_discovery,
        width=width, height=height, dpi=dpi,
        flank_whiskers=flank_whiskers
    )

    ax1, ax2 = fig.axes

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


def kmer_density_plot_base(df, region_tig, width=7, height=4, dpi=300, flank_whiskers=True):
    """
    Base k-mer density plot using a k-mer density DataFrame.

    :param df: Inversion call DataFrame.
    :param region_tig: Contig region plot is generated over.
    :param width: Figure width.
    :param height: Figure height.
    :param dpi: Figure DPI.
    :param flank_whiskers: If `True`, show whiskers above or below points indicating if they match the upstream or
        downstream flanking inverted duplication.

    :return: Plot figure object.
    """
    
    # Make figure
    fig = plt.figure(figsize=(width, height), dpi=dpi)

    ax1, ax2 = fig.subplots(2, 1)


    ## Smoothed state (top pane) ##

    # Flanking DUP whiskers
    if flank_whiskers:
        for index, subdf in df.loc[
            ~ pd.isnull(df['MATCH'])
        ].groupby(
            ['STATE', 'MATCH']
        ):

            # Whisker y-values
            ymin, ymax = sorted(
                [
                    -1 * (subdf.iloc[0]['STATE_MER'] - 1),
                    -1 * (subdf.iloc[0]['STATE_MER'] - 1) + (0.25 if index[1] == 'SAME' else -0.25)
                ]
            )

            # Plot whiskers
            ax1.vlines(
                x=subdf['INDEX'] + region_tig.pos,
                ymin=ymin, ymax=ymax,
                color='dimgray',
                linestyles='solid',
                linewidth=0.5
            )

    # Points
    ax1.scatter(
        x=np.asarray(df.loc[df['STATE_MER'] == 0, 'INDEX']) + region_tig.pos,
        y=np.repeat(1, np.sum(df['STATE_MER'] == 0)),
        color='blue',
        alpha=0.2
    )

    ax1.scatter(
        x=np.asarray(df.loc[df['STATE_MER'] == 1, 'INDEX']) + region_tig.pos,
        y=np.repeat(0, np.sum(df['STATE_MER'] == 1)),
        color='purple',
        alpha=0.2
    )

    ax1.scatter(
        x=np.asarray(df.loc[df['STATE_MER'] == 2, 'INDEX']) + region_tig.pos,
        y=np.repeat(-1, np.sum(df['STATE_MER'] == 2)),
        color='red',
        alpha=0.2
    )

    # Max density line (smoothed state call)
    ax1.plot(
        df['INDEX'] + region_tig.pos,
        (df['STATE'] - 1) * -1,
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
    return fig
