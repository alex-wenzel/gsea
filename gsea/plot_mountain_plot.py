from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import figure, setp, subplot
from numpy import argmax
from seaborn import swarmplot

from .plot.plot.decorate import decorate
from .plot.plot.get_ax_positions_relative_to_ax import \
    get_ax_positions_relative_to_ax
from .plot.plot.save_plot import save_plot


def plot_mountain_plot(cumulative_sums,
                       hits,
                       gene_score,
                       title='GSEA Mountain Plot',
                       file_path=None):
    """
    Plot mountain plot.
    Arguments:
        cumulative_sums (array):
        hits (array):
        gene_score (array):
        title (str):
        file_path (str):
    Returns:
        None
    """

    figure(figsize=(16, 8))

    gridspec = GridSpec(100, 1)

    dividers = (
        70,
        85, )
    ax = subplot(gridspec[:dividers[0], :])
    ax_swarmplot = subplot(gridspec[dividers[0]:dividers[1], :], sharex=ax)
    ax_gene_score = subplot(gridspec[dividers[1]:, :], sharex=ax)

    x_grid = range(cumulative_sums.size)
    linewidth = 2.6
    fill_between_alpha = 0.26

    color = '#20D9BA'
    ax.plot(x_grid, cumulative_sums, linewidth=linewidth, color=color)
    ax.fill_between(
        x_grid, cumulative_sums, color=color, alpha=fill_between_alpha)

    max_x, max_y = argmax(cumulative_sums.values), cumulative_sums.max()
    ax.plot(
        max_x,
        max_y,
        marker='o',
        markersize=18,
        markeredgewidth=1.8,
        color='#FF1968',
        markeredgecolor='#EBF6F7')
    ax_y_min, ax_y_max = get_ax_positions_relative_to_ax()[2:]
    ax.text(
        max_x,
        max_y + 0.026 * (ax_y_max - ax_y_min),
        '{:.3f}'.format(max_y),
        size=18,
        weight='bold',
        color='#FF1968',
        horizontalalignment='center')
    decorate(ax=ax, style='white', title=title, ylabel='Enrichment Score')
    setp(ax.get_xaxis(), visible=False)

    swarmplot(
        x=[i for i in x_grid if hits[i]],
        ax=ax_swarmplot,
        size=8,
        color='#4E40D9',
        clip_on=False)
    decorate(
        ax=ax_swarmplot,
        despine_kwargs={
            'left': True,
            'bottom': True,
        },
        yticks=())
    setp(ax_swarmplot.get_xaxis(), visible=False)

    color = '#9017E6'
    ax_gene_score.plot(
        range(len(gene_score)), gene_score, linewidth=linewidth, color=color)
    ax_gene_score.fill_between(
        range(len(gene_score)),
        gene_score,
        color=color,
        alpha=fill_between_alpha)
    decorate(ax=ax_gene_score, xlabel='Gene Rank', ylabel='Gene Score')

    ax.set_xticklabels(
        [int(float(text.get_text()) + 1) for text in ax.get_xticklabels()])

    if file_path:
        save_plot(file_path)
