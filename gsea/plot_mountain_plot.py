from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import figure, subplot
from numpy import argmax
from seaborn import swarmplot

from .plot.plot.decorate import decorate
from .plot.plot.get_ax_positions_relative_to_ax import \
    get_ax_positions_relative_to_ax
from .plot.plot.save_plot import save_plot


# TODO: plot gene scores
def plot_mountain_plot(cumulative_sums, in_, score, title, file_path):
    """
    Plot mountain plot.
    Arguments:
        cumulative_sums (array):
        in_ (array):
        score (number):
        title (str):
        file_path (str):
    Returns:
        None
    """

    figure(figsize=(16, 8))

    gridspec = GridSpec(100, 1)

    ax = subplot(gridspec[:80, :])
    ax_bottom = subplot(gridspec[80:, :], sharex=ax)

    ax.plot(range(in_.size), cumulative_sums, linewidth=2.6, color='#20D9BA')
    ax.fill_between(
        range(in_.size), cumulative_sums, color='#EBF6F7', alpha=0.8)

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
        'ES = {:.3f}'.format(float(score)),
        size=18,
        weight='bold',
        color='#FF1968',
        horizontalalignment='center')

    decorate(
        ax=ax,
        despine_kwargs={
            'bottom': True,
        },
        style='white',
        title=title,
        ylabel='Enrichment Score')

    ax.set_xticklabels(
        [int(float(text.get_text()) + 1) for text in ax.get_xticklabels()])

    swarmplot(
        x=[i for i in range(in_.size) if in_[i]],
        ax=ax_bottom,
        size=8,
        color='#9017E6',
        alpha=0.92)
    decorate(
        ax=ax_bottom,
        despine_kwargs=dict(left=True),
        xlabel='Gene Rank',
        yticks=[])

    if file_path:
        save_plot(file_path)
