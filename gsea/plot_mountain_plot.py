from matplotlib.pyplot import figure
from seaborn import rugplot

from .plot.plot.decorate import decorate
from .plot.plot.save_plot import save_plot


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

    ax = figure(figsize=(16, 8)).gca()

    ax.plot(range(in_.size), cumulative_sums, linewidth=2.6, color='#20D9BA')
    rugplot(
        [i for i in range(in_.size) if in_[i]],
        ax=ax,
        height=0.08,
        linewidth=1.08,
        color='#9017E6',
        alpha=1)

    decorate(
        title='{}\nEnrichment Score = {:.3f}'.format(title, float(score)),
        legend_loc=None)

    if file_path:
        save_plot(file_path)
