from matplotlib.pyplot import figure

from .plot.plot.decorate import decorate
from .plot.plot.save_plot import save_plot


def plot_mountain_plot(in_, cumulative_sums, score, title, file_path):
    """
    Plot mountain plot.
    Arguments:
        in_ (array):
        cumulative_sums (array):
        score (number):
        title (str):
        file_path (str):
    Returns:
        None
    """

    ax = figure(figsize=(12, 8)).gca()

    ax.plot(range(in_.size), in_, linewidth=1.8, color='#2E211B', alpha=0.8)
    ax.plot(range(in_.size), cumulative_sums, linewidth=2.6, color='#20D9BA')

    decorate(
        title='{}\nEnrichment Score = {}'.format(title, score),
        legend_loc=None)

    if file_path:
        save_plot(file_path)
