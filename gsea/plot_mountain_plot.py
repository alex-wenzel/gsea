from matplotlib.pyplot import figure

from .plot.plot.decorate import decorate
from .plot.plot.save_plot import save_plot


def plot_mountain_plot(in_, y, cumulative_sums, title, file_path):
    """
    Plot mountain plot.
    Arguments:
        in_ (array):
        y (array):
        cumulative_sums (array):
        title (str):
        file_path (str):
    Returns:
        None
    """

    ax = figure(figsize=(12, 8)).gca()

    ax.plot(range(in_.size), in_, linewidth=1.8, color='#808080', alpha=0.16)
    ax.plot(range(in_.size), y, linewidth=2.6, color='#9017E6')
    ax.plot(range(in_.size), cumulative_sums, linewidth=1.8, color='#20D9BA')

    decorate(title=title, legend_loc=None)

    if file_path:
        save_plot(file_path)
