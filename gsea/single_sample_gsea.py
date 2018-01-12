from numpy import abs, asarray, in1d, where
from pandas import Series

from .nd_array.nd_array.normalize_1d_array import normalize_1d_array
from .plot_mountain_plot import plot_mountain_plot


def single_sample_gsea(gene_score,
                       gene_set_genes,
                       normalization_method='rank',
                       power=1,
                       statistic='ks',
                       plot=False,
                       title=None,
                       plot_file_path=None):
    """
    Single-sample GSEA.
    Arguments:
        gene_score (Series): (n_gene_with_score)
        gene_set_genes (iterable): (n_gene)
        normalization_method (str): '-0-' | '0-1' | 'rank'
        power (number): power to raise gene_score
        statistic (str): 'ks' (Kolmogorov-Smirnov) | 'auc' (area under curve)
        plot (bool): whether to plot
        title (str):
        plot_file_path (str):
    Returns:
        float: score
    """

    if normalization_method:

        gene_score = Series(
            normalize_1d_array(gene_score, normalization_method),
            name=gene_score.name,
            index=gene_score.index)

    gene_score.sort_values(ascending=False, inplace=True)

    in_ = in1d(gene_score.index, gene_set_genes.dropna(), assume_unique=True)

    if power != 1:
        gene_score = abs(asarray(gene_score))**power

    in_int = in_.astype(int)
    hit = (gene_score * in_int) / gene_score[in_].sum()
    miss = (1 - in_int) / (in_.size - in_.sum())
    y = hit - miss

    cumulative_sums = y.cumsum()

    if statistic == 'ks':
        max_ = cumulative_sums.max()
        min_ = cumulative_sums.min()
        score = where(abs(min_) < abs(max_), max_, min_)

    elif statistic == 'auc':
        score = cumulative_sums.sum()

    else:
        raise ValueError('Unknown statistic: {}.'.format(statistic))

    if plot:
        plot_mountain_plot(cumulative_sums, in_, score, (
            title,
            gene_score.name, )[title is None], plot_file_path)

    return score
