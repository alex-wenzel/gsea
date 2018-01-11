from numpy import abs, asarray, in1d, where
from pandas import Series

from .nd_array.nd_array.normalize_1d_array import normalize_1d_array
from .plot_mountain_plot import plot_mountain_plot


def compute_enrichment_score(gene_scores,
                             gene_set_genes,
                             normalization_method='rank',
                             power=1,
                             statistic='ks',
                             plot=False,
                             title=None,
                             plot_file_path=None):
    """
    Compute how much gene scores enrich gene-set genes.
    Arguments:
        gene_scores (Series): (n_gene_with_score)
        gene_set_genes (iterable): (n_gene)
        normalization_method (str): '-0-' | '0-1' | 'rank'
        power (number): power to raise gene_scores
        statistic (str): 'ks' (Kolmogorov-Smirnov) | 'auc' (area under curve)
        plot (bool): whether to plot
        title (str):
        plot_file_path (str):
    Returns:
        float: enrichment score
    """

    if normalization_method:
        gene_scores = Series(
            normalize_1d_array(gene_scores, normalization_method),
            name=gene_scores.name,
            index=gene_scores.index)

    gene_scores = gene_scores.sort_values(ascending=False)

    in_ = in1d(gene_scores.index, gene_set_genes, assume_unique=True)

    if power != 1:
        gene_scores = abs(asarray(gene_scores))**power

    in_int = in_.astype(int)
    hit = (gene_scores * in_int) / gene_scores[in_].sum()
    miss = (1 - in_int) / (in_.size - in_.sum())
    y = hit - miss

    cumulative_sums = y.cumsum()

    if statistic == 'ks':
        max_ = cumulative_sums.max()
        min_ = cumulative_sums.min()
        enrichment_score = where(abs(min_) < abs(max_), max_, min_)

    elif statistic == 'auc':
        enrichment_score = cumulative_sums.sum()

    else:
        raise ValueError('Unknown statistic: {}.'.format(statistic))

    if plot:

        plot_mountain_plot(in_, y, cumulative_sums, (
            title,
            gene_scores.name, )[title is None], plot_file_path)

    return enrichment_score
