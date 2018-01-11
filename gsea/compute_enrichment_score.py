from matplotlib.pyplot import figure, gca, show
from numpy import abs, asarray, in1d, where
from pandas import Series

from .nd_array.nd_array.normalize_1d_array import normalize_1d_array


def compute_enrichment_score(gene_scores,
                             gene_set_genes,
                             normalization_method='rank',
                             power=1,
                             statistic='ks',
                             plot=False):
    """
    Compute how much gene scores enrich gene-set genes.
    Arguments:
        gene_scores (Series): (n_gene_with_score)
        gene_set_genes (iterable): (n_gene)
        normalization_method (str): '-0-' | '0-1' | 'rank'
        power (number): power to raise gene_scores
        statistic (str): 'ks' (Kolmogorov-Smirnov) | 'auc' (area under curve)
        plot (bool): whether to plot
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

    # TODO: modularize plot function
    if plot:
        figure(figsize=(8, 5))
        ax = gca()
        ax.plot(range(in_.size), in_, color='#808080', alpha=0.16)
        ax.plot(range(in_.size), y, color='#9017E6')
        ax.plot(range(in_.size), cumulative_sums, color='#20D9BA')
        show()

    return enrichment_score
