from matplotlib.pyplot import figure, gca, show
from numpy import abs, asarray, in1d, where


def compute_enrichment_score(gene_scores,
                             gene_set_genes,
                             power=1,
                             statistic='auc',
                             plot=False):
    """
    Compute how much gene scores enrich gene-set genes.
    Arguments:
        gene_scores (Series): (n_genes_with_score); sorted and indexed by gene
        gene_set_genes (iterable): (n_genes)
        power (number): power to raise gene_scores
        statistic (str): 'auc' (area under curve) | 'ks' (Kolmogorov-Smirnov)
        plot (bool): whether to plot the mountain plot
    Returns:
        float: enrichment score
    """

    gene_scores = gene_scores.sort_values(ascending=False)

    in_ = in1d(gene_scores.index, gene_set_genes, assume_unique=True)

    if power != 1:
        gene_scores = abs(asarray(gene_scores))**power

    in_int = in_.astype(int)
    hit = (gene_scores * in_int) / gene_scores[in_].sum()
    miss = (1 - in_int) / (in_.size - in_.sum())
    y = hit - miss

    cumulative_sums = y.cumsum()

    if statistic == 'auc':
        enrichment_score = cumulative_sums.sum()

    elif statistic == 'ks':
        max_ = cumulative_sums.max()
        min_ = cumulative_sums.min()
        enrichment_score = where(abs(min_) < abs(max_), max_, min_)

    else:
        raise ValueError('Unknown statistic: {}.'.format(statistic))

    if plot:
        figure(figsize=(8, 5))
        ax = gca()
        ax.plot(range(in_.size), in_, color='#808080', alpha=0.16)
        ax.plot(range(in_.size), y, color='#9017E6')
        ax.plot(range(in_.size), cumulative_sums, color='#20D9BA')
        show()

    return enrichment_score
