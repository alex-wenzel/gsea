from numpy import abs, asarray, in1d, where


def compute_enrichment_score(gene_scores,
                             gene_set_genes,
                             power=1,
                             statistic='AUC'):
    """
    Compute how much gene_scores enriches gene_set_genes.
    Arguments:
        gene_scores (Series): (n_genes_with_score); sorted and gene indexed
        gene_set_genes (iterable): (n_genes)
        power (number): Power to raise gene_scores
        statistic (str): 'AUC' (Area Under Curve) | 'KS' (Kolmogorov-Smirnov)
    Returns:
        float: Enrichment score
    """

    gene_scores = gene_scores.sort_values(ascending=False)

    # Check if gene_scores genes are in gene_set_genes
    in_ = in1d(gene_scores.index, gene_set_genes, assume_unique=True)
    in_int = in_.astype(int)

    if power != 1:
        gene_scores = abs(asarray(gene_scores))**power

    hit = (gene_scores * in_int) / gene_scores[in_].sum()
    miss = (1 - in_int) / (in_.size - in_.sum())
    y = hit - miss

    # Compute enrichment score
    cumulative_sums = y.cumsum()

    if statistic == 'AUC':
        enrichment_score = cumulative_sums.sum()

    elif statistic == 'KS':
        max_ = cumulative_sums.max()
        min_ = cumulative_sums.min()
        enrichment_score = where(abs(min_) < abs(max_), max_, min_)

    else:
        raise ValueError('Unknown statistic: {}.'.format(statistic))

    # TODO: Plot
    # import matplotlib as mpl
    # mpl.pyplot.figure(figsize=(8, 5))
    # ax = mpl.pyplot.gca()
    # ax.plot(range(in_.size), in_, color='#808080', alpha=0.16)
    # ax.plot(range(in_.size), y, color='#9017E6')
    # ax.plot(range(in_.size), cumulative_sums, color='#20D9BA')
    # mpl.pyplot.show()

    return enrichment_score
