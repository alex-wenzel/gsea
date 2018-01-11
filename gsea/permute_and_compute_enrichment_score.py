from numpy import empty
from numpy.random import shuffle

from .compute_enrichment_score import compute_enrichment_score


def permute_and_compute_enrichment_score(gene_scores,
                                         gene_set_genes,
                                         n_permutation,
                                         normalization_method='rank',
                                         power=1,
                                         statistic='ks'):
    """
    Compute how much permuted gene scores enrich gene-set genes.
    Arguments:
        gene_scores (Series): (n_gene_with_score)
        gene_set_genes (iterable): (n_gene)
        n_permutation (int):
        normalization_method (str): '-0-' | '0-1' | 'rank'
        power (number): power to raise gene_score
        statistic (str): 'ks' (Kolmogorov-Smirnov) | 'auc' (area under curve)
    Returns:
        array: (n_permutation)
    """

    enrichment_scores = empty(n_permutation)

    for i in range(n_permutation):

        shuffle(gene_scores)

        enrichment_scores[i] = compute_enrichment_score(
            gene_scores,
            gene_set_genes,
            normalization_method=normalization_method,
            power=power,
            statistic=statistic)

    return enrichment_scores
