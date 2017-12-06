from numpy import empty
from numpy.random import shuffle

from .compute_enrichment_score import compute_enrichment_score


def permute_and_compute_enrichment_score(gene_scores,
                                         gene_set_genes,
                                         n_permutations,
                                         power=1,
                                         statistic='auc'):
    """
    Compute how much permuted gene scores enrich gene-set genes.
    Arguments:
        gene_scores (Series): (n_genes_with_score); sorted and indexed by gene
        gene_set_genes (iterable): (n_genes)
        power (number): power to raise gene_scores
        statistic (str): 'auc' (area under curve) | 'ks' (Kolmogorov-Smirnov)
    Returns:
        array: (n_permutations)
    """

    enrichment_scores = empty(n_permutations)

    for i in range(n_permutations):

        shuffle(gene_scores)

        enrichment_scores[i] = compute_enrichment_score(
            gene_scores, gene_set_genes, power=power, statistic=statistic)

    return enrichment_scores
