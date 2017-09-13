from numpy import empty
from numpy.random import shuffle

from .compute_enrichment_score import compute_enrichment_score


def permute_and_compute_enrichment_score(gene_scores,
                                         gene_set_genes,
                                         n_permutations,
                                         power=1,
                                         statistic='AUC'):
    """
    Compute how much permuted gene_scores enriches gene_set_genes.
    Arguments:
        gene_scores (Series): (n_genes_with_score); sorted and gene indexed
        gene_set_genes (iterable): (n_genes)
        power (number): Power to raise gene_scores
        statistic (str): 'AUC' (Area Under Curve) | 'KS' (Kolmogorov-Smirnov)
    Returns:
        array: (n_permutations)
    """

    enrichment_scores = empty(n_permutations)

    for i in range(n_permutations):

        # Permute
        shuffle(gene_scores)

        # Compute ES
        enrichment_scores[i] = compute_enrichment_score(
            gene_scores, gene_set_genes, power=power, statistic=statistic)

    return enrichment_scores
