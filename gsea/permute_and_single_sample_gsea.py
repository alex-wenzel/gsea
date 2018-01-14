from numpy import empty
from numpy.random import shuffle

from .single_sample_gsea import single_sample_gsea


def permute_and_single_sample_gsea(gene_score,
                                   gene_set_genes,
                                   n_permutation,
                                   normalization_method='rank',
                                   power=1,
                                   statistic='ks'):
    """
    Permute and single-sample GSEA.
    Arguments:
        gene_score (Series): (n_gene_with_score)
        gene_set_genes (iterable): (n_gene)
        n_permutation (int):
        normalization_method (str): '-0-' | '0-1' | 'rank'
        power (number): power to raise gene_score
        statistic (str): 'ks' (Kolmogorov-Smirnov) | 'auc' (area under curve)
    Returns:
        array: (n_permutation)
    """

    permutation_score = empty(n_permutation)

    for i in range(n_permutation):

        shuffle(gene_score)

        permutation_score[i] = single_sample_gsea(
            gene_score,
            gene_set_genes,
            normalization_method=normalization_method,
            power=power,
            statistic=statistic)

    return permutation_score
