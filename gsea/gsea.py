from random import sample, shuffle

from numpy import empty
from pandas import Series

from .compute_gene_scores import compute_gene_scores
from .single_sample_gsea import single_sample_gsea


def gsea(gene_x_sample,
         phenotypes,
         gene_sets,
         method='log_ratio',
         normalization_method='rank',
         power=1,
         statistic='ks',
         n_permutation=0,
         permuting='phenotype',
         directory_path=None):
    """
    GSEA (with 2 phenotypes).
    Arguments:
        gene_x_sample (DataFrame): (n_gene, n_sample)
        phenotypes (iterable):
        gene_sets (DataFrame): (n_gene_set, max_gene_set_size)
        method (str): 'log_ratio' | 'signal_to_noise_ratio' | 't_test' | 'ic'
        normalization_method (str): '-0-' | '0-1' | 'rank'
        power (number): power to raise gene_score (Gene-x-Sample column)
        statistic (str): 'ks' (Kolmogorov-Smirnov) | 'auc' (area under curve)
        n_permutation (int):
        permuting (str): 'phenotype' | 'gene'
        directory_path (str):
    Returns:
        Series: (n_gene_set); scores
        Series: (n_gene_set); p-values
    """

    gene_set_score = Series(
        name='Enrichment Score', index=gene_sets.index, dtype=float)
    gene_set_score.index.name = 'Gene Set'

    gene_set_p_value = Series(
        name='P-Value', index=gene_sets.index, dtype=float)
    gene_set_p_value.index.name = 'Gene Set'

    gene_score = compute_gene_scores(gene_x_sample, phenotypes, method)

    for gene_set, gene_set_genes in gene_sets.iterrows():

        score = single_sample_gsea(
            gene_score,
            gene_set_genes,
            normalization_method=normalization_method,
            power=power,
            statistic=statistic,
            plot=True,
            title=gene_set)

        gene_set_score[gene_set] = score

        if n_permutation:

            permutation_score = empty(n_permutation)

            permuting__gene_x_sample = gene_x_sample.copy()
            permuting__phenotypes = list(phenotypes)

            for i in range(n_permutation):

                if permuting == 'phenotype':
                    shuffle(permuting__phenotypes)

                elif permuting == 'gene':
                    permuting__gene_x_sample.index = sample(
                        permuting__gene_x_sample.index.tolist(),
                        permuting__gene_x_sample.shape[0])
                else:
                    raise ValueError(
                        'Unknown permuting: {}.'.format(permuting))

                permutation_score[i] = single_sample_gsea(
                    compute_gene_scores(permuting__gene_x_sample,
                                        permuting__phenotypes, method),
                    gene_set_genes,
                    normalization_method=normalization_method,
                    power=power,
                    statistic=statistic)

            if 0 < score:
                p_value = sum(score <= permutation_score) / n_permutation
            else:
                p_value = sum(permutation_score <= score) / n_permutation
            gene_set_p_value[gene_set]
        else:
            gene_set_p_value[gene_set] = None

    return gene_set_score, gene_set_p_value
