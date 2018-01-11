from numpy import empty
from pandas import Series

from .compute_gene_scores import compute_gene_scores
from .singe_sample_gsea import singe_sample_gsea


def gsea(gene_x_sample,
         phenotypes,
         gene_sets,
         method='log_ratio',
         normalization_method='rank',
         power=1,
         statistic='ks',
         n_permutation=0,
         permuting='phenotype'):
    """
    Run GSEA with 2 phenotypes.
    Arguments:
        gene_x_sample (DataFrame): (n_gene, n_sample)
        phenotypes (iterable):
        gene_sets (DataFrame): (n_gene_set, max_gene_set_size)
        method (str): 'log_ratio' | 'signal_to_noise_ratio' | 't_test' | 'ic'
        normalization_method (str): '-0-' | '0-1' | 'rank'
        power (number): power to raise gene_scores (Gene-x-Sample column)
        statistic (str): 'ks' (Kolmogorov-Smirnov) | 'auc' (area under curve)
        n_permutation (int):
        permuting (str): 'phenotype' | 'gene'
    Returns:
        Series: (n_gene_set); enrichment socres
        Series: (n_gene_set); p-values
    """

    gene_set_score = Series(
        name='Enrichment Score', index=gene_sets.index, dtype=float)
    gene_set_score.index.name = 'Gene Set'

    gene_set_p_value = Series(
        name='P-Value', index=gene_sets.index, dtype=float)
    gene_set_p_value.index.name = 'Gene Set'

    gene_scores = compute_gene_scores(gene_x_sample, phenotypes, method)

    for gene_set, gene_set_genes in gene_sets.iterrows():

        enrichment_score = singe_sample_gsea(
            gene_scores,
            gene_set_genes,
            normalization_method=normalization_method,
            power=power,
            statistic=statistic,
            plot=True,
            title=gene_set)

        if n_permutation:

            permutation_enrichment_score = empty(n_permutation)

            for i in range(n_permutation):

                if permuting == 'phenotype':
                    raise NotImplementedError

                elif permuting == 'gene':
                    raise NotImplementedError

                else:
                    raise ValueError(
                        'Unknown permuting: {}.'.format(permuting))

                gene_scores = compute_gene_scores(gene_x_sample, phenotypes,
                                                  method)

                permutation_enrichment_score[i] = singe_sample_gsea(
                    gene_scores, gene_set_genes)

        gene_set_score[gene_set] = enrichment_score
        gene_set_p_value[gene_set] = sum(
            enrichment_score < permutation_enrichment_score) / n_permutation

    return gene_set_score
