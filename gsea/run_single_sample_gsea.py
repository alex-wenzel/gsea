from pandas import DataFrame

from .compute_enrichment_score import compute_enrichment_score
from .support.support.path import establish_path


def run_single_sample_gsea(gene_x_sample,
                           gene_sets,
                           power=1,
                           statistic='auc',
                           n_permutation=0,
                           plot=False,
                           file_path=None):
    """
    Gene-x-Sample ==> Gene-Set-x-Sample.
    Arguments:
        gene_x_sample (DataFrame): (n_genes, n_samples)
        gene_sets (DataFrame): (n_gene_sets, max_gene_set_size)
        power (number): power to raise gene_scores
        statistic (str): 'auc' (area under curve) | 'ks' (Kolmogorov-Smirnov)
        n_permutation (int): number of permutations to compute p-value
        plot (bool): whether to plot the mountain plot
        file_path (str):
    Returns:
        DataFrame: (n_gene_sets, n_samples)
    """

    score__gene_set_x_sample = DataFrame(
        index=gene_sets.index, columns=gene_x_sample.columns, dtype=float)

    p_value__gene_set_x_sample = gene_set_x_sample = DataFrame(
        index=gene_sets.index, columns=gene_x_sample.columns, dtype=float)

    if n_permutation:
        column_permuted__gene_x_sample = []

    for i, (gene_set, gene_set_genes) in enumerate(gene_sets.iterrows()):
        print('({}/{}) Computing {} enrichment ...'.format(
            i + 1, gene_sets.shape[0], gene_set))

        gene_set_genes.dropna(inplace=True)

        for sample, gene_scores in gene_x_sample.items():

            score = compute_enrichment_score(
                gene_scores,
                gene_set_genes,
                power=power,
                statistic=statistic,
                plot=plot)

            score__gene_set_x_sample.loc[gene_set, sample] = score

            if n_permutation:
                for i in range(n_permutation):

                    permuted_gene_score = column_permuted__gene_x_sample[i][
                        sample]

                    permuted_scores = compute_enrichment_score(
                        permuted_gene_score,
                        gene_set_genes,
                        power=power,
                        statistic=statistic,
                        plot=plot)

                p_value_higher = (permuted_scores <= score) / n_permutation
                p_value_lower = (score <= permuted_scores) / n_permutation

                p_value__gene_set_x_sample.loc[gene_set, sample] = [
                    p_value_lower, p_value_higher
                ][p_value_higher < p_value_lower]

    if file_path:
        establish_path(file_path)
        gene_set_x_sample.to_csv(file_path, sep='\t')

    return score__gene_set_x_sample, p_value__gene_set_x_sample
