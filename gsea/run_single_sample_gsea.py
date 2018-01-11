from pandas import DataFrame

from .compute_enrichment_score import compute_enrichment_score
from .support.support.path import establish_path


def run_single_sample_gsea(gene_x_sample,
                           gene_sets,
                           power=1,
                           statistic='ks',
                           plot=False,
                           file_path=None):
    """
    Gene-x-Sample ==> Gene-Set-x-Sample.
    Arguments:
        gene_x_sample (DataFrame): (n_genes, n_samples); avoid having any
            negative score
        gene_sets (DataFrame): (n_gene_sets, max_gene_set_size)
        power (number): power to raise gene_scores
        statistic (str): 'ks' (Kolmogorov-Smirnov) | 'auc' (area under curve)
        plot (bool): whether to plot the mountain plot
        file_path (str):
    Returns:
        DataFrame: (n_gene_sets, n_samples)
    """

    score__gene_set_x_sample = DataFrame(
        index=gene_sets.index, columns=gene_x_sample.columns, dtype=float)
    score__gene_set_x_sample.index.name = 'Gene Set'

    for i, (sample, gene_scores) in enumerate(gene_x_sample.items()):
        print('\n({}/{}) Computing gene set enrichment in {} ...'.format(
            i + 1, gene_x_sample.shape[1], sample))

        min_score = gene_scores.min()
        if min_score < 0:
            print(
                '\tThe sample column has negative value, so adding the minimum of the column ({:.3f}) to it ...'.
                format(min_score))
            gene_scores += min_score

        for i, (gene_set, gene_set_genes) in enumerate(gene_sets.iterrows()):
            print('\t({}/{}) Computing the enrichment of {} ...'.format(
                i + 1, gene_sets.shape[0], gene_set))

            gene_set_genes = gene_set_genes.dropna()

            score = compute_enrichment_score(
                gene_scores,
                gene_set_genes,
                power=power,
                statistic=statistic,
                plot=plot)

            score__gene_set_x_sample.loc[gene_set, sample] = score

    if file_path:
        establish_path(file_path, 'file')
        score__gene_set_x_sample.to_csv(file_path, sep='\t')

    return score__gene_set_x_sample
