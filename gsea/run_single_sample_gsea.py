from pandas import DataFrame

from .compute_enrichment_score import compute_enrichment_score
from .support.support.path import establish_path


def run_single_sample_gsea(gene_x_sample,
                           gene_sets,
                           normalization=None,
                           power=1,
                           statistic='auc',
                           plot=False,
                           file_path=None):
    """
    Gene-x-Sample ==> Gene-Set-x-Sample.
    Arguments:
        gene_x_sample (DataFrame): (n_genes, n_samples)
        gene_sets (DataFrame): (n_gene_sets, max_gene_set_size)
        normalization (str): None | 'rank'
        power (number): power to raise gene_scores
        statistic (str): 'auc' (area under curve) | 'ks' (Kolmogorov-Smirnov)
        plot (bool): whether to plot the mountain plot
        file_path (str):
    Returns:
        DataFrame: (n_gene_sets, n_samples)
    """

    gene_set_x_sample = DataFrame(
        index=gene_sets.index, columns=gene_x_sample.columns, dtype=float)

    for i, (gene_set, gene_set_genes) in enumerate(gene_sets.iterrows()):
        print('({}/{}) Computing {} enrichment ...'.format(
            i + 1, gene_sets.shape[0], gene_set))

        gene_set_genes.dropna(inplace=True)

        for sample, gene_scores in gene_x_sample.items():

            gene_set_x_sample.ix[gene_set, sample] = compute_enrichment_score(
                gene_scores,
                gene_set_genes,
                power=power,
                statistic=statistic,
                plot=plot)

    if file_path:
        establish_path(file_path)
        gene_set_x_sample.to_csv(file_path, sep='\t')

    return gene_set_x_sample
