from pandas import DataFrame

from .array_nd.array_nd.normalize_2d_array import normalize_2d_array
from .compute_enrichment_score import compute_enrichment_score
from .support.support.path import establish_path


def run_single_sample_gsea(gene_x_sample,
                           gene_sets,
                           normalization=None,
                           power=1,
                           statistic='auc',
                           file_path=None):
    """
    Gene-x-Sample ==> Gene-Set-x-Sample.
    Arguments:
        gene_x_sample (DataFrame): (n_genes, n_samples)
        gene_sets (DataFrame): (n_gene_sets, max_gene_set_size)
        normalization (str): None | 'rank'
        power (number): power to raise gene_scores
        statistic (str): 'auc' (area under curve) | 'ks' (Kolmogorov-Smirnov)
        file_path (str):
    Returns:
        DataFrame: (n_gene_sets, n_samples)
    """

    # Rank normalize sample columns
    if normalization == 'rank':
        # TODO: change rank method from 'average' to 'dense'
        index = gene_x_sample.index
        columns = gene_x_sample.columns

        gene_x_sample = DataFrame(
            normalize_2d_array(
                gene_x_sample.values, 'rank', axis=0, rank_method='average'),
            index=index,
            columns=columns)
    else:
        gene_x_sample = gene_x_sample.copy()

    # Make Gene-Set-x-Sample place holder
    gene_set_x_sample = DataFrame(
        index=gene_sets.index, columns=gene_x_sample.columns, dtype=float)

    # For each gene set
    for gene_set, gene_set_genes in gene_sets.iterrows():
        print('Computing {} enrichment ...'.format(gene_set))

        gene_set_genes.dropna(inplace=True)

        # Compute enrichment score (ES) for each sample
        for sample, gene_scores in gene_x_sample.items():
            gene_set_x_sample.ix[gene_set, sample] = compute_enrichment_score(
                gene_scores, gene_set_genes, power=power, statistic=statistic)

    if file_path:
        establish_path(file_path)
        gene_set_x_sample.to_csv(file_path, sep='\t')

    return gene_set_x_sample
