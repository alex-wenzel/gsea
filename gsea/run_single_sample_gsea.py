from pandas import DataFrame

from .compute_enrichment_score import compute_enrichment_score
from .support.support.path import establish_path


def run_single_sample_gsea(gene_x_sample,
                           gene_sets,
                           normalization_method='rank',
                           power=1,
                           statistic='ks',
                           file_path=None):
    """
    Gene-x-Sample ==> Gene-Set-x-Sample.
    Arguments:
        gene_x_sample (DataFrame): (n_gene, n_sample)
        gene_sets (DataFrame): (n_gene_set, max_gene_set_size)
        normalization_method (str): '-0-' | '0-1' | 'rank'
        power (number): power to raise gene_scores (Gene-x-Sample column)
        statistic (str): 'ks' (Kolmogorov-Smirnov) | 'auc' (area under curve)
        file_path (str): file path to save Gene-Set-x-Sample as .tsv file
    Returns:
        DataFrame: (n_gene_set, n_sample)
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

            score__gene_set_x_sample.loc[
                gene_set, sample] = compute_enrichment_score(
                    gene_scores,
                    gene_set_genes.dropna(),
                    normalization_method=normalization_method,
                    power=power,
                    statistic=statistic)

    if file_path:
        establish_path(file_path, 'file')
        score__gene_set_x_sample.to_csv(file_path, sep='\t')

    return score__gene_set_x_sample
