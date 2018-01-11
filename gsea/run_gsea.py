from pandas import Series

from .compute_enrichment_score import compute_enrichment_score
from .nd_array.nd_array.compute_log_ratios import compute_log_ratios


def run_gsea(gene_x_sample, phenotypes, gene_sets):
    """
    Run GSEA with 2 phenotypes.
    Arguments:
        gene_x_sample (DataFrame): (n_gene, n_sample)
        phenotypes (iterable):
        gene_sets (DataFrame): (n_gene_set, max_gene_set_size)
    Returns:
        None
    """

    unique_phenotypes = sorted(set(phenotypes))

    if len(unique_phenotypes) != 2:
        raise ValueError('There should be 2 phenotypes.')

    phenotype_means = []

    for phenotype in unique_phenotypes:

        gene_x_phenotype_sample = gene_x_sample.loc[:, [
            phenotype_ == phenotype for phenotype_ in phenotypes
        ]]

        phenotype_mean = gene_x_phenotype_sample.mean(axis=1)
        phenotype_mean.name = phenotype

        phenotype_means.append(phenotype_mean)

    phenotype_mean_0, phenotype_mean_1 = phenotype_means

    gene_scores = Series(
        compute_log_ratios(phenotype_mean_0, phenotype_mean_1),
        name='{} vs {}'.format(*unique_phenotypes),
        index=gene_x_sample.index)

    for gene_set, gene_set_genes in gene_sets.iterrows():

        compute_enrichment_score(
            gene_scores, gene_set_genes.dropna(), plot=True, title=gene_set)
