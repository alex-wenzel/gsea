from pandas import Series

from .nd_array.nd_array.compute_log_ratios import compute_log_ratios


def compute_gene_scores(gene_x_sample, phenotypes, method):
    """
    Compute gene scores.
    Arguments:
        gene_x_sample (DataFrame): (n_gene, n_sample)
        phenotypes (iterable):
        method (str): 'log_ratio' | 'signal_to_noise_ratio' | 't_test' | 'ic'
    Returns:
        Series (n_gene)
    """

    unique_phenotypes = sorted(set(phenotypes))

    if len(unique_phenotypes) != 2:
        raise ValueError('There should be 2 phenotypes.')

    if method == 'log_ratio':

        phenotype_means = []

        for phenotype in unique_phenotypes:

            gene_x_phenotype_sample = gene_x_sample.loc[:, [
                phenotype_ == phenotype for phenotype_ in phenotypes
            ]]

            phenotype_mean = gene_x_phenotype_sample.mean(axis=1)
            phenotype_mean.name = phenotype

            phenotype_means.append(phenotype_mean)

        gene_scores = compute_log_ratios(*phenotype_means)

    elif method == 'signal_to_noise_ratio':
        raise NotImplementedError

    elif method == 't_test':
        raise NotImplementedError

    elif method == 'ic':
        raise NotImplementedError

    else:
        raise ValueError('Unknown method: {}.'.format(method))

    return Series(
        gene_scores,
        name='{} vs {}'.format(*unique_phenotypes),
        index=gene_x_sample.index)
