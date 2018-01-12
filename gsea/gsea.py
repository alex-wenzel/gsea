from os.path import join
from random import sample, shuffle

from numpy import empty, empty_like
from pandas import DataFrame

from .compute_gene_scores import compute_gene_scores
from .single_sample_gsea import single_sample_gsea
from .support.support.json_ import write_json
from .support.support.path import clean_name, establish_path


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
        DataFrame: (n_gene_set); columns=('Score', 'P-Value', 'FDR', )
    """

    scores = empty(gene_sets.shape[0])
    p_values = empty_like(scores)
    fdrs = empty_like(scores)

    gene_score = compute_gene_scores(gene_x_sample, phenotypes, method)

    if directory_path:
        establish_path(directory_path, 'directory')

        write_json({
            'method': method,
            'normalization_method': normalization_method,
            'power': power,
            'statistic': statistic,
            'n_permutation': n_permutation,
        }, join(directory_path, 'parameters.json'))

    for i, (gene_set, gene_set_genes) in enumerate(gene_sets.iterrows()):

        if directory_path:
            plot_file_path = join(directory_path,
                                  '{}.png'.format(clean_name(gene_set)))
        else:
            plot_file_path = None

        score = single_sample_gsea(
            gene_score,
            gene_set_genes,
            normalization_method=normalization_method,
            power=power,
            statistic=statistic,
            plot=True,
            title=gene_set,
            plot_file_path=plot_file_path)

        scores[i] = score

        if n_permutation:

            permutation_scores = empty(n_permutation)

            permuting__gene_x_sample = gene_x_sample.copy()
            permuting__phenotypes = list(phenotypes)

            for j in range(n_permutation):

                if permuting == 'phenotype':
                    shuffle(permuting__phenotypes)

                elif permuting == 'gene':
                    permuting__gene_x_sample.index = sample(
                        permuting__gene_x_sample.index.tolist(),
                        permuting__gene_x_sample.shape[0])
                else:
                    raise ValueError(
                        'Unknown permuting: {}.'.format(permuting))

                permutation_scores[j] = single_sample_gsea(
                    compute_gene_scores(permuting__gene_x_sample,
                                        permuting__phenotypes, method),
                    gene_set_genes,
                    normalization_method=normalization_method,
                    power=power,
                    statistic=statistic)

            if 0 < score:
                n_equal_more_extreme = sum(score <= permutation_scores)
            else:
                n_equal_more_extreme = sum(permutation_scores <= score)

            p_values[i] = max(1, n_equal_more_extreme) / n_permutation

        else:
            p_values[i] = None

    gene_set_score_p_value = DataFrame(index=gene_sets.index)
    gene_set_score_p_value['Score'] = scores
    gene_set_score_p_value['P-Value'] = p_values
    gene_set_score_p_value['FDR'] = fdrs

    gene_set_score_p_value.sort_values('Score', inplace=True)

    if directory_path:
        gene_set_score_p_value.to_csv(
            join(directory_path, 'gene_set_score_p_value.tsv'), sep='\t')

    return gene_set_score_p_value
