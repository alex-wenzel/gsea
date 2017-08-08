from numpy import abs, asarray, empty, in1d, where
from numpy.random import shuffle
from pandas import DataFrame

from .dataplay.dataplay.a2d import normalize
from .file.file.file import establish_path
from .file.file.gct import write_gct


def single_sample_gsea(gene_x_sample,
                       gene_sets,
                       normalization=None,
                       power=1,
                       statistic='AUC',
                       file_path=None):
    """
    Gene-x-Sample ==> Gene-Set-x-Sample.
    :param gene_x_sample: DataFrame; (n_genes, n_samples)
    :param gene_sets: DataFrame; (n_gene_sets, max_gene_set_size)
    :param normalization: None | 'rank'
    :param power: number; power to raise gene_scores
    :param statistic: str; 'AUC' (Area Under Curve) | 'KS' (Kolmogorov-Smirnov)
    :param file_path: str;
    :return: DataFrame; (n_gene_sets, n_samples)
    """

    # Rank normalize sample columns
    if normalization == 'raw':
        gene_x_sample = gene_x_sample.copy()

    elif normalization == 'rank':
        # TODO: Change rank method to 'dense'
        gene_x_sample = DataFrame(
            normalize(gene_x_sample, 'rank', axis=0, rank_method='average'),
            index=gene_x_sample.index,
            columns=gene_x_sample.columns)

    else:
        raise ValueError('Unknown statistic: {}.'.format(statistic))

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
        write_gct(gene_set_x_sample, file_path)

    return gene_set_x_sample


def permute_and_compute_enrichment_score(gene_scores,
                                         gene_set_genes,
                                         n_permutations,
                                         power=1,
                                         statistic='AUC'):
    """
    Compute how much permuted gene_scores enriches gene_set_genes.
    :param gene_scores: Series; (n_genes_with_score); sorted and gene indexed
    :param gene_set_genes: iterable; (n_genes)
    :param power: number; power to raise gene_scores
    :param statistic: str; 'AUC' (Area Under Curve) | 'KS' (Kolmogorov-Smirnov)
    :return: array; (n_permutations)
    """

    enrichment_scores = empty(n_permutations)

    for i in range(n_permutations):

        # TODO: Speed up
        # Permute
        shuffle(gene_scores)

        # Compute ES
        enrichment_scores[i] = compute_enrichment_score(
            gene_scores, gene_set_genes, power=power, statistic=statistic)

    return enrichment_scores


def compute_enrichment_score(gene_scores,
                             gene_set_genes,
                             power=1,
                             statistic='AUC'):
    """
    Compute how much gene_scores enriches gene_set_genes.
    :param gene_scores: Series; (n_genes_with_score); sorted and gene indexed
    :param gene_set_genes: iterable; (n_genes)
    :param power: number; power to raise gene_scores
    :param statistic: str; 'AUC' (Area Under Curve) | 'KS' (Kolmogorov-Smirnov)
    :return: float; enrichment score
    """
    gene_scores = gene_scores.sort_values(ascending=False)

    # Check if gene_scores genes are in gene_set_genes
    in_ = in1d(gene_scores.index, gene_set_genes, assume_unique=True)
    in_int = in_.astype(int)

    if power != 1:
        gene_scores = abs(asarray(gene_scores))**power

    hit = (gene_scores * in_int) / gene_scores[in_].sum()
    miss = (1 - in_int) / (in_.size - in_.sum())
    y = hit - miss

    # Compute enrichment score
    cumulative_sums = y.cumsum()

    if statistic == 'AUC':
        enrichment_score = cumulative_sums.sum()

    elif statistic == 'KS':
        max_ = cumulative_sums.max()
        min_ = cumulative_sums.min()
        enrichment_score = where(abs(min_) < abs(max_), max_, min_)

    else:
        raise ValueError('Unknown statistic: {}.'.format(statistic))

    # # TODO: Plot
    # import matplotlib as mpl
    # mpl.pyplot.figure(figsize=(8, 5))
    # ax = mpl.pyplot.gca()
    # ax.plot(range(in_.size), in_, color='#808080', alpha=0.16)
    # ax.plot(range(in_.size), y, color='#9017E6')
    # ax.plot(range(in_.size), cumulative_sums, color='#20D9BA')
    # mpl.pyplot.show()

    return enrichment_score
