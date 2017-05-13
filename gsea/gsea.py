from numpy import asarray, empty, in1d, where
from numpy.random import shuffle
from pandas import DataFrame

from .dataplay.dataplay.a2d import normalize
from .file.file.file import establish_path
from .file.file.gct import write_gct


def single_sample_gsea(gene_x_sample,
                       gene_sets,
                       power=1,
                       statistic='Kolmogorov-Smirnov',
                       n_permutations=0,
                       file_path=None):
    """
    Gene-x-Sample ==> Gene-Set-x-Sample.
    :param gene_x_sample: DataFrame; (n_genes, n_samples)
    :param gene_sets: DataFrame;
    :param power: number; power to raise gene_scores
    :param statistic: str; 'Kolmogorov-Smirnov' | 'Cumulative Area'
    :param n_permutations: int;
    :return: DataFrame; (n_gene_sets, n_samples)
    """

    # Rank normalize columns
    # TODO: Check if multiplying by 10000 is needed
    g_x_s = normalize(gene_x_sample, 'rank', axis=0)

    # Make Gene-Set-x-Sample place holder
    gs_x_s = DataFrame(index=gene_sets.index, columns=g_x_s.columns)

    # For each gene set
    for gs_i, gs in gene_sets.iterrows():
        print('Computing {} enrichment ...'.format(gs_i))

        gs.dropna(inplace=True)

        # For each sample
        for s_i, s in g_x_s.items():

            # Compute enrichment score (ES)
            es = compute_enrichment_score(
                s, gs, power=power, statistic=statistic)

            if 0 < n_permutations:  # Score is permutation-normalized ES

                # ESs computed with permuted sample gene scores
                ess = empty(n_permutations)

                for i in range(n_permutations):
                    # Permute
                    shuffle(s)
                    # Compute ES
                    ess[i] = compute_enrichment_score(
                        s, gs, power=power, statistic=statistic)

                # Compute permutation-normalized enrichment score
                gs_x_s.ix[gs_i, s_i] = es / ess.mean()

            else:  # Score is ES
                gs_x_s.ix[gs_i, s_i] = es

    if file_path:
        establish_path(file_path)
        write_gct(gs_x_s, file_path)

    return gs_x_s


def compute_enrichment_score(gene_scores,
                             gene_set_genes,
                             power=1,
                             statistic='Kolmogorov-Smirnov'):
    """
    Compute how much gene_scores enriches gene_set_genes.
    :param gene_scores: Series; (n_genes_with_score); index must be genes
    :param gene_set_genes: iterable; (n_genes)
    :param power: number; power to raise gene_scores
    :param statistic: str; 'Kolmogorov-Smirnov' | 'Cumulative Area'
    :return: float; enrichment score
    """

    # Copy and sort gene_scores
    gss = gene_scores.sort_values(ascending=False)**power

    # Check if gene_scores genes are in gene_set_genes
    in_ = in1d(gss.index, gene_set_genes, assume_unique=True)
    in_int = in_.astype(int)

    # Score: values-at-hits / sum(values-at-hits) - is-miss' / number-of-misses
    gss = asarray(gss)
    y = (gss * in_int / gss[in_].sum()) - (1 - in_int) / (in_.size - in_.sum())

    # Compute enrichment score
    if statistic == 'Kolmogorov-Smirnov':

        cs = y.cumsum()

        max_es = cs.max()
        min_es = cs.min()

        es = where(abs(min_es) < abs(max_es), max_es, min_es)

    elif statistic == 'Cumulative Area':
        pass

    else:
        raise ValueError('Unknown statistic: {}.'.format(statistic))

    # TODO: Plot results
    # mpl.pyplot.figure(figsize=(8, 5))
    #     ax = mpl.pyplot.gca()
    #     ax.plot(range(in_.size), in_, color='black', alpha=0.16)
    #     ax.plot(range(in_.size), s)
    #     ax.plot(range(in_.size), cs)
    #     mpl.pyplot.show()

    return es
