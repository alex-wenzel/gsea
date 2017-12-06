# GSEA

Library for gene-set enrichment analysis :zap:

## Install

```sh
pip install gsea
```

## Update

```
pip install gsea --upgrade
```

## Get started

```python
from gsea.run_single_sample_gsea import run_single_sample_gsea

gene_set_x_sample = run_single_sample_gsea(gene_x_sample, gene_sets)
```

Look at `notebook/run_single_sample_gsea.ipynb`.

## Development

If you find a bug or have any trouble, please [submit an issue](https://github.com/KwatME/gsea/issues) and I'll take care of it.
