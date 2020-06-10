# MinardoModel

An R package for ordering of events identified from high-throughput time-series biomolecular data sets.

This package provides three major analysis techniques for clusters:
1. Statistically evaluate change occurring in the time intervals of each cluster.
2. Strategy for defining events (eg. phosphorylation and dephosphorylation).  
3. Temporally order events, and visualize this ordering.

The methods presented here can be applied to a wide variety of time-series high-throughput molecular biology datasets. We demonstrate application to a high-throughput time series [phosphoproteomics data set](./phospho.md), a time series [gene expression data set](./ge.md) and a combined [multiomics data set ordering](./multiomics.md).

### Overview (inputs, outputs)

The input required is

#### Prerequisite - Mfuzz (for clustering)

The MinardoModel package builds on clustered time-profiles. If you have not clustered your time profiles, Mfuzz can be used. It implements the commonly utilized k-means and c-means algorithms. Mfuzz is available through Bioconductor. If you don't have it installed, follow the instructions below, or follow the official website ([link](https://doi.org/doi:10.18129/B9.bioc.Mfuzz)).

```R
# Installing Mfuzz
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("Mfuzz")
```


#### Contact us

Any issues or additional suggestions are welcome. Please use the GitHub issue tracker or contact us at sandeep.kaur@unsw.edu.au.
