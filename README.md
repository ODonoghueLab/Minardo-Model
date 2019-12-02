# MinardoModel

An R package for ordering of time-series clusters based on events.

This package provides three major analysis techniques for clusters:
1. Statistically evaluate change occurring in the time intervals of each clusters.
2. Strategy for defining events (eg. phosphorylation and dephosphorylation).  
3. Order and visualise clusters based on the order of occurrence of their first events.


The methods presented here can be applied to a wide variety of time-series high-throughput molecular biology datasets. We demonstrate application to a high-throughput time series [phosphoproteomics data set](./phopho.md), and a time series [gene expression data set](./ge.md).



#### Prerequisite - Mfuzz

The MinardoModel package builds on clustered time-profiles. High throughput time series datasets have been largely clustered using the FCM algorithm, which is implemented in R in the Mfuzz package. Mfuzz is available through Bioconductor. If you don't have it installed, follow the instructions below, or follow the official website ([link](https://doi.org/doi:10.18129/B9.bioc.Mfuzz)).

```R
# Installing Mfuzz
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("Mfuzz")
```




#### API functions are available here.
