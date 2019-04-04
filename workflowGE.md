# Demonstration on Gene Expression dataset

In this document, we demonstrate how our package can also be applied to a wider variety of time series datasets, such as gene expression datasets.

### Dataset

We utilise the dataset made available by Ma *et al.* [1]. We downloaded the raw data from GEO (Gene Expression Omnibus) with the identifier GSE40565 and considered only the time-series arrays. We then obtained the differentially expressed genes from the seven time points (including basal) by following the section 'Microarray analysis' in the paper. This resulted in 2566 profiles, which we make available as a sample dataset.  

### Workflow

Overall, the work flow is very similar to the time-series phosphoproteomics datasets.

#### 1. Load the dataset and standardise it

```
# Load the data.
data(ge)


# Standardisation
tmp <- sweep(ge.temporal.changed, 1, apply(ge.temporal.changed, 1, mean), FUN="-")
ge.stand <- sweep(tmp, 1, apply(ge.temporal.changed, 1, sd), FUN="/")
ge.stand <- as.matrix(ge.stand)
remove(tmp)
```

#### 2. Clustering using Mfuzz

```
library(Mfuzz)

clustered <- cmeans(ge.stand, centers=20,  iter.max=200, m=1.25) # paper, 2014_ma_etal, says 20 clusters.

# Plotting the clusters
ge.stand.eset <- new("ExpressionSet", exprs=ge.stand)
mfuzz.plot2(ge.stand.eset, cl=clustered, mfrow=c(4,5), centre=TRUE)
```

#### 3. Quantifying change within clusters

```
# Create and run, for each cluster, a generalised linear model and carry out tukey post-hoc evaluations.
glmTukeyForEachClus <- calcClusterChng(ge.stand, clustered)

# Extract z-scores and p-values.
glmTukeyForEachClus.summary <- summaryGetZP(glmTukeyForEachClus, totalTimePoints=7)

## Plot the z-scores as a heatmap.
resWithOnlySignif <- plotZP(glmTukeyForEachClus.summary)
```



#### 4. Determine events

```
# Calculate 50% crossings (time and direction at the 50% abundance)
mat_fiftyPoints <- calc50crossing(clustered)


# Overlaying the 50% crossings on the heat map.
plotZP_fifty(glmTukeyForEachClus.summary, mat_fiftyPoints, 0.5)

# Overlaying the 50% crossings on the cluster plot.
plotClusters_fifty(ge.stand, clustered, mat_fiftyPoints)
```



#### 5. Plotting events
In the following function, an ordering is calculated (i.e. those occurring at statistically significantly different times) and a figure is generated where the clusters are ordered by occurrence of first event.

```
# Non-parametric test based ordering
orderTheEvents(ge.stand, clustered, mat_fiftyPoints, test="wilcox")

# Parametric test based ordering
orderTheEvents(ge.stand, clustered, mat_fiftyPoints, test="t-test")
```


To save the generated image as pdf (the width and height can be increased):
```
pdf("ge_dups_wilcox.pdf", width = 10, height=8)
orderTheEvents(ge.stand, clustered, mat_fiftyPoints)
dev.off()
```


References

1. Ma, X., Yang, P., Kaplan, W.H., Lee, B.H., Wu, L.E., Yang, J.Y.H., Yasunaga, M., Sato, K., Chisholm, D.J. and James, D.E., 2014. ISL1 Regulates Peroxisome Proliferator-Activated Receptor gamma Activation and Early Adipogenesis via Bone Morphogenetic Protein 4-Dependent and-Independent Mechanisms. *Molecular and cellular biology*, 34(19), pp.3607-3617.
