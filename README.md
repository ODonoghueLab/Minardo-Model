# MinardoModel

An R package for ordering of time-series clusters based on events.

This package provides three major analysis techniques for clusters:
1. Statistically evaluate change occurring in the time intervals within clusters.
2. Strategy for defining events (eg. phosphorylation and dephosphorylation) based on 50% abundance.  
3. Order and layout clusters based on the first occurrence of events.

Furthermore, this package also presents an alternative colour scheme of cluster plot, which is based on a constant low opacity color.


The methods presented here can be applied to a wide variety of time-series high-throughput molecular biology datasets. In the rest of this document, application to a time-series phosphoproteomics dataset is presented. For application to a gene expression dataset, follow [this link](workflowGE.md).





#### Prerequisite - Mfuzz

The MinardoModel package builds on clustered time-profiles. Time series phosphoproteomics datasets have been largely clustered using the FCM algorithm, which is implemented in R in the Mfuzz package. Mfuzz is available through Bioconductor. If you don't have it installed, follow the instructions below, or follow the official website ([link](https://doi.org/doi:10.18129/B9.bioc.Mfuzz)).
```R
# Installing Mfuzz
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("Mfuzz")
```


## Example workflow

### 1. Load the dataset and standardise it.
An example dataset taken from Humphrey *et al.* [1], contains phospho-proteomics profiles measured over 9 time-points (including basal).

```R

# Load the data
data(humphrey_noDup)
```

The loaded data contains 3,172 profiles. These profiles have been filtered from the originally published dataset [1] (which consisted of 37,248 profiles), as follows:

1. Phosphorylation changes were quantified at all measured time  points (i.e. values are present at all measured time points).
2. Sites were differentially altered (determined by performing empirical Bayes modelling and moderated t-tests, followed by FDR correction).
3. At at-least one quantified time-point there is a 2-fold change.
4. Profiles with duplicate identifiers were removed.

#### Standardisation

After loading the data, the next step is standardising the abundance ratios, as follows:

```R

# Standardise
tmp <- sweep(humphrey.noDup, 1, apply(humphrey.noDup, 1, mean), FUN="-")
humphrey.stand <- sweep(tmp, 1, apply(humphrey.noDup, 1, sd), FUN="/")
humphrey.stand <- as.matrix(humphrey.stand)
remove(tmp) # remove temporary variable
```

### 2. Generate clusters using Mfuzz
These standardised data can then be clustered using the Mfuzz package, as follows:

```R
# Load the Mfuzz library for clustering
library(Mfuzz)

# Do a clustering of the data (and specify the number of clusters using the 'centers' parameter below)
clustered <- cmeans(humphrey.stand, centers=17,  iter.max=200, m=1.25)

# The results can be plotted and visualised as follows
plotClusters(humphrey.stand, clustered)
```


![Mfuzz clustering](images/Humphrey/humphrey_clusters.png)
Fig. 1: Cluster plots of 3,172 profiles from Humphrey *et al.*, generated using the `plotClusters` function.





Alternatively the data can be converted into an expression set and plotted using Mfuzz's plot function:
```
# Note: the results are not shown here.

hum.stand.eset <- new("ExpressionSet", exprs=humphrey.stand)
mfuzz.plot2(hum.stand.eset, cl=clustered, mfrow=c(4,5), centre=TRUE)
```

### 3. Evaluate change in the clusters

This sections depicts how to apply generalised linear models (GLM) to your clusters in order to evaluate change at various time intervals.

```R
# Create and run, for each cluster, a generalised linear model and carry out tukey post-hoc evaluations.
glmTukeyForEachClus <- calcClusterChng(humphrey.stand, clustered)
# This step can take a few minutes.


# Extract z-scores and p-values (each cluster vs time intervals).
glmTukeyForEachClus.summary <- summaryGetZP(glmTukeyForEachClus, totalTimePoints=9)

# Plot the z-scores as a heat map (as seen in Fig. 2).
resWithOnlySignif <- plotZP(glmTukeyForEachClus.summary)

```

![Heat map](images/Humphrey/humphrey_heatmap.png)
Fig. 2: Heat map showing z-scores (differences in mean) for each of the clusters (x-axis) at time-intervals (y-axis) with significant p-values. Z-scores at non-significant intervals are in gray.

The `plotZP` function also returns a matrix, ``resWithOnlySignif``, which consists of z-scores where p-value < 0.5 (or less than the p-value supplied as input to the function).



### 4. Determine events

Phosphorylation and dephosphorylation events are defined when 50% abundance is crossed in the increasing or decreasing direction, respectively.

```R
# Calculate the 50% crossings

mat_fiftyPoints <- calc50crossing(clustered)

```
The `calc50crossing` returns a matrix (`mat_fiftyPoints`) with four columns, where each row contains information about a phosphorylation or dephosphorylation event for the cluster centroid. The columns are:
* Cluster - the cluster containing the event
* Time - the at which the event occurs (i.e when 50% abundance is crossed)
* Abundance - 50% abundance
* Dir - the event type: 1 for phosphorylation and -1 for dephosphorylation.

In our example 24 events are present which can be visualised on the cluster plot as:

```R
plotClusters_fifty(humphrey.stand, clustered, mat_fiftyPoints)

```
![Clusters](images/Humphrey/humphrey_clusters_50.png)
Fig. 3: Cluster plots with phosphorylation and dephosphorylation events, shown via red and blue dots, respectively. The black dashed horizonal line indiates the 50% abundance of the centroid (the bold back line at the center of a cluster).



To plot the events on the z-score heat map:

```R

plotZP_fifty(glmTukeyForEachClus.summary, mat_fiftyPoints)


```

![Heat map](images/Humphrey/humphrey_heatmap_50.png)
Fig. 4: The time-interval change indicating heat map (see Fig. 2), which now also depicts phosphorylation and dephosphorylation events, via red and blue lines, respectively.





### 5. Ordered events

In the following function, an ordering of clusters is calculated, based on events (i.e. phosphorylation and dephosphorylation). Briefly, the method is that, first a time distribution is generated for each cluster. These distributions are compared, either parametrically (t-tests) or non-parametrically (wilcox-test) to determine if the time at which these events occur are different, and if so, which occurs earlier. This gives us an ordering of the events which is visualised as follows.


```R
# Non-parametric test based ordering
mat_fiftyPts_withOrder <- orderTheEvents(humphrey.stand, clustered, mat_fiftyPoints, test="wilcox")
```
![Clusters](images/Humphrey/humphrey_nonParam.png)
Fig. 5: Clusters ordered by first event (where the event ordering was calculated non-parametrically). The events (depicted by dots) which are connected via red dashed lines do not occur at significantly different times.

```R
# Parametric test based ordering
mat_fiftyPts_withOrder <- orderTheEvents(humphrey.stand, clustered, mat_fiftyPoints, test="t-test")

```

![Clusters](images/Humphrey/humphrey_param.png)
Fig. 6: Clusters ordered by first event (where the event ordering was calculated parametrically). Similarly to Fig. 5, those events which are connected via red dashed lines do not occur at significantly different times.

The function ``orderTheEvents`` returns the ordering of the various events in the clusters appended to the matrix `mat_fiftyPoints`.



### References

1. Humphrey SJ, Yang G, Yang P, Fazakerley DJ, Stockli J, Yang JY, James DE. Dynamic adipocyte phosphoproteome reveals that Akt directly regulates mTORC2. Cell metabolism. 2013 Jun 4;17(6):1009-20.
2. O'Donoghue SI, Baldi BF, Clark SJ, Darling AE, Hogan JM, Kaur S, Maier-Hein L, McCarthy DJ, Moore WJ, Stenau E, Swedlow JR. Visualization of Biomedical Data. Annual Review of Biomedical Data Science. 2018 Jul 20;1:275-304.
