# Demonstration on a phosphoproteomics data set

### Data set

We utilise the data set made available by Humphrey *et al.* [1]. This data set measured the global cellular phosphoproteomics response to insulin stimulation. In total, 37,248 phosphosites were quantified, of which 3,172 were significantly changed in response to insulin. These profiles are provided in the dataset `humphrey`.


#### 1. Load the data set and standardise it

The first step is to load and standardize these data.

```R

# Loading the data
data(humphrey)

# Standardising it
tmp <- sweep(humphrey, 1, apply(humphrey, 1, mean), FUN="-")
humphrey.stand <- sweep(tmp, 1, apply(humphrey, 1, sd), FUN="/")
humphrey.stand <- as.matrix(humphrey.stand)
remove(tmp) # remove temporary variable
```

#### 2. Generate clusters
The next step is to cluster these data. For this data set, Yang et al. [2] showed that 17 clusters optimally partitioned these data.

The Mfuzz package provides a fuzzy c-means (FCM) algorithm to partition these data.

```R

library(Mfuzz)

# generate the clusters
clustered <- cmeans(humphrey.stand, centers=17, iter.max=200, m=1.25)

# plot the clusters
plotClusters(humphrey.stand, clustered) ## See Fig. 1
```

![Mfuzz clustering](images/Ge/ge_clusters.png)
Fig. 1: Clusters of gene expression dataset.  

#### 3. Quantifying change within clusters

To quantify which intervals are significantly changing and in which direction, regression can be made use of. The following functions utilize the `glm` function in the `stats` package, and model each cluster as a general linear model, followed by post-hoc Tukey contrasts for every interval.  



```R
# For each cluster, a linear model is formulated (where standardized ratio is the response) and time points along with profiles are predictors; the results of the post hoc tukey contrasting timepoints are presented.
glmTukeyForEachClus <- calcClusterChng(humphrey.stand, clustered)

# Summarizes and returns two matrices containing post hoc tukey z-scores and p-values for consecutive time intervals.
glmTukeyForEachClus.summary <- summaryGetZP(glmTukeyForEachClus, humphrey.stand)

# Returns heatmap plot and a matrix. The matrix contains z-scores for consecutive time points, where values above the significanceTh are set to NA
resWithOnlySignif <- plotZP(glmTukeyForEachClus.summary, significanceTh=0.001) ## See Fig. 2
```

#### 4. Determine events

To calculate events, the first step involves calculating time regions. Time regions are calculated from Tukey contrasts calculated in the previous step. They are non-overlapping intervals with the maximal change, such that change in subintervals must also be in the the direction of the maximal change. The time regions can be filtered by either the z-score or the p-value or both.

It is within these time-regions that events are defined. The threshold can be set to any value between 0 (start of interval) and 1 (end of interval) and reflects on the users intuition or understanding of the system, such as what percent saturation must be attained for the involved phosphosites to contribution to the underlying mechanisms. For example, 30% phosphorylation of a site may be considered enough for the site's corresponding protein to affect downstream proteins.

In the example below however, 50% threshold is set for both a phosphorylation and dephosphorylation event.

```R
# Returns a list of matrices containing the computed time regions for each cluster.
timeRegions <- getTimeRegionsWithMaximalChange(glmTukeyForEachClus, 9, 0.05, phosZscoreTh=15, dephosZscoreTh=-15)


# Returns a matrix containing event information, computed for the cluster centroids.  
mat_events <- calcEvents(timeRegions, clustered)


# Running the function below is optional, its results may assist in updating thresholds to exclude/include events; It returns a matrix containing the number and percentage of profiles which get removed from an event's distribution because these profiles may not follow the general 'differences in means' direction.
mat_missingStats <- missingStats(humphrey.stand, clustered, mat_events, 1, 1)

 # Plot clusters with events overlaid (the events are computed for the cluster centroids).
plotClusters_withEvents(humphrey.stand, clustered, mat_events) ## See Fig. 3.

# Produces a heatmap with events marked, and returns a matrix (similar to the `plotZP` function).
resWithOnlySignif <- plotZP_withEvents(glmTukeyForEachClus.summary, mat_events, 0.001) ## See Fig. 4.
```



#### 5. Order filtered(/unfiltered events).
Order events (using mean or median) & plot.

```R
# Returned is an object containing information regarding the event and cluster order.
theOrder <- calculateOrder(humphrey.stand, clustered, mat_events, "wilcox")

#The order is can then be plotted using event map and event sparkline.
visualizeOrder(theOrder) ## See Fig. 5.

```

#### 6. Rearrange clusters (requires ordering).

Once the order is generated, the clusters can then be rearranged. For this purpose we make available the function `rearrangeCluster`. This function mainly rearranges only two attributes of the `fclust` `clustered`
object, namely the `centers` and the `cluster`.



```R

# Returns the clustered object but with centers and cluster rearranged.
rearranged <- rearrangeClusters(clustered, theOrder)
```

Once rearranged,

```R
# Then, simply recompute everything with the new ordering, and clusters with the rearranged ordering can be visualized.
glmTukeyForEachClus <- calcClusterChng(humphrey.stand, rearranged)
glmTukeyForEachClus.summary <- summaryGetZP(glmTukeyForEachClus, humphrey.stand)

timeRegions <- getTimeRegionsWithMaximalChange(glmTukeyForEachClus, 9, phosZscoreTh=15, dephosZscoreTh=-15)
mat_events <- calcEvents(timeRegions, rearranged)


theOrder <- calculateOrder(humphrey.stand, rearranged, mat_events, "wilcox")
visualizeOrder(theOrder$mat_events_withOrder, theOrder$individEventOrder, theOrder$signifs, theOrder$test) ## See Fig. 6

# Additionally, cluster plots and heatmap plots can also be generated.
```


References

1.
2.
