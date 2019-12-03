# Demonstration on a phosphoproteomics data set

### Data set

We utilise the data set made available by Humphrey *et al.* [1]. This data set measured the global cellular phosphoproteomics response to insulin stimulation. In total, 37,248 phosphosites were quantified, of which 3,172 were significantly changed in response to insulin. These profiles are provided in the dataset humphrey_noDup.


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
The next step is to cluster these data. For this data set, Yang et al. showed that 17 clusters optimally partitioned these data.

The Mfuzz package provides a fuzzy c-means (FCM) algorithm to partition these data.

```R

library(Mfuzz)

# generate the clusters
clustered <- cmeans(humphrey.stand, centers=17, iter.max=200, m=1.25)

# plot the clusters
plotClusters(humphrey.stand, clustered)
```

#### 3. Quantifying change within clusters

To quantify which intervals are significantly changing and in which direction, regression can be made use of. The following functions utilize the `glm` function in the stats package, and model each cluster as a general linear model, followed by post-hoc Tukey contrasts for every interval.  



```R
# For each cluster, a linear model is formulates (where standardized ratio is the response) and time point and profile are predictions; the results of the post hoc tukey contrasting timepoints is presented.
glmTukeyForEachClus <- calcClusterChng(humphrey.stand, clustered)

# Two matrices containing post hoc tukey z-scores and p-values for consecutive time intervals.
glmTukeyForEachClus.summary <- summaryGetZP(glmTukeyForEachClus, humphrey.stand)

# Heatmap; a matrix is returned containing z-scores for consecutive time points, where values above the significanceTh are set to NA
resWithOnlySignif <- plotZP(glmTukeyForEachClus.summary, significanceTh=0.001)
```

#### 4. Determine events

To calculate events, the first step involves calculating time regions. Time regions are calculated from Tukey contrasts calculated in the previous step. They are non-overlapping intervals with the maximal change, such that change in subintervals must also be in the the direction of the maximal change. The time regions can be filtered by both either the z-score or the p-value or both.

It is within these time-regions that events are defined. The threshold can be set to any value between 0 (start of interval) and 1 (end of interval) and reflects on the users intuition or understanding of the system, that what percent saturation must be attained for the involved phosphosites to contribution to the underlying mechanisms. For example, 30% phosphorylation of a site may be considered enough for the site's corresponding protein to affect downstream proteins.

In the example below however, 50% threshold is set for both a phosphorylation and dephosphorylation event.

```R
timeRegions <- getTimeRegionsWithMaximalChange(glmTukeyForEachClus, 9, phosZscoreTh=15, dephosZscoreTh=-15) # a list of matrices containing the computed time regions for each cluster.

mat_events <- calcEvents(timeRegions, clustered) # a matrix containing the event positions computed for the cluster centroids.

mat_missingStats <- missingStats(humphrey.stand, clustered, mat_events, 1, 1) # just an informational function, which may assist in updating thresholds to exclude/include events; presents the number and percentage of profiles which get removed from a distribution because these profiles may not follow the general 'differences in means' direction.

plotClusters_withEvents(humphrey.stand, clustered, mat_events) # produces the clusters with the time regions, 50% of standardized abundance within a time region, and events marked.

resWithOnlySignif <- plotZP_withEvents(glmTukeyForEachClus.summary, mat_events, 0.001) # produces the heatmap with events marked, and returns a matrix, which is also returned by plotZP.
```



#### 5. Order filtered(/unfiltered events).
Order events (using mean or median) & plot.

```R

theOrder <- calculateOrder(humphrey.stand, clustered, mat_events, "wilcox")
visualizeOrder(theOrder$mat_events_withOrder, theOrder$individEventOrder, theOrder$signifs, theOrder$test)

```

##### 6. Rearrange clusters (requires ordering).


```R
## orig
glmTukeyForEachClus_orig <- calcClusterChng(humphrey.stand, clustered)
glmTukeyForEachClus.summary_orig <- summaryGetZP(glmTukeyForEachClus_orig, humphrey.stand)
timeRegions_orig <- getTimeRegionsWithMaximalChange(glmTukeyForEachClus_orig, 9, phosZscoreTh=15, dephosZscoreTh=-15)
mat_events_orig <- calcEvents(timeRegions_orig, clustered)
theOrder_orig <- calculateOrder(humphrey.stand, clustered, mat_events_orig, "wilcox")
visualizeOrder(theOrder_orig$mat_events_withOrder, theOrder_orig$individEventOrder, theOrder_orig$signifs, theOrder_orig$test)

rearranged <- rearrangeClusters(clustered, theOrder)


# ordered
# then simply recompute everything with the new ordering.
glmTukeyForEachClus <- calcClusterChng(humphrey.stand, rearranged)
glmTukeyForEachClus.summary <- summaryGetZP(glmTukeyForEachClus, humphrey.stand)
timeRegions <- getTimeRegionsWithMaximalChange(glmTukeyForEachClus, 9, phosZscoreTh=15, dephosZscoreTh=-15)
mat_events <- calcEvents(timeRegions, rearranged)


theOrder <- calculateOrder(humphrey.stand, rearranged, mat_events, "wilcox")
visualizeOrder(theOrder$mat_events_withOrder, theOrder$individEventOrder, theOrder$signifs, theOrder$test)
```


Then repeat steps 1-5 to generate figures with clusters rearranged if required.




=========================================================
Get back heat maps for each cluster.
