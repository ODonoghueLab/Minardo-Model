# Demonstration on a phosphoproteomics data set

### Data set


### Workflow


#### 1. Load the data set and standardise it
```R

# Loading the data
data(humphrey_noDup)

# Standardising it
tmp <- sweep(humphrey.noDup, 1, apply(humphrey.noDup, 1, mean), FUN="-")
humphrey.stand <- sweep(tmp, 1, apply(humphrey.noDup, 1, sd), FUN="/")
humphrey.stand <- as.matrix(humphrey.stand)
remove(tmp) # remove temporary variable
```

#### 2. Generate clusters

We use the FCM algorithm made available in Mfuzz.
```R

library(Mfuzz)

# generate the clusters
clustered <- cmeans(humphrey.stand, centers=17, iter.max=200, m=1.25)

# plot the clusters
plotClusters(humphrey.stand, clustered)

```

#### 3. Quantifying change within clusters

```R

glmTukeyForEachClus <- calcClusterChng(humphrey.stand, clustered)
glmTukeyForEachClus.summary <- summaryGetZP(glmTukeyForEachClus, humphrey.stand)
resWithOnlySignif <- plotZP(glmTukeyForEachClus.summary)

```
#### 4. Determine events based on 3. (i.e. use z-scores to set threshold for calculating events.)

```R
timeRegions <- getTimeRegionsWithMaximalChange(glmTukeyForEachClus, 9, 0.05)
mat_fiftyPoints <- calc50Crossing_v2(timeRegions, clustered) # for centroid and plots
plotClusters_fifty_v2
plotZP_fifty_v2
```

#### 5. Order filtered(/unfiltered events).
Order events (using mean or median) & plot.


#### 6. Publication ready:
	Do steps 1-5 in the background (and then generate figures rearranged, for output).


=========================================================
Get back heat maps for each cluster.
