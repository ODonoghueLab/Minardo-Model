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
resWithOnlySignif <- plotZP(glmTukeyForEachClus.summary, 0.05)



glmTukeyForEachClus <- calcClusterChng(ge.stand, clustered)
glmTukeyForEachClus.summary <- summaryGetZP(glmTukeyForEachClus, ge.stand)
resWithOnlySignif <- plotZP(glmTukeyForEachClus.summary, 0.05)

```
#### 4. Determine events based on 3. (i.e. use z-scores to set threshold for calculating events.)

```R
timeRegions <- getTimeRegionsWithMaximalChange(glmTukeyForEachClus, 9, phosZscoreTh=15, dephosZscoreTh=-15 )
mat_fiftyPoints <- calc50crossing_v3(timeRegions, clustered) # for centroid and plots
  # rearrange mat_fiftyPoints (else causes issues in orderTheEvents())
plotClusters_fifty_v2(humphrey.stand, clustered, mat_fiftyPoints)
plotZP_fifty(glmTukeyForEachClus.summary, mat_fiftyPoints, 0.05)




timeRegions <- getTimeRegionsWithMaximalChange(glmTukeyForEachClus, 7,phosZscoreTh=15, dephosZscoreTh=-15 )
mat_fiftyPoints <- calc50crossing_v3(timeRegions, clustered)

plotClusters_fifty_v2(ge.stand, clustered, mat_fiftyPoints)
plotZP_fifty(glmTukeyForEachClus.summary, mat_fiftyPoints, 0.05)
```



#### 5. Order filtered(/unfiltered events).
Order events (using mean or median) & plot.

```R
mat_fiftyPoints <- mat_fiftyPoints[order(mat_fiftyPoints[,1],mat_fiftyPoints[,2]),]


orderTheEvents(humphrey.stand, clustered, mat_fiftyPoints)

orderTheEvents(ge.stand, clustered, mat_fiftyPoints)
```

#### 6. Publication ready:
	Do steps 1-5 in the background (and then generate figures rearranged, for output).


=========================================================
Get back heat maps for each cluster.
