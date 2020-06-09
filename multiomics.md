# Demonstration on a multiomics data set

This is still underdevelopment.

## Data set

The Yang *et al.* 2019 data set comprises time series measurements of phosphoproteomics, transcriptomics and proteomics. Yang *et al.* clustered phosphoproteomics and proteomics data using FCM and made available these clusters. We clustered the transcriptomics data using the STEM tool, and make the clusters available here.


## The measured time intervals
The following variable indicates the time points at which each of the three data sets were measured (see Yang *et al.* 2019 Fig. ). It will be required for the final combined ordering step.

```R
phosTimeMap = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
mRnaTimeMap = c(1, 5, 7, 8, 9, 10, 11, 12)
protTimeMap = c(1, 4, 5, 7, 8, 9, 10, 11, 12)

```

## 1. Load the multiomics data sets and clusters

```R
# Loading the data
data(multiomics)

```

## 2. Quantifying change within clusters for each subset
```R
# phosphoproteomics
phos_glmTukey <- calcClusterChng(phosTsData.stand, phosClusVec)
phos_glmTukey_summary <- summaryGetZP(phos_glmTukey, phosTsData.stand)


# Transcriptomics
rna_glmTukey <- calcClusterChng(rnaSeqTsData.stand, rnaSeqClusVec)
rna_glmTukey_summary <- summaryGetZP(rna_glmTukey, rnaSeqTsData.stand)

# Proteomics
prot_glmTukey <- calcClusterChng(protTsData.stand, protClusVec)
prot_glmTukey_summary <- summaryGetZP(prot_glmTukey, protTsData.stand)

```


## 3. Determine events for each subset
```R
# Phosphoproteomics
phos_eventWindows <- getTimeRegionsWithMaximalChange(phos_glmTukey, ncol(phosTsData.stand), 0.05, phosZscoreTh=15, dephosZscoreTh=-15)
phos_events <- calcEvents(phos_eventWindows, phosClusVec, phosTsData.stand)


# Transcriptomics
rna_eventWindows <- getTimeRegionsWithMaximalChange(rna_glmTukey, ncol(rnaSeqTsData.stand), 0.05, phosZscoreTh=15, dephosZscoreTh=-15)
rnaSeq_events <- calcEvents(rna_eventWindows, rnaSeqClusVec, rnaSeqTsData.stand)


# Proteomics
prot_eventWindows <- getTimeRegionsWithMaximalChange(prot_glmTukey, ncol(protTsData.stand), 0.05, phosZscoreTh=15, dephosZscoreTh=-15)
prot_events <- calcEvents(prot_eventWindows, protClusVec, protTsData.stand)

```


## 4. Order all events combined
```R
list_tsData = list(phosTsData.stand, rnaSeqTsData.stand, protTsData.stand)
list_clusters = list(phosClusVec, rnaSeqClusVec, protClusVec)
list_events = list(phos_events, rnaSeq_events, prot_events)

list_timeMap = list(phosTimeMap, mRnaTimeMap, protTimeMap)


res <- calculateOrder_combined(list_tsData, list_clusters, list_events, list_timeMap, 't-test')

pdf("combined_order.pdf", width=10, height=10)
visualizeOrder_combined(res)
dev.off()
```
