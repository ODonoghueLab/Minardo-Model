#' @title
#' Calculate changes at all intervals within each cluster
#'
#' @description
#' Create and run generalised linear models and post-hoc tukey contrasts at all time points for each cluster.
#'
#' @param Tc A matrix containing time course data.
#' @param clusters A vector with same length as number of profiles. The cluster labels are assumed to be numbers, starting from 1.
#'
#' @return A list (with same length as number of clusters) containing the Tukey summaries.
#'
#'
#' @importFrom multcomp glht mcp
#' @importFrom stats glm
#' @importFrom methods is
#'
#'
#' @seealso \code{\link[e1071]{cmeans}} for clustering.
#'
#' @export
calcClusterChng <- function(Tc, clusters){


	# 1. split into matrix for each cluster...
	# 2. convert to df. of required format for glm
	# 3. for each cluster run glm.
	# 4. for each glm result run tukey contrast. (which can give the summary directly or which can initially give the full matrix (between each time point))
	# another function to do the plotting.

	stopifnot(is(Tc, "matrix"), (length(clusters) == nrow(Tc)))

	list_matrices <- splitIntoSubMatrices(Tc, clusters)
	list_dfsGlm <- convertMatToDfForGlmFormat(list_matrices)

	list_lm <- runGlmForEachCluster(list_dfsGlm)
	# print("Note: finished creating glm's")

	list_resSummaries <- runPostHocTukeyForEachClus(list_lm)
	# print("Note: finished post-hoc tukey")

	# list_resSummaries <- createGlmTukeyForEachClus(list_dfsGlm)
	# return (list_resSummaries)


	class(list_resSummaries) <- "clusterChange"
	return (list_resSummaries)
}

#' @title
#' Summarize the Tukey constrasting results for consecutive intervals.
#'
#' @description
#' Summarize results obtained by running the \code{calcClusterChng} function. Specifically this function extracts z-scores and p-values for consecutive intervals.
#'
#' @param list_resSummaries A list returned by the \code{calcClusterChng} function.
#' @param Tc The time course data
#'
#' @return A list containing two items, where the first item is a z-score matrix, and the second item is a p-value matrix.
#'
#' @importFrom methods is
#'
#' @seealso \code{\link{calcClusterChng}}
#'
#' @export
summaryGetZP <- function(list_resSummaries, Tc){

	stopifnot(is(list_resSummaries,"clusterChange"), is(Tc, "matrix"))


	list_concTpSummary <- list()
	colNames = c()

	zScores <- matrix(nrow=length(list_resSummaries), ncol=ncol(Tc)-1)
	pValues <- matrix(nrow=length(list_resSummaries), ncol=ncol(Tc)-1)



	for(clustNum in 1:length(list_resSummaries)){
		for(tp in 2:ncol(Tc)){

			selName = paste(tp, "-", (tp-1))
			if (!is.null(colnames(Tc))){
				name = paste(colnames(Tc)[tp], "-", colnames(Tc)[(tp-1)])
			}
			else{
				name = selName
			}


			# print(name)
			if(clustNum == 1){
				colNames <- c(colNames, name)
			}


			zScores[clustNum, tp-1] <- list_resSummaries[[clustNum]]$test$tstat[selName]

			idx = which(names(list_resSummaries[[clustNum]]$test$tstat) == selName)

			pValues[clustNum, tp-1] <- list_resSummaries[[clustNum]]$test$pvalues[[idx]]

		}
	}


	colnames(zScores) <- colNames
	colnames(pValues) <- colNames

	list_concTpSummary[[1]] <- zScores
	list_concTpSummary[[2]] <- pValues

	class(list_concTpSummary) <- "summary.clusterChange"

	return(list_concTpSummary)
}

#' @title
#' Print changes in consecutive contrasting intervals
#'
#' @description
#' Prints two matrices, the first containing z-scores, and the second, p-values, where each row corresponds to each cluster, and each column corresponds to the time-interval.
#'
#' @param list_concTpSummary Z-scores and p-values extracted from tukey-constrasts for glm's of each cluster.
#'
#'
#' @seealso \code{\link{summaryGetZP}}
#'
#' @export
printZP <- function(list_concTpSummary){
	cat('===== list_[[1]] - zScores\n')
	print(list_concTpSummary[[1]])
	cat('\n\n')

	cat('===== list_[[2]] - pValues\n')
	print(list_concTpSummary[[2]])
	cat('\n\n')
}

#' @title
#' Plot a Z-score heatmap
#'
#' @description
#' Plot a heatmap of z-scores, masking the z-scores at insignificant p-values. Rows in the heatmap correspond to clusters, and columns in the heatmap correspond to time intervals.
#'
#' @param list_concTpSummary Z-scores and p-values extracted from tukey-constrasts of consecutive time intervals for each cluster. Obtained by running the \code{summaryGetZP} function.
#' @param significanceTh Tukey's p-value cutoff to be applied for significance.
#'
#' @return Displays a plot and returns a matrix, with z-scores masked (NA) according to the signficance (p-value) threshold specified.
#'
#'
#' @seealso \code{\link{summaryGetZP}}
#'
#' @importFrom gplots heatmap.2
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is
#'
#' @export
plotZP <- function(list_concTpSummary, significanceTh=0.05){

	stopifnot(is(list_concTpSummary, "summary.clusterChange"), (significanceTh >0 && significanceTh <= 1))

	theSignifs <- list_concTpSummary[[2]] < significanceTh

	matToPlot <- list_concTpSummary[[1]]

	rwb <- grDevices::colorRampPalette(colors = c("blue", "#cfcfcf", "red"))(n=299)


	if (!all(theSignifs)){
		zscore_largest <- ceiling(max(c(abs(min(matToPlot)), abs(max(matToPlot)))))
		theBreaks <- seq(-zscore_largest, zscore_largest, length.out=299)

		zscore_largest_noSignif <- max(c(abs(max(matToPlot[!theSignifs])), abs(min(matToPlot[!theSignifs]))))
		idx_forGray50 <- theBreaks <= zscore_largest_noSignif & theBreaks >= -zscore_largest_noSignif


		rwb[idx_forGray50] <- "#cfcfcf"

		matToPlot[!theSignifs] <- NA
	}


	titleTxt = paste("Significant changes\n(z-values of intervals,\n with significance p \u003C", significanceTh, ")", sep="")

	gplots::heatmap.2(matToPlot, main=titleTxt, xlab="Time interval", ylab="Cluster", Rowv=FALSE, Colv=FALSE, dendrogram="none", col=rwb, na.color="#cfcfcf", tracecol=NA, density.info="none", sepcolor="#cfcfcf", sepwidth=c(0.001, 0.001), colsep=0:ncol(matToPlot), rowsep=0:nrow(matToPlot), srtCol=45, cexCol=0.8)

	return (matToPlot)
}

plotZP_eachCluster <- function(list_resSummaries, Tc, glmTukeyForEachClus){


	# stopifnot(is(list_resSummaries,"clusterChange"), is(Tc, "matrix"))

	list_matsToPlot <- list()


	for (clusNum in 1:length(list_resSummaries)){
		mat_aClus <- matrix(nrow=ncol(Tc), ncol=ncol(Tc))

		# names_r <- c()
		# names_c <- c()
		# for (i in length(glmTukeyForEachClus[[clusNum]]$test$tstat)){

		# }

		for (tp2 in 1:ncol(Tc)){
			for (tp1 in tp2:ncol(Tc)){
				if (tp2 != tp1){
					# print(paste(tp1, '-', tp2))
					selName = paste(tp1, "-", (tp2))

					# names_r <- c(names_r, paste(colnames(Tc)[tp1], '-', colnames(Tc)[tp2]))

					mat_aClus[tp1, tp2] = as.numeric(glmTukeyForEachClus[[clusNum]]$test$tstat[selName])


					mat_aClus[tp2, tp1] = -1 * as.numeric(glmTukeyForEachClus[[clusNum]]$test$tstat[selName])
				}
				else{
					mat_aClus[tp1, tp2] = 0
				}
			}
		}
		# print(names_r)
		colnames(mat_aClus) <- colnames(Tc)
		rownames(mat_aClus) <- colnames(Tc)
		list_matsToPlot[[clusNum]] <- mat_aClus

	}
	doThePlotting(list_matsToPlot)
	return (list_matsToPlot)
}

getAbsMaxValue <- function(list_matsToPlot){

	largestAbsVal <- 0
	for (clusNum in 1:length(list_matsToPlot)){
		newMax <- max(abs(list_matsToPlot[[clusNum]]))

		if (newMax > largestAbsVal){
			largestAbsVal <- newMax
		}
	}

	return (largestAbsVal)
}

doThePlotting <- function(list_matsToPlot){

	pdf("allByAllHeatmaps.pdf", height=10, width=10)

	par(mfrow=c(ceiling(length(list_matsToPlot)/8), 8))

	rwb <- grDevices::colorRampPalette(colors = c("blue", "#cfcfcf", "red"))(n=99)

	largestAbsVal <- getAbsMaxValue(list_matsToPlot)
	zscore_largest <- ceiling(largestAbsVal)

	theBreaks <- sort(seq(-zscore_largest, zscore_largest, length.out=100))

	print(theBreaks)

	for (clusNum in 1:length(list_matsToPlot)){




		# zscore_largest_noSignif <- max(c(abs(max(matToPlot[!theSignifs])), abs(min(matToPlot[!theSignifs]))))
		# idx_forGray50 <- theBreaks <= zscore_largest_noSignif & theBreaks >= -zscore_largest_noSignif


		# rwb[idx_forGray50] <- "#cfcfcf"

		# matToPlot[!theSignifs] <- NA

		titleTxt = paste("Cluster ", clusNum, sep="")

		gplots::heatmap.2(list_matsToPlot[[clusNum]], main=titleTxt, xlab="Time interval", ylab="Time interval", Rowv=FALSE, Colv=FALSE, dendrogram="none", col=rwb, na.color="#cfcfcf", tracecol=NA, density.info="none", sepcolor="#cfcfcf", sepwidth=c(0.001, 0.001), colsep=0:ncol(list_matsToPlot[[clusNum]]), rowsep=0:nrow(list_matsToPlot[[clusNum]]), srtCol=45, cexCol=0.8, breaks=theBreaks)
	}
	dev.off()
}




##################################### PRIVATE FUNCTIONS
#' Run tukey post-hoc evaluations for the glm's.
#'
#' @param list_lms A list containing a lm for each cluster.
#' @return A list of tukey evaluations
#'
runPostHocTukeyForEachClus <- function(list_lms){
	list_tukeys <- list()

	for (clustNum in 1:length(list_lms)){
		# list_tukeys[[clustNum]] <- emmeans(list_lms[[clustNum]], pairwise ~ fac_timepoint, adjust="tukey")

		list_tukeys[[clustNum]] <- summary(multcomp::glht(list_lms[[clustNum]], multcomp::mcp(fac_timepoint="Tukey")))
	}

	return (list_tukeys)
}

#' Run glm for each cluster.
#'
#' @param list_dfsGlm A list of dataframes (for each cluster) of the time series data.
#' @return A list of glm's for each cluster.
#'
runGlmForEachCluster <- function(list_dfsGlm){

	list_lm <- list()

	for (dfNum in 1:length(list_dfsGlm)){
		aGlm <- list_dfsGlm[[dfNum]]
		list_lm[[dfNum]] <- stats::glm(formula=experimentalObs ~ fac_profileNum + fac_timepoint, data=aGlm)
	}

	return (list_lm)
}


#' Convert a list of matrices of each cluster, to a list of data frames (in the required glm format) for each cluster. Timepoint, and profileNum are factors, and experimentalObs (ratio) is y.
#'
#' @param list_matrices A list of matrices for each cluster
#' @return A list containing glm's for each cluster.
#'
convertMatToDfForGlmFormat <- function(list_matrices){
	list_dfsGlm <- list()

	# converting each matrix to a data_frame for GLM.
	for (clustNum in 1:length(list_matrices)){

		fac_timepoint = c()
		fac_profileNum = c()
		experimentalObs = c()

		for (tp in 1:ncol(list_matrices[[clustNum]])){

			fac_timepoint <- append(fac_timepoint, c(rep(tp, nrow(list_matrices[[clustNum]]))))


			fac_profileNum <- append(fac_profileNum, 1:nrow(list_matrices[[clustNum]]))

			experimentalObs <- append(experimentalObs, list_matrices[[clustNum]][,tp])

		}

		glmDf <- data.frame(fac_timepoint, fac_profileNum, experimentalObs)
		# add this df to the list
		# print(glmDf$fac_timepoint)
		glmDf$fac_timepoint <- factor(glmDf$fac_timepoint, ordered=TRUE)
		glmDf$fac_profileNum <- factor(glmDf$fac_profileNum) # , ordered=TRUE)

		list_dfsGlm[[clustNum]] <- glmDf

	}

	return(list_dfsGlm)
}


#' Split a matrix into cluster based submatrices.
#'
#' @param Tc A matrix containing time course data.
#' @param clusters A vector with same length as number of profiles. The cluster labels are assumed to be numbers, starting from 1.
#'
#' @return A list of matrices
#'
splitIntoSubMatrices <- function(Tc, clusters){

	list_matrices <- list()

	for (i in 1:max(clusters)){
		list_matrices[[i]] <- matrix(nrow=0, ncol=ncol(Tc))
	}

	for (i in 1:nrow(Tc)){
		rn <- rownames(Tc)[i]
		# print(rn)
		clusterOfRow <- clusters[[rn]]

		list_matrices[[clusterOfRow]] <- rbind(list_matrices[[clusterOfRow]], Tc[i,])

		rownames(list_matrices[[clusterOfRow]])[nrow(list_matrices[[clusterOfRow]])] <- rownames(Tc)[i]
	}

	return(list_matrices)
}
