#### MAIN FUNCTION 3.


#' @title
#' Overlays the 50\% abundance crossings on the Z-score heatmap.
#'
#' @param list_concTpSummary Z-scores and p-values extracted from tukey-constrasts for glm's of each cluster. Obtained by running the "summaryGetZP" function.
#' @param significanceTh P-value cutoff. This value ranges between 0 and 1. Generally, the threshold of p<0.5 is considered significant. This threshold can be reduced to, for example 0.01, 0.05 or 0.001, for only plotting the highly significant z-scores. To plot all z-scores, the treshold can be set to 1.
#' @param mat_fiftyPoints A matrix containing the event information generated by running the "calc50crossing" function.
#'
#'
#'
#' @return Displays a plot and returns a matrix, with z-scores masked (NA) according to the signficance threshold specified.
#'
#' @importFrom methods is
#' @importFrom gplots heatmap.2
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics segments
#' @importFrom shape Arrows
#'
#' @seealso \code{\link{plotZP}}, \code{\link{summaryGetZP}}, \code{\link{calc50crossing}}
#'
#' @export
plotZP_withEvents <- function(list_concTpSummary, mat_fiftyPoints, significanceTh=0.5){

	stopifnot( (significanceTh >0 && significanceTh <= 1), is(mat_fiftyPoints, "matrix"))

	segs <- calcSegsOfCentroids(mat_fiftyPoints) # need to add to global, for graphics::segments to find.
	assign("segs", segs, envir=.GlobalEnv)
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
	# print(titleTxt)

	idx_up = segs[[5]] == "#8B0000"
	assign("idx_up", idx_up, envir=.GlobalEnv)
	idx_down = segs[[5]] == "#00008B"
	assign("idx_down", idx_down, envir=.GlobalEnv)


	gplots::heatmap.2(matToPlot, main=titleTxt, xlab="Time interval", ylab="Cluster", Rowv=FALSE, Colv=FALSE, dendrogram="none", col=rwb, na.color="#cfcfcf", tracecol=NA, density.info="none", sepcolor="#cfcfcf", sepwidth=c(0.001, 0.001), colsep=0:ncol(matToPlot), rowsep=0:nrow(matToPlot), srtCol=45, add.expr=c(shape::Arrows(x0=segs[[1]][idx_up], x1=segs[[1]][idx_up], y0=segs[[2]][idx_up]-0.1,  y1=segs[[2]][idx_up]+0.5, arr.type="triangle", arr.length=0.2, arr.width=0.2, col="#8B0000", arr.adj=0, segment=FALSE), shape::Arrows(x0=segs[[1]][idx_down], x1=segs[[1]][idx_down], y0=segs[[2]][idx_down]+2,  y1=segs[[2]][idx_down]+0.5, arr.type="triangle", arr.length=0.2, arr.width=0.2, col="#00008B", arr.adj=0, segment=FALSE), graphics::segments(x0=segs[[1]], y0=segs[[2]]+0.3, x1=segs[[1]], y1=segs[[4]]-0.3, col=segs[[5]], lwd=3, lend=2)), cexCol=0.8)


	return (matToPlot)
}

addToMatrix <- function(theMat, clusNum, x_50, y_50, direction){
	theMat <- rbind(theMat, c(clusNum, x_50, y_50, direction))
	return (theMat)
}


calcSegsOfCentroids <- function(mat_fiftyPoints){

	totalNumCluster <- max(mat_fiftyPoints[,cols_matFifty$col_clus]) + 0.5
	y_val = totalNumCluster

	x <- c(); y <- c(); x1 <- c(); y1 <- c(); colors <- c();

	for (rowNum in 1:nrow(mat_fiftyPoints)){
		x[rowNum] <- mat_fiftyPoints[rowNum, cols_matFifty$col_x] - 0.5
		x1[rowNum] <- mat_fiftyPoints[rowNum, cols_matFifty$col_x] - 0.5
		y[rowNum] <- ((totalNumCluster - (mat_fiftyPoints[rowNum, cols_matFifty$col_clus])) * 1)
		y1[rowNum] <- ((totalNumCluster - (mat_fiftyPoints[rowNum, cols_matFifty$col_clus] - 1)) * 1)

		if (mat_fiftyPoints[rowNum, cols_matFifty$col_dir] == 1){
			colors[rowNum] <- "#8B0000"
		}
		else{
			colors[rowNum] <- "#00008B"
		}
	}

	return(list(x, y, x1, y1, colors))
}


filterEvents <- function(mat_fiftyPoints, resWithOnlySignif, phosZscoreTh=10, dephosZscoreTh=10, asPercent=TRUE){

	if (asPercent == TRUE){
		# convert resWithOnlySignif to percentages

		minVal = min(resWithOnlySignif, na.rm = T)
		maxVal = max(resWithOnlySignif, na.rm = T)

		resWithOnlySignif <- (resWithOnlySignif - minVal)/ (maxVal - minVal) * 100
	}

	for (rowNum in 1:nrow(resWithOnlySignif)){
		for (colNum in 1:ncol(resWithOnlySignif)){

			if (!is.na(resWithOnlySignif[rowNum, colNum]) && resWithOnlySignif[rowNum, colNum] > 0 && resWithOnlySignif[rowNum, colNum] < phosZscoreTh){
				resWithOnlySignif[rowNum, colNum] <- 0
				removeIfAnyEventAffected(mat_fiftyPoints, rowNum, colNum)
			}

			if (!is.na(resWithOnlySignif[rowNum, colNum]) && resWithOnlySignif[rowNum, colNum] < 0 && resWithOnlySignif[rowNum, colNum] > dephosZscoreTh){
				resWithOnlySignif[rowNum, colNum] <- 0
				removeIfAnyEventAffected(mat_fiftyPoints, rowNum, colNum)

				print(resWithOnlySignif[rowNum, colNum])
			}
		}
	}




	return (resWithOnlySignif)
}

removeIfAnyEventAffected <- function(){

}

getMidX <- function(x1, y1, x2, y2, y_50){
	m <- (y2 - y1)/(x2 - x1)
	x_50 <- x1 + ((y_50 - y1)/m)
	return (x_50)
}


#### MAIN FUNCTION 2.
calcEvents <- function(timeRegions, clustered, phosEventTh=0.5, dephosEventTh=0.5){

	stopifnot(is(clustered, "fclust"), (phosEventTh >= 0 && phosEventTh <= 1), (dephosEventTh >= 0 && dephosEventTh <= 1))

	mat_fiftyPoints <- matrix(ncol=6) # cluster, time, abundance, dir, startTp, endTp

	for (clusNum in 1:nrow(clustered$center)){
		for (i in 1:nrow(timeRegions[[clusNum]])){
			startTp = as.integer(min(timeRegions[[clusNum]][i, cols_timeReg$tpStart], timeRegions[[clusNum]][i, cols_timeReg$tpEnd]))
			endTp = as.integer(max(timeRegions[[clusNum]][i, cols_timeReg$tpStart], timeRegions[[clusNum]][i, cols_timeReg$tpEnd]))

			crossPt <- calcCrossing_v3(clustered$center[clusNum, startTp:endTp], timeRegions[[clusNum]][i, cols_timeReg$dir], (startTp -1), phosEventTh, dephosEventTh)

			# if (crossPt[1] == -1){
				# print(clustered$center[clusNum, startTp:endTp])
			# }
			mat_fiftyPoints <- addToMatrixWithRegionInfo(mat_fiftyPoints, clusNum, crossPt[1], crossPt[2], timeRegions[[clusNum]][i, cols_timeReg$dir], startTp, endTp)
		}
	}

	mat_fiftyPoints <- mat_fiftyPoints[-1,,drop=FALSE]

	mat_fiftyPoints <- mat_fiftyPoints[order(mat_fiftyPoints[,1],mat_fiftyPoints[,2]),] # order by cluster, then time; else issue in downstream function (as only subsequent events are checked over there)

	colnames(mat_fiftyPoints) <- c("Cluster", "Time", "Th. abundance", "Direction", "Start interval", "End interval")

	return (mat_fiftyPoints)
}

addToMatrixWithRegionInfo <- function(theMat, clusNum, x_50, y_50, direction, tpStart, tpEnd){
	theMat <- rbind(theMat, c(clusNum, x_50, y_50, direction, tpStart, tpEnd))
	return (theMat)
}



calcCrossing_v3 <- function(region, dir, offset, phosTh, dephosTh){

	x_50 <- -1 # initial time
	y_50 <- -1 # initial th_point

	y_max <- max(region)
	y_min <- min(region)

	if (dir == 1){
		y_10Percent <- (y_max - y_min) * phosTh
		y_50 <- y_min + y_10Percent
	}
	else{
		y_10Percent <- (y_max - y_min) * (1- dephosTh)
		y_50 <- y_min + y_10Percent
	}


	e = 0.0001 # due to calculation some round error when threshold is 0.
	# find intervals and store
	for (j in 2:length(region)){

		if (dir == 1 && ( (y_50 >= region[j-1] && y_50 <= region[j]) || abs(y_50 - region[j-1]) < e  || abs(y_50 - region[j]) < e )){
			# 50 is crossed here in phos direction
			x_50 <- getMidX((offset+j-1), region[j-1], offset+j, region[j], y_50)
			# mat_fiftyPoints <- addToMatrix(mat_fiftyPoints, clusNum, x_50, y_50, 1)
		}
		else if (dir == -1 && ( (y_50 <= region[j-1] && y_50 >=region[j]) || abs(y_50 - region[j-1]) < e  || abs(y_50 - region[j]) < e ) ){
			# 50 is crossed here in dephos direction
			# getMidX()
			x_50 <- getMidX((offset+j-1), region[j-1], offset+j, region[j], y_50)
			# mat_fiftyPoints <- addToMatrix(mat_fiftyPoints, clusNum, x_50, y_50, -1)

		}
	}

	return (c(x_50, y_50))
}

missingStats <- function(Tc, clustered, mat_fiftyPoints, phosTh, dephosTh){
	list_matrices <- splitIntoSubMatrices(Tc, clustered)

	list_distributions <- getDistOfAllEvents_v2(mat_fiftyPoints, list_matrices, phosTh, dephosTh)


	list_distributions
	mat_missingStats <- matrix(ncol=7) #cluster, time, th.abun., dir, numNa, numTotal, percentNa

	for (eventNum in 1:length(list_distributions)){

		numNa <-length(which(is.na(list_distributions[[eventNum]][,cols_matFifty$col_x])))

		numTotal <- nrow(list_distributions[[eventNum]])

		percentNa <- (numNa * 100) / numTotal

		mat_missingStats <- rbind(mat_missingStats, c(mat_fiftyPoints[eventNum, cols_matFifty_v2$clus:cols_matFifty_v2$dir], numNa, numTotal, percentNa))

		# print(paste(eventNum, numNa, numTotal, percentNa, sep=" "))
	}

	colnames(mat_missingStats) <- c("Cluster", "Time", "Th. abun.", "Dir.", "Num. NA", "Num. total", "Percent NA")

	mat_missingStats <- mat_missingStats[-1,,drop=FALSE]
	return(mat_missingStats)
}


#### MAIN FUNCTION 1.
getTimeRegionsWithMaximalChange <- function(glmTukeyForEachClus, numTps, signifTh=0.05, phosZscoreTh=0, dephosZscoreTh=0){

	# Add checks
	stopifnot(is(glmTukeyForEachClus, "clusterChange"), numTps > 0, (signifTh > 0 && signifTh <= 1))

	list_timeRegNoOvrlp <- list()

	for (i in 1:length(glmTukeyForEachClus)){ # for each cluster

		# convert z-scores to matrix (by time * time)
		combined <- convertToMatAndAddNAs(glmTukeyForEachClus[[i]], numTps, signifTh)


		# Convert regions to 0 which cannot statisfy zscoreTh.
		combined <- filterByZscore(combined, phosZscoreTh, dephosZscoreTh, asPercent)

		# find time regions
		timeRegions <- findMonotonicRegions(combined, i)

		timeRegions_noOverlaps<- removeAnyOverlaps(timeRegions)


		list_timeRegNoOvrlp[[i]] <- timeRegions_noOverlaps
	}



	return(list_timeRegNoOvrlp)
}

filterByZscore <- function(combined, phosZscoreTh, dephosZscoreTh, asPercent){
	# combined_local <- combined


	for (rowNum in 1:nrow(combined)){
		for (colNum in 1:ncol(combined)){

			if (!is.na(combined[rowNum, colNum]) && combined[rowNum, colNum] > 0 &&combined[rowNum, colNum] < phosZscoreTh){
				combined[rowNum, colNum] <- 0

			}

			if (!is.na(combined[rowNum, colNum]) && combined[rowNum, colNum] < 0 && combined[rowNum, colNum] > dephosZscoreTh){

				combined[rowNum, colNum] <- 0

			}
		}
	}



	return (combined)
}

isAnyOverlap <- function(noOverlaps, currTp1, currTp2){

	if (nrow(noOverlaps) == 0){
		return (FALSE)
	}

	for (i in 1: nrow(noOverlaps)) {
		if (noOverlaps[i,cols_timeReg$tpStart] > min(c(currTp1,currTp2)) && noOverlaps[i,cols_timeReg$tpStart]  < max(c(currTp1, currTp2))) {
			return (TRUE)
		}
		if ((noOverlaps[i,cols_timeReg$tpEnd]  > min(c(currTp1,currTp2)) && noOverlaps[i,cols_timeReg$tpEnd]  <  max(c(currTp1, currTp2))) || (noOverlaps[i,cols_timeReg$tpEnd]  > currTp1 && noOverlaps[i,cols_timeReg$tpEnd]  < currTp2)){
			return (TRUE)
		}
		if (currTp1 > min(c(noOverlaps[i,cols_timeReg$tpStart],noOverlaps[i,cols_timeReg$tpEnd]))  && currTp2 < max(c(noOverlaps[i,cols_timeReg$tpStart], noOverlaps[i,cols_timeReg$tpEnd]) )) {
			return (TRUE)
		}
		if (currTp2 > min(c(noOverlaps[i,cols_timeReg$tpStart],noOverlaps[i,cols_timeReg$tpEnd]))  && currTp2 < max(c(noOverlaps[i,cols_timeReg$tpStart], noOverlaps[i,cols_timeReg$tpEnd]) ) ){
			return (TRUE)
		}
	}

	return (FALSE)
}

removeAnyOverlaps <- function(timeRegions){
	noOverlaps_phos <- matrix(ncol=5)
	noOverlaps_dephos <- matrix(ncol=5)
	isFirstAdded_phos <- FALSE
	isFirstAdded_dephos <- FALSE
	for (i in 1:nrow(timeRegions)){
		if (nrow(timeRegions) == 0){
			break;
		}
		if ( any(timeRegions[,cols_timeReg$zScore] > 0)){
			maxIdx <- which.max(timeRegions[,cols_timeReg$zScore])
			if (isFirstAdded_phos == FALSE){
				noOverlaps_phos <- rbind(noOverlaps_phos, timeRegions[maxIdx,])
				isFirstAdded_phos = TRUE
				noOverlaps_phos <- noOverlaps_phos[-1,,drop=FALSE ]
				# print(timeRegions)
			}
			else{
				if (isAnyOverlap(noOverlaps_phos, timeRegions[maxIdx,cols_timeReg$tpStart], timeRegions[maxIdx,cols_timeReg$tpEnd]) == FALSE){
					noOverlaps_phos <- rbind(noOverlaps_phos, timeRegions[maxIdx,])
				}
			}
			timeRegions <- timeRegions[-maxIdx,,drop=F]
		}

		if (nrow(timeRegions) > 0 && any(timeRegions[,cols_timeReg$zScore] < 0)){
			minIdx <- which.min(timeRegions[,cols_timeReg$zScore])
			if (isFirstAdded_dephos == FALSE){
				noOverlaps_dephos <- rbind(noOverlaps_dephos, timeRegions[minIdx,])
				isFirstAdded_dephos = TRUE
				noOverlaps_dephos <- noOverlaps_dephos[-1,,drop=FALSE]
			}
			else{
				if (isAnyOverlap(noOverlaps_dephos, timeRegions[minIdx,cols_timeReg$tpStart], timeRegions[minIdx,cols_timeReg$tpEnd]) == FALSE){
					noOverlaps_dephos <- as.matrix(rbind(noOverlaps_dephos, timeRegions[minIdx,]))
				}
			}
			timeRegions <- timeRegions[-minIdx,,drop=FALSE]
		}
	}

	if (is.na(noOverlaps_phos[1,1] )){
		return (noOverlaps_dephos)
	}
	if (is.na(noOverlaps_dephos[1,1])){
		return (noOverlaps_phos)
	}

	overlaps <- rbind(noOverlaps_phos, noOverlaps_dephos)

	return (overlaps)

}

printTimeRegions <- function(timeRegionsMat){
	if (length(timeRegionsMat) == 0 || nrow(timeRegionsMat) == 0){
		return ();
	}
	else{
		for (i in 1:nrow(timeRegionsMat)){
			print(timeRegionsMat[i,])
		}
	}
}


traverseRow <- function(combined, row, col, currMaxMin){
	aRow = c()

	while(col > 0){
		# print(combined[row, col])
		col = col - 1
	}

	return (aRow)
}


addToTimeRegionsMat <- function(timeRegions, cluster, zScore, tpStart, tpEnd, direction){
	# dir = -1
	# if (zScore > 0){
	#	dir = 1
	# }
	tpStart = tpStart + 1
	# startReg = min(tpStart, tpEnd)
	# endReg = max(tpStart, tpEnd)

	timeRegions <- rbind(timeRegions, c(cluster, zScore, tpStart, tpEnd, direction))

	return(timeRegions)
}



########################### AUX
findMonotonicRegions <- function(combined, cluster){
	timeRegions <- matrix(ncol=5) # cluster, z-score, tpStart, tpEnd, direction, z-score

	direction = 0
	currMaxMin = 0
	rowStart = 0
	colStart = 0
	rowEnd = 0
	colEnd = 0


	# trial 3
	rowNum =1
	while (rowNum <= nrow(combined)){
		colNum = 1
		isBreak = FALSE
		while (colNum <= rowNum){
			if (combined[rowNum, colNum] > 0 && direction != 1) {
				# save any prev value ...
				if (currMaxMin < 0){
					timeRegions <- addToTimeRegionsMat(timeRegions, cluster, combined[rowEnd, colEnd], rowEnd, colEnd, direction)

					# print("Set 1: ")
					# print(paste(cluster, rowStart, colStart, rowEnd+1, colEnd, combined[rowEnd, colEnd], direction))
				}
				# new one encountered in + dir.
				direction = 1
				currMaxMin = combined[rowNum, colNum]
				rowStart = rowNum; rowEnd = rowNum; colStart = colNum; colEnd = colNum;


			}
			else if (combined[rowNum, colNum] > 0 && direction == 1 ){
				# encountered in same direction

				if (combined[rowNum, colNum] > currMaxMin){
					currMaxMin = combined[rowNum, colNum]
					rowEnd = rowNum; colEnd = colNum;
				}
			}


			if (combined[rowNum, colNum] < 0 && direction != -1) {
				# save any prev value ...
				if (currMaxMin > 0){
					# print(paste(rowStart, colStart, rowEnd+1, colEnd, combined[rowEnd, colEnd]))
					timeRegions <- addToTimeRegionsMat(timeRegions, cluster, combined[rowEnd, colEnd], rowEnd, colEnd, direction)

					# print("Set 2: ")
					# print(paste(cluster, rowStart, colStart, rowEnd+1, colEnd, combined[rowEnd, colEnd], direction))
				}
				# new one encountered in + dir.
				direction = -1
				currMaxMin = combined[rowNum, colNum]
				rowStart = rowNum; rowEnd = rowNum; colStart = colNum; colEnd = colNum;
				#print("here! ")
				#print(combined[rowNum, colNum])
			}
			else if (combined[rowNum, colNum] < 0 && direction == -1 ){
				# encountered in same direction
				if (combined[rowNum, colNum] < currMaxMin){
					currMaxMin = combined[rowNum, colNum]
					rowEnd = rowNum; colEnd = colNum;
				}
			}

			colNum = colNum + 1
		}


		rowNum = rowNum + 1


	}
	if (combined[rowEnd, colEnd] != 0){
		timeRegions <- addToTimeRegionsMat(timeRegions, cluster, combined[rowEnd, colEnd], rowEnd, colEnd, direction)

		# print(paste(rowStart, colStart, rowEnd+1, colEnd, combined[rowEnd, colEnd]))
	}




	timeRegions <- timeRegions[-1,, drop=FALSE]
	return (timeRegions)
}



printZPmat <- function(zscore){
	for (row in 1:nrow(zscore)){
		for (col in 1:ncol(zscore)){

			cat(zscore[row, col])

			if (col < ncol(zscore)){
				cat('\t')
			}
		}
		cat('\n')
	}

}

convertToMatAndAddNAs <- function(glmTukeyForAClus, numTps, signifTh){
	mat_zscore <- matrix(nrow=(numTps-1), ncol=(numTps-1))
	mat_pvalue <- matrix(nrow=(numTps-1), ncol=(numTps-1))
	combined <- matrix(nrow=(numTps-1), ncol=(numTps-1)) # the row ranges from (2, 9), the col ranges from (1,8)


	for (val1 in seq(1, numTps)){
		if ((val1+1) <= numTps){
	#		for (val2 in seq((val1+1), (numTps))){
			for (val2 in seq((val1+1), numTps)) {

				name = paste(val2, "-", (val1))
				# print(name)

				idx = which(names(glmTukeyForAClus$test$tstat) == name )
				# print(paste(idx , name))
				# print(glmTukeyForAClus$test$tstat[idx])
		# 		print(glmTukeyForAClus$contrasts[idx,]$z.ratio)


				mat_zscore[val2-1, val1] <- as.numeric(glmTukeyForAClus$test$tstat[idx]) # * -1

				# idx = which(names(zscores) == name)

				mat_pvalue[val2-1, val1] <- glmTukeyForAClus$test$pvalues[[idx]]

				if(as.numeric(mat_pvalue[val2-1, val1]) < signifTh){
					combined[val2-1, val1] <- mat_zscore[val2-1, val1]
				}
				else{
					combined[val2-1, val1] <- 0
				}
				# print(name)
				# print(idx)

			}
		}
	}



	return (combined)
}