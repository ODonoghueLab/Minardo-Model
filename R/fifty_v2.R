calc50crossing_v2 <- function(timeRegions, clustered, phosTh=0.5, dephosTh=0.5){

}


#### MAIN FUNCTION 1.
getTimeRegionsWithMaximalChange <- function(glmTukeyForEachClus, numTps, signifTh=0.05){
	list_timeRegNoOvrlp <- list()

	for (i in 1:17){
		print("    ")
		print(paste("Cluster ", i))
		combined <- convertToMatAndAddNAs(glmTukeyForEachClus[[i]]$test$tstat, glmTukeyForEachClus[[i]]$test$pvalues, numTps, signifTh)


		printZPmat(combined)

		timeRegions <- findMonotonicRegions(combined, i)
		print("time regions (with overlaps)")
		print(timeRegions)
		# printTimeRegions(timeRegions)
		timeRegions_noOverlaps<- removeAnyOverlaps(timeRegions)
		# addToList(timeRegions_noOverlaps, somelistOfGlobals)
		print("time regions without overlaps")
		print(timeRegions_noOverlaps)

		list_timeRegNoOvrlp[[i]] <- timeRegions_noOverlaps
	}



	return(list_timeRegNoOvrlp)
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


	# print("overlaps PHOS")
	# print(noOverlaps_phos)

	# print("overlaps DEPHOS")
	# print(noOverlaps_dephos)
	# print("time regions")
	# print(timeRegions)
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
		print(combined[row, col])
		col = col - 1
	}

	return (aRow)
}


addToTimeRegionsMat <- function(timeRegions, cluster, zScore, tpStart, tpEnd, direction){
	tpStart = tpStart + 1
	timeRegions <- rbind(timeRegions, c(cluster, zScore, tpStart, tpEnd, direction))

	return(timeRegions)
}



########################### AUX
findMonotonicRegions <- function(combined, cluster){
	timeRegions <- matrix(ncol=5) # cluster, z-score, tpStart, tpEnd, direction

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

					# print(paste(rowStart, colStart, rowEnd+1, colEnd, combined[rowEnd, colEnd]))
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

convertToMatAndAddNAs <- function(zscores, pvalues, numTps, signifTh){
	mat_zscore <- matrix(nrow=(numTps-1), ncol=(numTps-1))
	mat_pvalue <- matrix(nrow=(numTps-1), ncol=(numTps-1))
	combined <- matrix(nrow=(numTps-1), ncol=(numTps-1)) # the row ranges from (2, 9), the col ranges from (1,8)


	for (val1 in seq(1, numTps)){
		if ((val1+1) <= numTps){
			for (val2 in seq((val1+1), (numTps))){
				name = paste(val2, "-", (val1))

				mat_zscore[val2-1, val1] <- zscores[name]
				idx = which(names(zscores) == name)
				mat_pvalue[val2-1, val1] <- pvalues[idx]

				if(as.numeric(pvalues[idx]) < signifTh){
					combined[val2-1, val1] <- zscores[name]
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
