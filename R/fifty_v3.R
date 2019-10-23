#### MAIN FUNCTION 2.
calc50crossing_v3 <- function(timeRegions, clustered, phosTh=0.5, dephosTh=0.5){

	# phosTh = as.numeric(phosTh)

	# stopifnot(is(phosTh>0 && phosTh<=1), is(dephosTh>0 && dephosTh<=1))
	mat_fiftyPoints <- matrix(ncol=6) # cluster, time, abundance, dir, startTp, endTp

	for (clusNum in 1:nrow(clustered$center)){
		for (i in 1:nrow(timeRegions[[clusNum]])){
			startTp = as.integer(min(timeRegions[[clusNum]][i, cols_timeReg$tpStart], timeRegions[[clusNum]][i, cols_timeReg$tpEnd]))
			endTp = as.integer(max(timeRegions[[clusNum]][i, cols_timeReg$tpStart], timeRegions[[clusNum]][i, cols_timeReg$tpEnd]))


			#print("Value here..?")
			print(paste(startTp, endTp))
			print(clustered$center[clusNum, startTp:endTp])

			crossPt <- calcCrossing_v3(clustered$center[clusNum, startTp:endTp], timeRegions[[clusNum]][i, cols_timeReg$dir], (startTp -1), phosTh, dephosTh)
			# print(crossPt[1])
			mat_fiftyPoints <- addToMatrixWithRegionInfo(mat_fiftyPoints, clusNum, crossPt[1], crossPt[2], timeRegions[[clusNum]][i, cols_timeReg$dir], startTp, endTp)
		}
	}

	mat_fiftyPoints <- mat_fiftyPoints[-1,,drop=FALSE]
	return (mat_fiftyPoints)
}

addToMatrixWithRegionInfo <- function(theMat, clusNum, x_50, y_50, direction, tpStart, tpEnd){
	theMat <- rbind(theMat, c(clusNum, x_50, y_50, direction, tpStart, tpEnd))
	return (theMat)
}



calcCrossing_v3 <- function(region, dir, offset, phosTh, dephosTh){
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

	print(paste("Value here is: ", y_max, y_min, y_50))

	# y_50 <- (y_max + y_min)/2

	x_50 <- -1 # initial time
	# print(c(y_max, y_50, y_min))

	# find intervals and store
	for (j in 2:length(region)){
		if (dir == 1 && y_50 >= region[j-1] && y_50 <= region[j]){
			# 50 is crossed here in phos direction
			x_50 <- getMidX((offset+j-1), region[j-1], offset+j, region[j], y_50)
			# mat_fiftyPoints <- addToMatrix(mat_fiftyPoints, clusNum, x_50, y_50, 1)
			# print(mat_fiftyPoints)
		}
		else if (dir == -1 && y_50 <= region[j-1] && y_50 >=region[j]) {
			# 50 is crossed here in dephos direction
			# getMidX()
			x_50 <- getMidX((offset+j-1), region[j-1], offset+j, region[j], y_50)
			# mat_fiftyPoints <- addToMatrix(mat_fiftyPoints, clusNum, x_50, y_50, -1)
		}
	}

	return (c(x_50, y_50))
}


#### MAIN FUNCTION 1.
getTimeRegionsWithMaximalChange <- function(glmTukeyForEachClus, numTps, signifTh=0.05){
	list_timeRegNoOvrlp <- list()

	for (i in 1:17){
		print("    ")
		print(paste("Cluster ", i))
		combined <- convertToMatAndAddNAs(glmTukeyForEachClus[[i]], numTps, signifTh)


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

convertToMatAndAddNAs <- function(glmTukeyForAClus, numTps, signifTh){
	mat_zscore <- matrix(nrow=(numTps-1), ncol=(numTps-1))
	mat_pvalue <- matrix(nrow=(numTps-1), ncol=(numTps-1))
	combined <- matrix(nrow=(numTps-1), ncol=(numTps-1)) # the row ranges from (2, 9), the col ranges from (1,8)


	for (val1 in seq(1, numTps)){
		if ((val1+1) <= numTps){
			for (val2 in seq((val1+1), (numTps))){
				name = paste(val1, "-", (val2))

				idx = which(glmTukeyForAClus$contrasts[1] == name )

				mat_zscore[val2-1, val1] <- as.numeric(glmTukeyForAClus$contrasts[idx,]$z.ratio) * -1

				# idx = which(names(zscores) == name)
				mat_pvalue[val2-1, val1] <- glmTukeyForAClus$contrasts[idx,]$p.value

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
