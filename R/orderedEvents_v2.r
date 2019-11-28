#' @title
#' Orders clusters based on time of the 50\% percent abundance crossings
#'
#' @description
#' Here, the events identified from the centroids are used to generate time-distributions using profiles in the clusters. From these distributions, a matrix (either by running t-test or wilcox) is generated of the p-values and the t-statistic (or the 95\% confidence interval, respectively). The p-values are FDR corrected and are used to define a graph, where events with significant time differences are ordered according to the t-statistic (for t-test) or the 95\% confidence interval (for wilcox test). This matrix is transitively reduced and is used to generate the ordering and the resulting figure.
#'
#' @param Tc A matrix containing time course data.
#' @param clustered A fclust object containing the clustering details.
#' @param mat_fiftyPoints A matrix containing the event information generated by running the "calc50crossing" function.
#' @param test The test to run to determine the orders. The choice is between a parametric or a non-paramtric test:
#' \describe{
#' 		\item{\code{wilcox}}{Choosing this will result in the computation of 2 NxN matrices containing the 95\% confidence interval value for direction and the p-values, where N is the number of events}
#' 		\item{\code{t-test}}{Choosing this will result in the computation of 2 NxN matrices containing the t-statistic for direction and the p-values}
#' }
#' @param fdrSignif The significant cutoff for FDR correction (performed on a p-value matrix (see above)).
#'
#'
#' @return A figure containing the ordered clusters by occurance of first event is generated. Along with this, returned is a matrix, similar to the mat_fiftyPoints matrix, where an additional column is appended containing the order of occurance of the corresponding event.
#'
#'
#' @importFrom methods is
#' @importFrom stats p.adjust wilcox.test t.test
#' @importFrom nem transitive.reduction
#' @importFrom igraph graph_from_adjacency_matrix adjacent_vertices V are.connected vertex_connectivity
#' @importFrom graphics plot segments axis rect text lines
#' @importFrom shape Arrows
#'
#' @seealso \code{\link[e1071]{cmeans}} for clustering time profiles, \code{\link{calc50crossing}} for identifying events.
#'
#' @export
orderTheEvents <- function(Tc, clustered, mat_fiftyPoints, test="wilcox", fdrSignif=0.05, phosTh=0.5, dephosTh=0.5){ # } exclTh=0){
	test_param = "t-test"
	test_nonParam = "wilcox"

	title_testType = ""

	stopifnot(is(Tc, "matrix"), is(clustered, "fclust"), is(mat_fiftyPoints, "matrix"))

	if (!(test == test_param) && !(test == test_nonParam)){
		stop(paste("Test ", test, " not recognized.", sep=""))
	}

	if (fdrSignif < 0 || fdrSignif > 1){
		stop(paste("fdrsignif value is out of range.", sep=""))
	}



	list_matrices <- splitIntoSubMatrices(Tc, clustered)

	list_distributions <- getDistOfAllEvents_v2(mat_fiftyPoints, list_matrices, phosTh, dephosTh)

	# mat_missingStats <- missingStats(list_distributions) # exclude events based on percentage missing. mat_fiftyPoints, list_distributions
	# return(mat_missingStats)

	# return (list_distributions)


	list_pVal_stat <- list()

	if (test==test_nonParam){
		print ("Running wilcoxon test.")
		list_pVal_stat <- performWilcoxonSignedRankTabular(list_distributions)


		title_testType = "(Non-parametric)"
	}
	else {
		print("Running t-test.")
		list_pVal_stat <- performTtestsTabular(list_distributions)

		title_testType = "(Parametric)"
	}

	# FDR adjustment for multiple testing.
	mat_pValsFDR <- matrix(stats::p.adjust(as.vector(list_pVal_stat[[1]]), method='fdr'), ncol=nrow(mat_fiftyPoints))

	# Convert to adjacency matrix.
	signifs <- mat_pValsFDR < fdrSignif

	statistic_t <- list_pVal_stat[[2]]
	statistic_t[!signifs] <- NA


	statistic_t[statistic_t > 0] <- NA
	statistic_t[is.na(statistic_t)] <- 0
	statistic_t[statistic_t < 0] <- 1

	# library(nem) # for trasitive reduction
	# library(igraph) # for adjacency matrix
	reducted <- nem::transitive.reduction(statistic_t)
	net = igraph::graph_from_adjacency_matrix(reducted, mode="directed")

	list_eventsOrder <- getEventsOrderInGraph(reducted, net) # computing the order of the events in the mat_fiftyPoints based on a bfs search.
	list_eventsOrder <- adjustEvents(list_eventsOrder, signifs)


	mat_fiftyPts_withOrder <- appendOrder(mat_fiftyPoints, list_eventsOrder)



	## Getting other things ready for plotting

	# plotting constants
	sigEventsDiff_x = 0.5
	nonSigEventsDiff_x = 0.25

	eventPtGraph_y = -1
	yLabelInit = -1.8
	yLabelSpace = 0.5

	lineDiff_y = 1

	# graphXStart = -1
	eventStart_x = 0

	blockSpacingFmEvent_x = sigEventsDiff_x/2
	blockStart_x = eventStart_x - blockSpacingFmEvent_x


	#x_lineStart, y_lineStart, y_lineDecr, x_incr_sig, x_incr_nonSig

	yPos_eventPts = -1

	# plotting values
	list_blocks <- getRectBlock(list_eventsOrder, signifs)
	mat_eventPoints <- getTheClusLines(mat_fiftyPoints, list_eventsOrder, signifs, list_blocks, eventStart_x, sigEventsDiff_x, nonSigEventsDiff_x, lineDiff_y)

	eventEnd_x = max(as.numeric(mat_eventPoints[,cols_matFifty$col_x])) #  + blockSpacingFmEvent_x
	blockEnd_x = eventEnd_x + blockSpacingFmEvent_x

	mat_bgRects <- getRectPoints(mat_eventPoints, list_blocks, blockStart_x, blockSpacingFmEvent_x)


	# mat_eventPoints <- adjustByBgRects(mat_bgRects, list_eventsOrder, mat_eventPoints)


	mat_grayLines <- getGrayLines(mat_eventPoints, signifs)
	mat_grayLines_v2 <- getGrayLines_v2(list_eventsOrder, mat_eventPoints, signifs)

	# print("The max x is ")
	# print(max(as.numeric(mat_eventPoints[,cols_matFifty$col_x])))
	mat_clusConnLines <- getClusConnLines(mat_eventPoints, blockStart_x, blockEnd_x)
	print("The mat clus conn lines")
	print(mat_clusConnLines)
	vec_labels <- getYaxisLabels(mat_eventPoints)


	mat_xLabelsClus <- getXaxisLabels(list_eventsOrder, mat_eventPoints, signifs, list_blocks, eventStart_x, sigEventsDiff_x, nonSigEventsDiff_x, yLabelInit, yLabelSpace);

	# constants


	# a set of points to set up the graph.
	graphics::plot(mat_eventPoints[,cols_matFifty$col_x], mat_eventPoints[,cols_matFifty$col_y], asp=NA, yaxt="n", lwd=0.25, col=colors_orderedEvents$incr, pch=".", xlab="Temporal order", ylab="Clusters",  main=paste("Temporal order of events in clusters", title_testType, sep=""), xlim=c(eventStart_x, eventEnd_x), ylim=c(-4, max(as.numeric(mat_eventPoints[,cols_matFifty$col_clus]))), xaxt="n", bty="n", cex.main = 0.8) #,

	# background gray rectangles.
	if (nrow(mat_bgRects) > 0){
		for (rowNum in 1:nrow(mat_bgRects)){
			if (rowNum%%2 == 1){
				graphics::rect(mat_bgRects[rowNum,1], lineDiff_y/2, mat_bgRects[rowNum,2], max(as.numeric(mat_eventPoints[,cols_matFifty$col_clus])) + (lineDiff_y/2), col = "#ededed", border=NA)
			}

		}
	}

	# axis labels on sides 2 and 4.
	graphics::axis(side=2, las=1, at=seq(1,length(vec_labels)), labels=rev(vec_labels), cex=0.02,  col = NA ) #, col.ticks = 1)

	graphics::axis(side=4, las=1, at=seq(1,length(vec_labels)), labels=rev(vec_labels), cex=0.02,  col = NA ) #, col.ticks = 1)


	graphics::segments(x0=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_x0]), y0=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_y0]), x1=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_x1]), y1=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_y1]), col="#808080", lty="dotted", lwd=0.5)
	# C99999



	# rect(-0.5, 0, 0.5, 17.5, col = "yellow", border="#ededed")

	# the cluster lines
	graphics::segments(x0=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_x0]), y0=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_y0]), x1=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_x1]), y1=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_y1]), col=mat_clusConnLines[,cols_clusPlotObjs$col_col], lwd=16, lend=2)


	# the events (arrows) within a graph
	mat_upEvents = mat_eventPoints[mat_fiftyPoints[,cols_matFifty$col_dir] == 1,]
	mat_downEvents = mat_eventPoints[mat_fiftyPoints[,cols_matFifty$col_dir] == -1,]


	# print(mat_upEvents)


	## arrow bottoms
	shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01-0.2,  y1=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01, arr.type="triangle", arr.length=0.01, arr.width=0.01, col=color_events$up, lwd=2, lend=2, arr.adj=0)

	shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25,  y1=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25-0.2, arr.type="triangle", arr.length=0.01, arr.width=0.01, col=color_events$down, lwd=2, lend=2, arr.adj=0)


	## sharp arrow heads
	shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01-0.2,  y1=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01, arr.type="triangle", arr.length=0.13, arr.width=0.15, col=color_events$up, arr.adj=0, segment=FALSE)

	shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25,  y1=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25-0.2, arr.type="triangle", arr.length=0.13, arr.width=0.15, col=color_events$down, arr.adj=0, segment=FALSE)


	# gray lines in events points maps
	straightGrayLines = mat_grayLines_v2[mat_grayLines_v2[,cols_grayLines_v2$col_isSemiCirc] == FALSE,]
	graphics::segments(x0=as.numeric(straightGrayLines[,cols_grayLines_v2$col_x0]), y0=rep(yPos_eventPts, length(straightGrayLines)), x1=as.numeric(straightGrayLines[,cols_grayLines_v2$col_x1]), y1=rep(yPos_eventPts, length(straightGrayLines)), lwd=1, col="#A2A2A2") #, lty="dotted")

	# gray curves in events points map
	mat_curvePoints = mat_grayLines_v2[mat_grayLines_v2[,cols_grayLines_v2$col_isSemiCirc] == TRUE,,drop=FALSE]

	if (length(mat_curvePoints) > 0 && nrow(mat_curvePoints) > 0){
		for (i in 1:nrow(mat_curvePoints)){
			list_pts = calcCurves(mat_curvePoints[i, cols_grayLines_v2$col_x0], mat_curvePoints[i, cols_grayLines_v2$col_x1], -1)
			graphics::lines(list_pts[[1]], list_pts[[2]], col="#A2A2A2")
		}
	}




	# printing up and down events (in arc events map)
	y_adjust = 0.15

	## arrow bottom lines
	shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=rep(-1.2, nrow(mat_upEvents)), y1=rep(-1, nrow(mat_upEvents)), arr.type="triangle", arr.length=0.01, arr.width=0.01, col=color_events$up, lwd=2, lend=2, arr.adj=0)

	shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=rep(-0.8, nrow(mat_downEvents)), y1=rep(-1, nrow(mat_downEvents)), arr.type="triangle", arr.length=0.01, arr.width=0.01, col=color_events$down, lwd=2, lend=2, arr.adj=0)

	## sharp arrow heads
	shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=rep(-1.2, nrow(mat_upEvents)), y1=rep(-1, nrow(mat_upEvents)), arr.type="triangle", arr.length=0.13, arr.width=0.15, col=color_events$up, arr.adj=0, segment=FALSE)

	shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=rep(-0.8, nrow(mat_downEvents)), y1=rep(-1, nrow(mat_downEvents)), arr.type="triangle", arr.length=0.13, arr.width=0.15, col=color_events$down, arr.adj=0, segment=FALSE)


	# vec_clusEvent = mat_eventPoints[,4]
	# print(vec_clusEvent)
	# vec_clusEvent[mat_eventPoints[mat_fiftyPoints[,cols_matFifty$col_dir] == 1]] = color_events$up
	# vec_clusEvent[mat_eventPoints[mat_fiftyPoints[,cols_matFifty$col_dir] == -1] = color_events$down

	# x label text.
	graphics::text(x=as.numeric(mat_xLabelsClus[,1]), y=as.numeric(mat_xLabelsClus[,2]), labels=mat_xLabelsClus[,3], offset=0, cex=0.5, col=mat_xLabelsClus[,4])

	return (mat_fiftyPts_withOrder)
}

excludeEvents <- function(mat_missingStats, list_distributions, mat_fiftyPoints, exclTh){
	list_ <- list()

	idxToExcl <- which(mat_missingStats[,cols_missingStats$percentNa] > exclTh)

	list_[[1]] <- mat_missingStats[-idxToExcl,]
	list_[[2]] <- list_distributions[-idxToExcl]
	list_[[3]] <- mat_fiftyPoints[-idxToExcl,]

	return(list_)
}

missingStats <- function(list_distributions){
	mat_missingStats <- matrix(ncol=3) #numNa, numTotal, percentNa

	for (eventNum in 1:length(list_distributions)){

		numNa <-length(which(is.na(list_distributions[[eventNum]][,cols_matFifty$col_x])))

		numTotal <- nrow(list_distributions[[eventNum]])

		percentNa <- (numNa * 100) / numTotal

		mat_missingStats <- rbind(mat_missingStats, c(numNa, numTotal, percentNa))

		# print(paste(eventNum, numNa, numTotal, percentNa, sep=" "))
	}

	mat_missingStats <- mat_missingStats[-1,,drop=FALSE]
	return(mat_missingStats)
}


calcCurves <- function(x0, x1, yCentPos){

	x_10pts <- seq(x0, x1, length.out=50)

	y_10pts <- topSemiCircleEq(x_10pts, yCentPos)

	list_pts <- list(x_10pts, y_10pts)
	return (list_pts)

}

topSemiCircleEq <- function(xPts, yCentPos){
	yPts = c()


	r = abs(xPts[1] - xPts[length(xPts)])/2
	xCentPos = (xPts[1] + xPts[length(xPts)]) /2

	for (xPt in xPts){
		yPt = sqrt(abs(r^2 - (xPt-xCentPos)^2)) + yCentPos + 0.4
		yPts = c(yPts, yPt)
	}

	return (yPts)
}


appendOrder <- function(mat_fiftyPoints, list_eventsOrder){

	mat_withOrder <- cbind(mat_fiftyPoints, rep(0,nrow(mat_fiftyPoints)))

	for (i in 1:length(list_eventsOrder)){
		eventNum <- list_eventsOrder[[i]]
		mat_withOrder[eventNum, cols_orderedEvents$col_order] = i
	}

	return(mat_withOrder)
}

################################## AUXILLARY
adjustEvents <- function(list_eventsOrder, signifs){
	list_newEventsOrder <- list()
	i_new = 1

	for (i in 1:length(list_eventsOrder)){
		# print(list_eventsOrder[[i]])
		isAddedVec = FALSE
		for (j in 1:length(list_eventsOrder[[i]])){

			if (isAddedVec == FALSE){ # add every row atleast once
				# add to list first
				list_newEventsOrder[[length(list_newEventsOrder) + 1]] <- list_eventsOrder[[i]]
				isAddedVec = TRUE
				# print(list_newEventsOrder)
			}

			if (length(list_eventsOrder[[i]]) == 1){ # if only one element break
				break
			}

			if (isAnyOneNonSignifFromNextAll(list_eventsOrder, list_eventsOrder[[i]][j], i, signifs)){ # remove that element, and add to end.
				# add the jth one to a new next layer
				# print("here!")
				# print(list_eventsOrder[[i]][j])


				list_newEventsOrder <- rmElemAndAddSepElseNew(list_newEventsOrder, list_eventsOrder[[i]][j])
				# print("now list is")
				# print(list_newEventsOrder)

			}
		}

	}

	return (list_newEventsOrder)
}

checkAndAddEventAtEnd <- function(anEvent, list_newEventsOrder){
	isFound = FALSE
	for (i in 1:length(list_newEventsOrder)){
		if (anEvent %in% list_newEventsOrder){
			isFound = TRUE
		}
	}

	if (isFound == FALSE){
		list_newEventsOrder <- append(list_newEventsOrder, c(anEvent))
	}

	return(list_newEventsOrder)
}

rmElemAndAddSepElseNew <- function(list_newEventsOrder, theElem, theVec){
	# print("The eleme is ")
	# print(theElem)
	for (i in 1:length(list_newEventsOrder)){
		if (theElem %in% list_newEventsOrder[[i]]){
			idx = match(theElem, list_newEventsOrder[[i]])
			# print(" the idx is ")
			# print(idx)
			if (idx > 0){
				list_newEventsOrder[[i]] <- list_newEventsOrder[[i]][-idx]
			}


		}
	}

	list_newEventsOrder <- append(list_newEventsOrder, c(theElem))

	return (list_newEventsOrder)
}

isVecInList <- function(theVec, theList, elemToRm){

	for (subVec in theList){
		if (all(subVec==theVec)){
			return (TRUE)
		}
	}

	return (FALSE)
}

isAnyOneNonSignifFromNextAll <- function(list_eventsOrder, currEvent, currentLayer, signifs){
	if ((currentLayer + 1) >= length(list_eventsOrder)){
		return (FALSE)
	}

	for (i in seq((currentLayer + 1), length(list_eventsOrder))){

		isAny = isAnyNonSignifFromPrev(c(currEvent), list_eventsOrder[[i]], signifs)

		if (isAny == TRUE){
			return (TRUE)
		}
	}

	return (FALSE)
}



getEventsOrderInGraph <- function(reducted, net){
	vec_roots <- getRoots(reducted) # get root nodes
	vec_linCurrNodes <- vec_roots

	# list_currNodes <- as.list(vec_roots)


	vec_visited <- c()

	list_eventsOrder <- list()


	while(length(vec_linCurrNodes) > 0){

		# list_prevNodes = list_currNodes


		# print("vec_linCurrNodes")
		# print (vec_linCurrNodes)

		vec_newForVisit <- getNotVisitedBefore(vec_linCurrNodes, vec_visited)

		vec_newForVisit <- onlyPatron(vec_newForVisit, net)


		# print("to visit")
		# print (vec_newForVisit)


		if (length(vec_newForVisit) > 0){
			# print("to visit")
			# print (vec_newForVisit)

			list_eventsOrder <- addNodesToOrderList(vec_newForVisit, list_eventsOrder)
			vec_visited <- visitNodes(vec_newForVisit, vec_visited)
		}



		list_currNodes <- getAllInNextHop(net, vec_newForVisit, "out")
		vec_linCurrNodes <- convertToUniqueVec(list_currNodes)

		if (length(vec_linCurrNodes) == 0 && length(igraph::V(net)) != length(vec_visited)){
			vec_linCurrNodes <- setdiff(igraph::V(net), vec_visited)
		}
	}
	return (list_eventsOrder)
}





cols_clusPlotObjs$col_x0 <- 1
cols_clusPlotObjs$col_y0 <- 2
cols_clusPlotObjs$col_x1 <- 3
cols_clusPlotObjs$col_y1 <- 4
cols_clusPlotObjs$col_col <- 5





getNumOfCrossingAndDir <- function(mat_fiftyPoints, clusNum){

	totalNumOfCrossings <- 0
	list_dir = list()

	for (rowNum in 1:nrow(mat_fiftyPoints)){
		if (mat_fiftyPoints[rowNum, cols_matFifty$col_clus] == clusNum){
			totalNumOfCrossings = totalNumOfCrossings + 1
			list_dir = append(list_dir, mat_fiftyPoints[rowNum, cols_matFifty$col_dir])
		}
	}

	return (list(totalNumOfCrossings, list_dir))
}


getDistOfAllEvents_v2 <- function(mat_fiftyPoints, list_matrices, phosTh, dephosTh){
	list_distributions <- list()

	for (eventNum in 1:nrow(mat_fiftyPoints)){
		dist <- genDistForEvent(list_matrices[[mat_fiftyPoints[eventNum,cols_matFifty_v2$clus]]], mat_fiftyPoints[eventNum, cols_matFifty_v2$dir], mat_fiftyPoints[eventNum, cols_matFifty_v2$startTp], mat_fiftyPoints[eventNum, cols_matFifty_v2$endTp], phosTh, dephosTh)

		list_distributions[[length(list_distributions) + 1]] <-  dist
	}

	return(list_distributions)
}

genDistForEvent <- function(clusMat, dir, startTp, endTp, phosTh, dephosTh){
	mat_distribution <- matrix(ncol=4)


	for (profNum in 1:nrow(clusMat)){

		if ((dir == 1 && clusMat[profNum, startTp] < clusMat[profNum, endTp]) || (dir == -1 && clusMat[profNum, startTp] > clusMat[profNum, endTp])){
			# direction holds; do calculation, and save.
			# print(clusMat[profNum, startTp:endTp])

			crossPt <- calcCrossing_v3(clusMat[profNum, startTp:endTp], dir, (startTp -1), phosTh, dephosTh)

			mat_distribution <- addToMatrix(mat_distribution, profNum, crossPt[1], crossPt[2], dir)

		}
		else{
			mat_distribution <- addToMatrix(mat_distribution, profNum, NA, NA, dir)
		}
	}

	mat_distribution <- mat_distribution[-1,, drop=FALSE]
	return (mat_distribution)
}

getDistOfAllEvents <- function(mat_fiftyPoints, list_matrices){
	list_distributions <- list()

	for (clusNum in 1:length(list_matrices)){
		totalCrossingsAndDir <- getNumOfCrossingAndDir(mat_fiftyPoints, clusNum)

		list_distributions <- genDistributions(list_matrices[[clusNum]], totalCrossingsAndDir[[1]], totalCrossingsAndDir[[2]], list_distributions)
		# print(dist)


	}

	return (list_distributions)
}

genDistributions <- function(cluster, numCrossings, list_dir, list_distributions){


	for (crossNum in 1:numCrossings){

		mat_crossings <- matrix(ncol=4)

		for (rowNum in 1:nrow(cluster)){

			eventCounter = 0

			y_max <- max(cluster[rowNum,])
			y_min <- min(cluster[rowNum,])

			y_50 <- (y_max + y_min)/2

			for (j in  2:ncol(cluster)){

				if ((y_50 >= cluster[rowNum,j-1] && y_50 <= cluster[rowNum,j]) || y_50 <= cluster[rowNum,j-1] && y_50 >= cluster[rowNum,j]){ # if phos or dephos
					eventCounter = eventCounter + 1
					# print('here!')
					if (y_50 >= cluster[rowNum,j-1] && y_50 <= cluster[rowNum,j] && list_dir[crossNum] == 1 && eventCounter == crossNum){ # 50 is crossed here in phos direction
						x_50 <- getMidX((j-1), cluster[rowNum,j-1], j, cluster[rowNum,j], y_50)
						mat_crossings <- addToMatrix(mat_crossings, rowNum, x_50, y_50, 1)
						break
					}
					else if (y_50 <= cluster[rowNum,j-1] && y_50 >= cluster[rowNum,j] && list_dir[crossNum] == -1 && eventCounter == crossNum) { # 50 is crossed here in dephos direction
						x_50 <- getMidX((j-1), cluster[rowNum,j-1], j, cluster[rowNum,j], y_50)
						mat_crossings <- addToMatrix(mat_crossings, rowNum, x_50, y_50, -1)
						break
					}

				}
			}

		}

		mat_crossings <- mat_crossings[-1,]
		list_distributions[[length(list_distributions)+1]] <-  mat_crossings
	}

	return (list_distributions)
}


isComparisonPossible <- function(list1, list2){
	list1 <- list1[!is.na(list1)]
	list2 <- list2[!is.na(list2)]

	# print(list1)

	if(all(abs(list1 - mean(list1)) < 0.00001) && all(abs(list2 - mean(list2)) < 0.00001) && abs(mean(list1) - mean(list2)) < 0.00001) {
		return (FALSE)
		# print("false!!!)")
	}


	return (TRUE)
}


performWilcoxonSignedRankTabular <- function(list_distributions){
	mat_pVal = matrix(nrow=length(list_distributions), ncol=length(list_distributions))
	mat_statistic = matrix(nrow=length(list_distributions), ncol=length(list_distributions))

	for (i in 1:length(list_distributions)){
		for (j in 1:length(list_distributions)){

			if (i != j){

				if (!isComparisonPossible(list_distributions[[i]][,cols_matFifty$col_x], list_distributions[[j]][,cols_matFifty$col_x])){ # add as pval = 1 (being the exact same from a uniform distribution).
					mat_pVal[i, j] <- 1
					mat_statistic[i, j] <- NA
				}
				else{ # do the comparison
					res <- stats::wilcox.test(list_distributions[[i]][,cols_matFifty$col_x], list_distributions[[j]][,cols_matFifty$col_x], paired=FALSE, conf.int=TRUE)
					# return (res)
					# print (res)
					mat_pVal[i, j] <- res$p.value
					mat_statistic[i, j] <- res$estimate
				}
			}
		}
	}

	list_pVal_stat <- list()
	list_pVal_stat[[1]] <- mat_pVal
	list_pVal_stat[[2]] <- mat_statistic
	return (list_pVal_stat)
}



performTtestsTabular <- function(list_distributions){
	mat_pVal = matrix(nrow=length(list_distributions), ncol=length(list_distributions))
	mat_statistic = matrix(nrow=length(list_distributions), ncol=length(list_distributions))

	for (i in 1:length(list_distributions)){
		for (j in 1:length(list_distributions)){

			if (i != j){
				res <- stats::t.test(list_distributions[[i]][,cols_matFifty$col_x], list_distributions[[j]][,cols_matFifty$col_x], paired=FALSE)
				mat_pVal[i, j] <- res$p.value
				mat_statistic[i, j] <- res$statistic
			}
		}
	}

	list_pVal_stat <- list()
	list_pVal_stat[[1]] <- mat_pVal
	list_pVal_stat[[2]] <- mat_statistic
	return (list_pVal_stat)
}


################################ GRAPH_TRAVERSAL

getRoots <- function(reducted){
	vec_roots = c()
	for (colNum in 1:nrow(reducted)){

		isAllZeros = TRUE
		for (rowNum in 1:ncol(reducted)){
			if (reducted[rowNum, colNum] == 1){
				isAllZeros = FALSE
				break
			}
		}

		if (isAllZeros == TRUE){
			vec_roots <- append(vec_roots, colNum)
		}

	}

	return (vec_roots)
}



getNotVisitedBefore <- function(vec_currNodes, vec_visited){
	vec_newForVisit <- c()

	for (i in 1:length(vec_currNodes)){

		if (!(vec_currNodes[i] %in% vec_visited) && !(vec_currNodes[i] %in% vec_newForVisit)){
			vec_newForVisit <- c(vec_newForVisit, vec_currNodes)
		}
	}

	return (vec_newForVisit)
}


onlyPatron <- function(vec_nodes, net){

	# only add a node if it is not a child of anything else in this set.

	vec_parentsOnly <- c()

	for (i in 1:length(vec_nodes)){
		isChild = FALSE;

		vec_withoutCurr <- setdiff(vec_nodes, vec_nodes[i])

		if (length(vec_withoutCurr) > 0) {
			list_children <- getAllInNextHop(net, vec_withoutCurr, 'out')
			vec_children <- convertToUniqueVec(list_children)

			# print(vec_children)
			# print(vec_nodes[i])

			if (vec_nodes[i] %in% vec_children){
				isChild = TRUE
			}
			else { # make sure that this node is also not in the downstream lineage.
				if (length(vec_children) > 0){
					for (j in 1:length(vec_children)){
						if (igraph::are.connected(net, vec_children[j], vec_nodes[i]) == TRUE || igraph::vertex_connectivity(net, vec_children[j], vec_nodes[i]) == 1){
							isChild = TRUE
						}
					}
				}
			}

		}

		if (isChild == FALSE){
			vec_parentsOnly <- c(vec_parentsOnly, vec_nodes[i])
		}

	}

	return (vec_parentsOnly)
}


addNodesToOrderList <- function(vec_nodesInNewLayer, list_eventsOrder){
	pos <- length(list_eventsOrder) + 1
	list_eventsOrder[[pos]] <- vec_nodesInNewLayer
	return (list_eventsOrder)
}



visitNodes <- function(vec_newForVisit, vec_visited){

	vec_visited <- c(vec_visited, vec_newForVisit)

	return (vec_visited)
}



getAllInNextHop <- function(network, list_nodes, edgeDir){
	list_ <- igraph::adjacent_vertices(network, list_nodes, mode="out")


	return (list_)
}


convertToUniqueVec <- function(list_nodes){
	vec_nodes <- c()
	# print("lvl1")
	if (length(list_nodes) > 0){
		for (i in 1:length(list_nodes)){
			# print("lvl2")
			if (length(list_nodes[[i]]) > 0){
				for(j in 1:length(list_nodes[[i]])){
					# print(c(i,j), sep="\t")
					# print(list_nodes)
					if (! (list_nodes[[i]][j] %in% vec_nodes)){
						vec_nodes <- c(vec_nodes, list_nodes[[i]][j])
					}
				}
			}
		}
	}
	# print(list_nodes)
	return (vec_nodes)
}




###################### FOR PLOTTING
getXaxisLabels <- function(list_eventsOrder, mat_eventPoints, signifs, list_blocks, eventStart_x, xSigIncr, xNonSigIncr, yLabelInit, yLabelSpace){
	mat_xLabelsClus = matrix(ncol=4) #x0, y0, labelNames, color
	x0 = eventStart_x;
	eventsVisited = c()

	isFirst = TRUE

	for (i in 1:length(list_eventsOrder)){

		eventsVisited = c(eventsVisited, list_eventsOrder[[i]])

		if (i > 1 && isAnyNonSignifFromPrevAll(list_eventsOrder, i, signifs) == TRUE || isInSameBlock(list_eventsOrder[[i]], eventsVisited, list_blocks)){
			if (isFirst == FALSE){
				x0 = x0 + xNonSigIncr;
			}
			else{
				isFirst = FALSE
			}
		}
		else {
			if (isFirst == FALSE){
				x0 = x0 + xSigIncr;
			}
			else {
				isFirst = FALSE
			}

		}

		y0 = yLabelInit

		sortedEvents <- getSortedEventsByY(list_eventsOrder[[i]], mat_eventPoints)

		for (event in sortedEvents){
			col = color_events$up

			if (mat_eventPoints[event, cols_matFifty$col_dir] != colors_orderedEvents$incr){
				col = color_events$down
			}
			rowVal = c(x0, y0, mat_eventPoints[event, cols_matFifty$col_clus], col)
			mat_xLabelsClus <- rbind(mat_xLabelsClus, rowVal)

			y0 = y0 - yLabelSpace;
		}


	}

	mat_xLabelsClus <- mat_xLabelsClus[-1,]
	return (mat_xLabelsClus)
}

getSortedEventsByY <- function(events, mat_eventPoints){
	theOrder <- rev(order(as.numeric(mat_eventPoints[events, cols_matFifty$col_y])))
	return(events[theOrder])
}


getGrayLines_v2 <- function(list_eventsOrder, mat_eventPoints, signifs){
	theGrayLines = matrix(ncol=4) # x0, x1, isSemiCircle, color

	for (i in 1:nrow(signifs)){
		for (j in i:ncol(signifs)){
			if (i != j){
				if (signifs[j, i] == FALSE){
					# if in same layer no need for line

					j_layer = getLayerNum(list_eventsOrder, j)
					i_layer = getLayerNum(list_eventsOrder, i)

					if (j_layer != i_layer){
						if (i_layer - 1 == j_layer || j_layer - 1 == i_layer){
							rowVal <- c(mat_eventPoints[i, cols_matFifty$col_x], mat_eventPoints[j, cols_matFifty$col_x], FALSE, "#808080")
						}
						else{
							rowVal <- c(mat_eventPoints[i, cols_matFifty$col_x], mat_eventPoints[j, cols_matFifty$col_x], TRUE, "#808080")
						}
						theGrayLines <- rbind(theGrayLines, rowVal)
					}
					# if in diff layer
						# if in consecutive layer, draw straight line

						# if not, draw semicircle.


					# draw a gray50 line between the i and j events




				}
			}
		}
	}
	theGrayLines <- theGrayLines[-1,]

	# print(theGrayLines)
	return (theGrayLines)
}

getLayerNum <- function(list_eventsOrder, event){
	for (i in 1:length(list_eventsOrder)){
		if (event %in% list_eventsOrder[[i]]){
			return (i)
		}
	}

	return (-1)
}

getGrayLines <- function(mat_eventPoints, signifs){
	theGrayLines = matrix(ncol=5)


	for (i in 1:nrow(signifs)){
		for (j in i:ncol(signifs)){
			if (i != j){
				if (signifs[j, i] == FALSE){
					# draw a gray50 line between the i and j events
					rowVal <-  c(mat_eventPoints[i, cols_matFifty$col_x], mat_eventPoints[i, cols_matFifty$col_y], mat_eventPoints[j, cols_matFifty$col_x], mat_eventPoints[j, cols_matFifty$col_y], "red")
					theGrayLines <- rbind(theGrayLines, rowVal)
				}
			}
		}
	}
	theGrayLines <- theGrayLines[-1,]

	return (theGrayLines)
}



getRectPoints <- function(mat_eventPoints, list_blocks, rectStart, rectSpacingFmEvent){ #xStart, incr
	mat_rectPts = matrix(ncol=2) # [x0, x1]

	for (i in 1:length(list_blocks)){
		xVal <- max(as.numeric(mat_eventPoints[list_blocks[[i]], cols_matFifty$col_x])) + rectSpacingFmEvent

		if (i == 1){
			x0 <- rectStart
		}
		else{
			x0 <- mat_rectPts[nrow(mat_rectPts), cols_rectPoints$x1]
		}

		mat_rectPts <- rbind(mat_rectPts, c(x0, xVal))

	}
	mat_rectPts <- mat_rectPts[-1, ]
	return (mat_rectPts)
}


getRectBlock <- function(list_eventsOrder, signifs){
		list_blocks <- list()

		prevBlock <- c()
		currBlock <- list_eventsOrder[[1]]


		if(length(list_eventsOrder) > 1) {
			for (i in 2:length(list_eventsOrder)){


				isAny = isAnyNonSignifFromPrev(list_eventsOrder[[i]], currBlock, signifs)

				if (isAny == TRUE){
					currBlock <- c(currBlock, list_eventsOrder[[i]])

				}
				else{

					prevBlock <- c(prevBlock, currBlock)
					list_blocks[[length(list_blocks)+1]] <- currBlock
					currBlock <- list_eventsOrder[[i]]
					# rectPoints <- addToRectPoints(rectPoints, list_eventsOrder[[i]], mat_eventPoints)

					 # can improve further - i.e. take max of prevBlock listEvents x value (done)
				}

				# if (i == length(list_eventsOrder)){
			#		if (list_blocks[[length(list_blocks)]])
			#	}

				if (length(prevBlock) > 0 && length(currBlock) > 0 && isAnyNonSignifFromPrev(currBlock, prevBlock, signifs) == TRUE){

					prevBlock <- prevBlock[!prevBlock == list_blocks[[length(list_blocks)]]]
					# prevBlock <- prevBlock[!prevBlock == joiningEvents]
					currBlock <- c(currBlock, list_blocks[[length(list_blocks)]])
					list_blocks <- list_blocks[-length(list_blocks)]

				}
			}

			# print(currBlock)
			if (any(currBlock %in% list_blocks[[length(list_blocks)]]) == TRUE && all(currBlock %in% list_blocks[[length(list_blocks)]]) == FALSE ){
				list_blocks[[length(list_blocks)]] <- currBlock
			}
			else if (all(currBlock %in% list_blocks[[length(list_blocks)]]) == FALSE){
				list_blocks[[length(list_blocks)+1]] <- currBlock
			}

		}
		# print("list blocks in getRectBlock")
		# print(list_blocks)

		# print(rectPoints)
		return (list_blocks)

}


addToRectPoints <- function(rectPoints, events, mat_eventPoints){
	# print(max(as.numeric(mat_eventPoints[events, cols_matFifty$col_x])))
	maxX <- max(as.numeric(mat_eventPoints[events, cols_matFifty$col_x]))

	rowVal = c()
	if (nrow(rectPoints) <= 1){
		rowVal <- c(-1, maxX-0.5)

	}
	else{
		rowVal <- c(rectPoints[nrow(rectPoints), 2], maxX-0.5)
	}
	# print("here")
	# print(rowVal)

	rectPoints <- rbind(rectPoints, rowVal)
	return (rectPoints)
}



getClusConnLines <- function(mat_eventPoints, xStart, xEnd){
	mat_clusConnLines = matrix(ncol=5,nrow=nrow(mat_eventPoints))
	mat_basalConnLines = matrix(ncol=5,nrow=max(as.numeric(mat_eventPoints[,cols_matFifty$col_clus]))) # only need 1 for each cluster

	basalCounter = 1
	prevClus = -1

	for(i in 1:nrow(mat_eventPoints)){


		# adding start point
		mat_clusConnLines[i, cols_clusPlotObjs$col_x0] <- mat_eventPoints[i, cols_matFifty$col_x]
		mat_clusConnLines[i, cols_clusPlotObjs$col_y0] <- mat_eventPoints[i, cols_matFifty$col_y]

		# adding end y
		mat_clusConnLines[i, cols_clusPlotObjs$col_y1] <- mat_eventPoints[i, cols_matFifty$col_y]


		# adding end x
		mat_clusConnLines[i, cols_clusPlotObjs$col_x1] <- getTheEndX(mat_eventPoints, i, xEnd)


		# adding color
		mat_clusConnLines[i, cols_clusPlotObjs$col_col] <- mat_eventPoints[i, cols_matFifty$col_dir]


		## handling the adding of the basal
		if (prevClus == -1 ){
			# just store everything
			mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_x0] = xStart;
			mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_x1] = mat_clusConnLines[i, cols_clusPlotObjs$col_x0];
			mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_y0] = mat_clusConnLines[i, cols_clusPlotObjs$col_y0];
			mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_y1] = mat_clusConnLines[i, cols_clusPlotObjs$col_y1];

			mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_col]  <- getOppColor(mat_eventPoints[i, cols_matFifty$col_dir])

			prevClus = mat_eventPoints[i,cols_matFifty$col_clus]
			# print(prevClus)
		}
		else if (prevClus == as.integer(mat_eventPoints[i,cols_matFifty$col_clus])){
			# update x if smaller (and then also color to be opposite)

			if (as.numeric(mat_clusConnLines[i, cols_clusPlotObjs$col_x0]) < as.numeric(mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_x1]) ){
				mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_x0] = mat_clusConnLines[i, cols_clusPlotObjs$col_x0];
				mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_col]  <- getOppColor(mat_eventPoints[i, cols_matFifty$col_dir])
			}
		}
		else if (prevClus !=  mat_eventPoints[i,cols_matFifty$col_clus]){
			prevClus = mat_eventPoints[i,cols_matFifty$col_clus]
			basalCounter = basalCounter + 1
			# store everything again
			mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_x0] = xStart;
			mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_x1] = mat_clusConnLines[i, cols_clusPlotObjs$col_x0];
			mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_y0] = mat_clusConnLines[i, cols_clusPlotObjs$col_y0];
			mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_y1] = mat_clusConnLines[i, cols_clusPlotObjs$col_y1];

			mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_col]  <- getOppColor(mat_eventPoints[i, cols_matFifty$col_dir])


		}

	}

	mat_clusConnLines = rbind(mat_clusConnLines, mat_basalConnLines)
	# print(mat_clusConnLines)
	# print(mat_eventPoints)
	return (mat_clusConnLines)

}

getOppColor <- function (thisCol){
	if (thisCol == colors_orderedEvents$incr){
		return (colors_orderedEvents$decr)
	}

	return (colors_orderedEvents$incr)
}

getTheEndX <- function(mat_eventPoints, eventNum, defaultEndX){
	endX <- defaultEndX

	startX <- mat_eventPoints[eventNum, cols_matFifty$col_x]

	for (i in 1:nrow(mat_eventPoints)){
		if (i != eventNum){ # if not the current event
			if (mat_eventPoints[i, cols_matFifty$col_clus] == mat_eventPoints[eventNum, cols_matFifty$col_clus]){ # if same cluster
				if (as.numeric(mat_eventPoints[i, cols_matFifty$col_x]) > as.numeric(startX) && as.numeric(mat_eventPoints[i, cols_matFifty$col_x]) < as.numeric(endX)){ # if newX > startX && < endX; then update
					endX <- mat_eventPoints[i, cols_matFifty$col_x]
				}

			}
		}
	}

	return (endX)
}

getYaxisLabels <- function(mat_eventPoints){
	ord <- order(as.numeric(mat_eventPoints[,cols_matFifty$col_x]))
	clusNums <- as.numeric(mat_eventPoints[ord,1])
	clusNumLabels <- c()

	for(i in 1:length(clusNums)){
		if (! clusNums[[i]] %in% clusNumLabels){
			clusNumLabels <- append(clusNumLabels, clusNums[[i]])
		}
	}
	vec_labels =  as.vector(clusNumLabels)

	return (vec_labels)
}


getYIfClusPres <- function(currClusNum, mat_eventPoints){

	for (rowNum in 1:nrow(mat_eventPoints)){
		if (!is.na(mat_eventPoints[rowNum, cols_matFifty$col_clus]) && currClusNum == mat_eventPoints[rowNum, cols_matFifty$col_clus]){
			return (mat_eventPoints[rowNum, cols_matFifty$col_y])
		}
	}

	return (-1)
}

isAnyNonSignifFromPrevAll <- function(list_eventsOrder, currLayerNum, signifs){

	for (i in seq(1, (currLayerNum - 1))){
		isAny = isAnyNonSignifFromPrev(list_eventsOrder[[currLayerNum]], list_eventsOrder[[i]], signifs)
		if (isAny == TRUE){
			return (TRUE)
		}
	}

	return (FALSE)
}




getTheClusLines <- function(mat_fiftyPoints, list_eventsOrder, signifs, list_blocks, xStart, xIncrSig, xIncrNonSig, yDecr){ #xStart, yStart, y_decr, x_incr_sig, x_incr_nonSig


	mat_eventPoints <- matrix(nrow=nrow(mat_fiftyPoints), ncol=4)

	x <- xStart
	y <- max(mat_fiftyPoints[,1])

	eventsVisited = c()

	for (i in 1:length(list_eventsOrder)){
		for (j in 1:length(list_eventsOrder[[i]])){
			eventNum <- list_eventsOrder[[i]][j]

			# adding y coordinate
			prevY <- getYIfClusPres(mat_fiftyPoints[eventNum, cols_matFifty$col_clus], mat_eventPoints)
			if (prevY == -1){
				mat_eventPoints[eventNum, cols_matFifty$col_y] <- y
				y <- y - yDecr
			}
			else{
				mat_eventPoints[eventNum, cols_matFifty$col_y] <- prevY
			}

			# assigning color
			if (mat_fiftyPoints[eventNum, cols_matFifty$col_dir] == 1){
				mat_eventPoints[eventNum, cols_matFifty$col_dir] = colors_orderedEvents[['incr']]
			}
			else{
				mat_eventPoints[eventNum, cols_matFifty$col_dir] = colors_orderedEvents[['decr']]
			}

			# assigning cluster number
			mat_eventPoints[eventNum, cols_matFifty$col_clus] <- mat_fiftyPoints[eventNum, cols_matFifty$col_clus]

			eventsVisited <- c(eventsVisited, eventNum)
		}

		# adding x coordinate
		if (i == 1){
			#print("over here..?")
			mat_eventPoints <- addXForAll(list_eventsOrder[[i]], mat_eventPoints, x)
		}
		else{
			# if (isAnyNonSignifFromPrev(list_eventsOrder[[i]], list_eventsOrder[[i-1]], signifs) == TRUE){ # place events closer together
			if (isAnyNonSignifFromPrevAll(list_eventsOrder, i, signifs) == TRUE || isInSameBlock(list_eventsOrder[[i]], eventsVisited, list_blocks)){
				x <- x + xIncrNonSig
			}
			else {
				x <- x + xIncrSig
			}

			mat_eventPoints <- addXForAll(list_eventsOrder[[i]], mat_eventPoints, x)

		}
	}

	# print("The mat event points:")
	# print (mat_eventPoints)
	return (mat_eventPoints)

}

isInSameBlock <- function(events, eventsVisited, list_blocks){

	# 1. find block.
	# 2. exclude events from block, and from eventsVisited.
	# 3. See if others in block, which are also visited.
	# if so, return true, else return false.



	blockNum <- findBlock(events[1], list_blocks)


	eventsInBlock <-  as.vector(list_blocks[[blockNum]])
	otherEventsInBlock <- setdiff(eventsInBlock, events)

	otherEventsVisited <- setdiff(eventsVisited, events)

	if (length(otherEventsInBlock) > 0 && length(otherEventsVisited) > 0 && (any(otherEventsVisited %in% otherEventsInBlock) || any(otherEventsInBlock %in% otherEventsVisited))) {
		return (TRUE)
	}

	return (FALSE)
}

findBlock <- function(event, list_blocks){
	i = -1
	for (i in 1:length(list_blocks)){
		if (event %in% list_blocks[[i]]){
			return (i)
		}
	}

	return (i)
}



isAnyNonSignifFromPrev <- function(vec_curr, vec_prev, signifs){

	for (i in 1:length(vec_curr)){
		for (j in 1:length(vec_prev)){

			if (signifs[vec_curr[i], vec_prev[j]] == FALSE){
				return (TRUE)
			}

		}
	}

	return (FALSE)
}

addXForAll <- function(vec_curr, mat_eventPoints, x){

	for (i in 1:length(vec_curr)){
		mat_eventPoints[vec_curr[i], cols_matFifty$col_x] <- x
	}

	return (mat_eventPoints)
}
