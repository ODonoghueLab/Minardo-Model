#' @title
#' Calculate the ordering of events
#'
#' @description
#' For each event, we generate distributions of time of event occurance using profiles within clusters. These distributions are then pair wise compared and the results of all comparisons are represented in the graph space, where events are the graph nodes, and the directional edges represent that the source node event occurs significantly early than the target node. This graph is then utilized to generate the ordering.
#'
#' @param Tc A matrix containing time course data.
#' @param clusters A vector with same length as number of profiles. The cluster labels are assumed to be numbers, starting from 1.
#' @param mat_events A matrix containing the time regions (and events for the centroid) generated by running the \code{calcEvents} function.
#' @param test The test used to compare the event time distributions. The choice is between a parametric or a non-parametric test:
#' \describe{
#' 		\item{\code{wilcox}}{Choosing this will result in event times being compared using the Mann Whitney U test}
#' 		\item{\code{t-test}}{Choosing this will result in event times being compared using the Students t-test}
#' }
#' @param fdrSignif The significance cutoff for FDR correction, performed after pairwise comparisons of the event times.
#' @param phosEventTh Value between 0 and 1, to define the threshold at which a phosphorylation (or increasing) event occurs. Here, a threshold of 0 corresponds to the minimum value within an interval. (Such an event is designated as 1).
#' @param dephosEventTh Value between 0 and 1, to define the threshold at which a dephosphorylation (or decreasing) event occurs. A threshold of 0 corresponds to the maximum value within an interval. (Such an event is designated as -1).
#'
#'
#' @return A list with four named objects is returned, which encodes the order of the events and the clusters.
#'
#'
#' @importFrom methods is
#' @importFrom stats p.adjust wilcox.test t.test
#' @importFrom nem transitive.reduction
#' @importFrom igraph graph_from_adjacency_matrix adjacent_vertices V are.connected vertex_connectivity
#'
#' @seealso \code{\link[e1071]{cmeans}} for clustering time profiles, \code{\link{calcEvents}} for identifying events.
#'
#' @export
calculateOrder_combined <- function(list_Tc, list_clusters, list_events, list_timeMap, test="wilcox", fdrSignif=0.05, incrEventTh=0.5, decrEventTh=0.5){


	# stopifnot(is(list_Tc, "list"), length(list_Tc) > 0, (length(clusters) == nrow(Tc)), is(mat_events, "matrix"), (fdrSignif >0 && fdrSignif <= 1 ), (incrEventTh >= 0 && incrEventTh <= 1), (decrEventTh >= 0 && decrEventTh <= 1))

	if (!(test == eventOrderTest$param) && !(test == eventOrderTest$nonParam)){
		stop(paste("Test ", test, " not recognized.", sep=""))
	}

	list_list_distributions = list()
	events_combined = list_events[[1]]

	for (i in 1:length(list_Tc)){
		# 1. split into sub matricies
		list_splitMat = splitIntoSubMatrices(list_Tc[[i]], list_clusters[[i]])

		# 2. getDistOfAllEvents
		list_distributions = getDistOfAllEvents_v3(list_events[[i]], list_splitMat, incrEventTh, decrEventTh, list_timeMap[[i]])

		# 3. Combine all distributions
		list_list_distributions = c(list_list_distributions, list_distributions)

		# 4. Renumber and combine mat_events.
		if (i == 1){
			events_combined[,Col_events$combinedDatasetNum] <- rep(i, nrow(events_combined))
		}
		else {
			# add dataset number
			list_events[[i]][,Col_events$combinedDatasetNum] <-rep(i, nrow(list_events[[i]]))

			# adjust the cluster numbers
			list_events[[i]][,Col_events$clus] <- list_events[[i]][,Col_events$clus] + max(events_combined[,Col_events$clus])

			# add to the events matrix.
			events_combined <- rbind(events_combined, list_events[[i]])
		}
	}



	# 5. Run

	list_pVal_stat <- list()

	if (test==eventOrderTest$nonParam){
		print ("Running wilcoxon test.")
		list_pVal_stat <- performWilcoxonSignedRankTabular(list_list_distributions)


		title_testType = "(Non-parametric)"
	}
	else {
		print("Running t-test.")
		list_pVal_stat <- performTtestsTabular(list_list_distributions)

		title_testType = "(Parametric)"
	}

	# FDR adjustment for multiple testing.
	mat_pValsFDR <- matrix(stats::p.adjust(as.vector(list_pVal_stat[[1]]), method='fdr'), ncol=nrow(events_combined))

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

	list_eventsOrder <- getEventsOrderInGraph(reducted, net) # computing the order of the events in the mat_events based on a bfs search.
	list_eventsOrder <- adjustEvents(list_eventsOrder, signifs)


	mat_events_withOrder <- appendOrder(events_combined, list_eventsOrder)

	res <- list(mat_events_withOrder, list_eventsOrder, signifs, test)
	names(res) <- c("mat_events_withOrder", "individEventOrder", "signifs", "test")
	return (res)


}



#' @title
#' Visualize the order of events.
#'
#' @description
#' A plot depicting the temporal order of events is generated. The order is shown using two visualizations, event maps (top) and event sparklines (bottom). The event sparkline is a summary of the event map.
#'
#' @param theOrder A list containing information regarding the order of events and clusters generated via the \code{calculateOrder} function.
#'
#'
#' @return A plot containing the ordered clusters by occurance of first event is generated (see description).
#'
#'
#' @importFrom methods is
#' @importFrom graphics plot segments axis rect text lines
#' @importFrom shape Arrows
#'
#' @seealso \code{\link{calculateOrder}}
#'
#' @export
visualizeOrder_combined <- function(theOrder){
# Visualise the ordering of events

	stopifnot(is(theOrder, "list"), length(theOrder) == 4)


	mat_events_withOrder <- theOrder$mat_events_withOrder
	list_eventsOrder <- theOrder$individEventOrder
	signifs <- theOrder$signifs
	test <- theOrder$test

	title_testType = ""
	## Getting other things ready for plotting

	if (test==eventOrderTest$nonParam){
		title_testType = "(Non-parametric)"
	}
	else {
		title_testType = "(Parametric)"
	}



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
	# print(list_blocks)


	mat_eventPoints <- getTheClusLines(mat_events_withOrder, list_eventsOrder, signifs, list_blocks, eventStart_x, sigEventsDiff_x, nonSigEventsDiff_x, lineDiff_y)


	eventEnd_x = max(as.numeric(mat_eventPoints[,cols_matFifty$col_x])) #  + blockSpacingFmEvent_x
	blockEnd_x = eventEnd_x + blockSpacingFmEvent_x

	mat_bgRects <- getRectPoints(mat_eventPoints, list_blocks, blockStart_x, blockSpacingFmEvent_x)



	# mat_eventPoints <- adjustByBgRects(mat_bgRects, list_eventsOrder, mat_eventPoints)


	mat_grayLines <- getGrayLines(mat_eventPoints, signifs)





	mat_grayLines_v2 <- getGrayLines_v2(list_eventsOrder, mat_eventPoints, signifs)




	mat_clusConnLines <- getClusConnLines(mat_eventPoints, blockStart_x, blockEnd_x)

	# make all the down clusConnLines to opaque
	mat_clusConnLines <- makeDecrBarOpaque(mat_clusConnLines, TRUE)







	# a set of points to set up the graph.
	# print(mat_eventPoints)
	graphics::plot(mat_eventPoints[,cols_matFifty$col_x], mat_eventPoints[,cols_matFifty$col_y], asp=NA, yaxt="n", lwd=0.25, col=colors_orderedEvents$incr, pch=".", xlab="Temporal order", ylab="Clusters",  main=paste("Temporal order of events in clusters", title_testType, sep=""), xlim=c(eventStart_x, eventEnd_x), ylim=c(-4, max(as.numeric(mat_eventPoints[,cols_matFifty$col_clus]))), xaxt="n", bty="n", cex.main = 0.8) #,
	# print("over here..")
	# background gray rectangles.
	if (nrow(mat_bgRects) > 0){
		for (rowNum in 1:nrow(mat_bgRects)){
			if (rowNum%%2 == 1){
				graphics::rect(mat_bgRects[rowNum,1], lineDiff_y/2, mat_bgRects[rowNum,2], max(as.numeric(mat_eventPoints[,cols_matFifty$col_clus])) + (lineDiff_y/2), col = "#ededed", border=NA)
			}

		}
	}


	graphics::segments(x0=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_x0]), y0=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_y0]), x1=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_x1]), y1=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_y1]), col="#808080", lty="dotted", lwd=0.5)
	# C99999



	# rect(-0.5, 0, 0.5, 17.5, col = "yellow", border="#ededed")

	# the cluster lines
	graphics::segments(x0=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_x0]), y0=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_y0]), x1=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_x1]), y1=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_y1]), col=mat_clusConnLines[,cols_clusPlotObjs$col_col], lwd=16, lend=2)


	# the events (arrows) within a graph

	# print(mat_events_withOrder)
	# print(max(mat_events_withOrder[,Col_events$combinedDatasetNum]))
	# print(mat_eventPoints)
	# return ()



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

	## plotting arrows in event map & event sparkline.

	# printing up and down events (in arc events map)
	y_adjust = 0.15

	for (i in 1:length(Color_multiomics)) {

		mat_upEvents = mat_eventPoints[mat_eventPoints[,cols_matFifty$col_dir] == Color_multiomics[[i]]$incr,]

		mat_downEvents = mat_eventPoints[mat_eventPoints[,cols_matFifty$col_dir] == Color_multiomics[[i]]$decr,]

		### Event map.
		## arrow bottoms
		shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01-0.2,  y1=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01, arr.type="triangle", arr.length=0.01, arr.width=0.01, col=Color_multiomics[[i]]$incr, lwd=2, lend=2, arr.adj=0)

		shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25,  y1=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25-0.2, arr.type="triangle", arr.length=0.01, arr.width=0.01, col=Color_multiomics[[i]]$decr, lwd=2, lend=2, arr.adj=0)


		## arrow bottoms
		shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01-0.2,  y1=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01, arr.type="triangle", arr.length=0.01, arr.width=0.01, col=Color_multiomics[[i]]$incr, lwd=2, lend=2, arr.adj=0)

		shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25,  y1=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25-0.2, arr.type="triangle", arr.length=0.01, arr.width=0.01, col=Color_multiomics[[i]]$decr, lwd=2, lend=2, arr.adj=0)




		## sharp arrow heads
		shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01-0.2,  y1=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01, arr.type="triangle", arr.length=0.13, arr.width=0.15, col=Color_multiomics[[i]]$incr, arr.adj=0, segment=FALSE)

		shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25,  y1=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25-0.2, arr.type="triangle", arr.length=0.13, arr.width=0.15, col=Color_multiomics[[i]]$decr, arr.adj=0, segment=FALSE)



		### Event sparkline

		## arrow bottom lines
		shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=rep(-1.2, nrow(mat_upEvents)), y1=rep(-1, nrow(mat_upEvents)), arr.type="triangle", arr.length=0.01, arr.width=0.01, col=Color_multiomics[[i]]$incr, lwd=2, lend=2, arr.adj=0)

		shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=rep(-0.8, nrow(mat_downEvents)), y1=rep(-1, nrow(mat_downEvents)), arr.type="triangle", arr.length=0.01, arr.width=0.01, col=Color_multiomics[[i]]$decr, lwd=2, lend=2, arr.adj=0)

		## sharp arrow heads
		shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=rep(-1.2, nrow(mat_upEvents)), y1=rep(-1, nrow(mat_upEvents)), arr.type="triangle", arr.length=0.13, arr.width=0.15, col=Color_multiomics[[i]]$incr, arr.adj=0, segment=FALSE)

		shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=rep(-0.8, nrow(mat_downEvents)), y1=rep(-1, nrow(mat_downEvents)), arr.type="triangle", arr.length=0.13, arr.width=0.15, col=Color_multiomics[[i]]$decr, arr.adj=0, segment=FALSE)
	}


	vec_labels <- getYaxisLabels(mat_eventPoints)
	mat_labels <- adjustLabelsByDataset(mat_events_withOrder, vec_labels)


	# axis labels on sides 2 and 4.
	Map(axis, side=2, at=seq(1,nrow(mat_labels)), col.axis=rev(mat_labels[,Col_labels$color]), labels=rev(mat_labels[,Col_labels$label]), lwd=0, las=1)
	graphics::axis(2,at=seq(1,nrow(mat_labels)),labels=FALSE, col=NA)

	Map(axis, side=4, at=seq(1,nrow(mat_labels)), col.axis=rev(mat_labels[,Col_labels$color]), labels=rev(mat_labels[,Col_labels$label]), lwd=0, las=1)
	graphics::axis(4,at=seq(1,nrow(mat_labels)),labels=FALSE, col=NA)



	# graphics::axis(side=4, las=1, at=seq(1,length(vec_labels)), labels=rev(vec_labels), cex=0.02,  col = NA ) #, col.ticks = 1)

	# x label text.
	mat_xLabelsClus <- getXaxisLabels(list_eventsOrder, mat_eventPoints, signifs, list_blocks, eventStart_x, sigEventsDiff_x, nonSigEventsDiff_x, yLabelInit, yLabelSpace);


	graphics::text(x=as.numeric(mat_xLabelsClus[,1]), y=as.numeric(mat_xLabelsClus[,2]), labels=mat_xLabelsClus[,3], offset=0, cex=0.5, col=mat_xLabelsClus[,4])


}
adjustLabelsByDataset <- function(mat_events_withOrder, vec_labels){
	mat_labels = matrix(ncol=2, nrow=length(vec_labels))


	for (i in 1:length(vec_labels)){
		# print(vec_labels[i])
		datasetNums <- mat_events_withOrder[mat_events_withOrder[,Col_events$clus] == vec_labels[i],Col_events$combinedDatasetNum]


		dsNum = datasetNums[1]

		if (dsNum > 1) {
			theRows <- mat_events_withOrder[,Col_events$combinedDatasetNum] == (dsNum - 1)

			maxClusNumOfPrev <- max(mat_events_withOrder[theRows, Col_events$clus])


			adjustedClusNum <- vec_labels[i] - maxClusNumOfPrev

			mat_labels[i, Col_labels$label] = adjustedClusNum
			mat_labels[i, Col_labels$color] = Color_multiomics[[dsNum]]$incr

			# print(paste(vec_labels[i], dsNum, maxClusNumOfPrev, adjustedClusNum))


		}
		else{
			mat_labels[i, Col_labels$label] = vec_labels[i]
			mat_labels[i, Col_labels$color] = Color_multiomics[[dsNum]]$incr

		}

	}

	return(mat_labels)
}

blah <- function(){
	for (i in 1:length(vec_labels)){
		# get data set num of clusNum
		theRows <- mat_events_withOrder[,Col_events$clus] == vec_labels[i]

		# print(theRows)
		dataSetNum <- mat_events_withOrder[theRows,]

		# print(mat_events_withOrder[theRows,])

		if (length(dataSetNum) > 1){
			dataSetNum <- dataSetNum[1]
		}

		print(paste(vec_labels[i], dataSetNum))

		if (dataSetNum > 1){  # get max clus num from prev data set num; then adjust
			idx = mat_events_withOrder[,Col_events$combinedDatasetNum] = dataSetNum - 1

			maxClusNumOfPrevDataset = max(mat_events_withOrder[idx, Col_events$clus])
			print(maxClusNumOfPrevDataset)
			vec_labels[[i]] = vec_labels[[i]] - maxClusNumOfPrevDataset
		}
	}

}



#############################

getDistOfAllEvents_v3 <- function(mat_fiftyPoints, list_matrices, phosTh, dephosTh, timeMap){
	list_distributions <- list()

	for (eventNum in 1:nrow(mat_fiftyPoints)){
		# print(paste("Eventnum ", eventNum, "--------------------------------------------------------------------"))
		startTp <- mat_fiftyPoints[eventNum, cols_matFifty_v2$startTp]
		endTp <- mat_fiftyPoints[eventNum, cols_matFifty_v2$endTp]
		dist <- genDistForEvent_v2(list_matrices[[mat_fiftyPoints[eventNum,cols_matFifty_v2$clus]]], mat_fiftyPoints[eventNum, cols_matFifty_v2$dir], mat_fiftyPoints[eventNum, cols_matFifty_v2$startTp], mat_fiftyPoints[eventNum, cols_matFifty_v2$endTp], phosTh, dephosTh, timeMap)

		list_distributions[[length(list_distributions) + 1]] <-  dist
	}

	return(list_distributions)
}

genDistForEvent_v2 <- function(clusMat, dir, startTp, endTp, phosTh, dephosTh, timeMap){
	mat_distribution <- matrix(ncol=4)


	for (profNum in 1:nrow(clusMat)){

		if ((dir == 1 && clusMat[profNum, startTp] < clusMat[profNum, endTp]) || (dir == -1 && clusMat[profNum, startTp] > clusMat[profNum, endTp])){
			# direction holds; do calculation, and save.
			# print(clusMat[profNum, startTp:endTp])

			crossPt <- calcCrossing_v4(clusMat[profNum, startTp:endTp], dir, (startTp -1), phosTh, dephosTh, timeMap)

			mat_distribution <- addToMatrix(mat_distribution, profNum, crossPt[1], crossPt[2], dir)

		}
		else{
			mat_distribution <- addToMatrix(mat_distribution, profNum, NA, NA, dir)
		}
	}

	mat_distribution <- mat_distribution[-1,, drop=FALSE]
	return (mat_distribution)
}



calcCrossing_v4 <- function(region, dir, offset, phosTh, dephosTh, timeMap){

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
			# print(paste("Orig", (offset+j-1), " ", offset+j))
			x_50 <- getMidX_v2(timeMap[(offset+j-1)], region[j-1], timeMap[offset+j], region[j], y_50)
			# mat_fiftyPoints <- addToMatrix(mat_fiftyPoints, clusNum, x_50, y_50, 1)
		}
		else if (dir == -1 && ( (y_50 <= region[j-1] && y_50 >=region[j]) || abs(y_50 - region[j-1]) < e  || abs(y_50 - region[j]) < e ) ){
			# 50 is crossed here in dephos direction
			# getMidX()
			# print(paste("Orig", (offset+j-1), " " , offset+j))
			x_50 <- getMidX_v2(timeMap[(offset+j-1)], region[j-1], timeMap[offset+j], region[j], y_50)
			# mat_fiftyPoints <- addToMatrix(mat_fiftyPoints, clusNum, x_50, y_50, -1)

		}
	}

	return (c(x_50, y_50))
}


getMidX_v2 <- function(x1, y1, x2, y2, y_50){

	m <- (y2 - y1)/(x2 - x1)
	x_50 <- x1 + ((y_50 - y1)/m)

	# print(paste(x1, " ", x2, " ", x_50))
	return (x_50)
}
