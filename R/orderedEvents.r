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
#' @importFrom graphics plot segments axis rect
#' @importFrom shape Arrows
#'
#' @seealso \code{\link[e1071]{cmeans}} for clustering time profiles, \code{\link{calc50crossing}} for identifying events.
#'
#' @export
orderTheEvents <- function(Tc, clustered, mat_fiftyPoints, test="wilcox", fdrSignif=0.05 ){
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

	list_distributions <- getDistOfAllEvents(mat_fiftyPoints, list_matrices)

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

	# print("list new events order")
	# print(list_eventsOrder)

	mat_fiftyPts_withOrder <- appendOrder(mat_fiftyPoints, list_eventsOrder)

	## Getting other things ready for plotting
	# if (FALSE){
	mat_eventPoints <- getTheClusLines(mat_fiftyPoints, list_eventsOrder, signifs)
	mat_grayLines <- getGrayLines(mat_eventPoints, signifs)
	mat_grayLines_v2 <- getGrayLines_v2(list_eventsOrder, mat_eventPoints, signifs)

	mat_clusConnLines <- getClusConnLines(mat_eventPoints)
	vec_labels <- getYaxisLabels(mat_eventPoints)
	mat_bgRects <- getRectPoints(list_eventsOrder, mat_eventPoints, signifs)

	print(mat_eventPoints)
	mat_xLabelsClus <- getXaxisLabels(list_eventsOrder, mat_eventPoints, signifs);


	yPos_eventPts = -1

	# a set of points to set up the graph.
	graphics::plot(mat_eventPoints[,cols_matFifty$col_x], mat_eventPoints[,cols_matFifty$col_y], asp=NA, yaxt="n", lwd=0.25, col=colors_orderedEvents$incr, pch=".", xlab="Temporal order", ylab="Cluster",  main=paste("Clusters ordered by significant order of occurance of events ", title_testType, sep=""), xlim=c(0, max(as.numeric(mat_eventPoints[,cols_matFifty$col_x]))), ylim=c(-4, max(as.numeric(mat_eventPoints[,cols_matFifty$col_clus]))), xaxt="n", bty="n") # pch=19

	# background gray rectangles.
	for (rowNum in 1:nrow(mat_bgRects)){
		if (rowNum%%2 == 1){
			graphics::rect(mat_bgRects[rowNum,1], 0.5, mat_bgRects[rowNum,2], max(as.numeric(mat_eventPoints[,cols_matFifty$col_clus])) + 0.5, col = "#ededed", border=NA)
		}

	}


	# axis labels on sides 2 and 4.
	graphics::axis(side=2, las=1, at=seq(1,length(vec_labels)), labels=rev(vec_labels), cex=0.2,  col = NA ) #, col.ticks = 1)

	graphics::axis(side=4, las=1, at=seq(1,length(vec_labels)), labels=rev(vec_labels), cex=0.2,  col = NA ) #, col.ticks = 1)


	graphics::segments(x0=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_x0]), y0=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_y0]), x1=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_x1]), y1=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_y1]), col="#808080", lty="dotted", lwd=0.5)




	# rect(-0.5, 0, 0.5, 17.5, col = "yellow", border="#ededed")

	# the cluster lines
	graphics::segments(x0=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_x0]), y0=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_y0]), x1=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_x1]), y1=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_y1]), col=mat_clusConnLines[,cols_clusPlotObjs$col_col], lwd=16, lend=2)



	# the event points on the graph.
	# mat_eventPoints[,cols_matFifty$]
	mat_upEvents = mat_eventPoints[mat_fiftyPoints[,cols_matFifty$col_dir] == 1,]
	mat_downEvents = mat_eventPoints[mat_fiftyPoints[,cols_matFifty$col_dir] == -1,]
	# print("now for the event points")
	# print(mat_eventPoints)
	#graphics::points(mat_upEvents[,cols_matFifty$col_x], mat_upEvents[,cols_matFifty$col_y], asp=NA, yaxt="n", lwd=1, pch=24, col=color_events$up, bg=color_events$up, xlab="Temporal order", ylab="Cluster",  main=paste("Clusters ordered by significant order of occurance of events ", title_testType, sep="", cex=0.5), xlim=c(0, max(as.numeric(mat_eventPoints[,cols_matFifty$col_x]))), ylim=c(-1, max(as.numeric(mat_eventPoints[,cols_matFifty$col_clus]))), xaxt="n", bty="n")
	#adjusted <- as.numeric(mat_upEvents[,cols_matFifty$col_y]) - 0.15
	# print(adjusted)
	#graphics::points(mat_upEvents[,cols_matFifty$col_x], adjusted, pch = 15, col=color_events$up, bg=color_events$up)

	shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01-0.2,  y1=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01, arr.type="triangle", arr.length=0.01, arr.width=0.01, col=color_events$up, lwd=2, lend=2, arr.adj=0)

	shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25,  y1=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25-0.2, arr.type="triangle", arr.length=0.01, arr.width=0.01, col=color_events$down, lwd=2, lend=2, arr.adj=0)



	shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01-0.2,  y1=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.01, arr.type="triangle", arr.length=0.13, arr.width=0.15, col=color_events$up, arr.adj=0, segment=FALSE)

	shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25,  y1=as.numeric(mat_downEvents[,cols_matFifty$col_y])+0.25-0.2, arr.type="triangle", arr.length=0.13, arr.width=0.15, col=color_events$down, arr.adj=0, segment=FALSE)


	# graphics::points(mat_downEvents[,cols_matFifty$col_x], mat_downEvents[,cols_matFifty$col_y], asp=NA, yaxt="n", lwd=1, pch=25, col=color_events$down, bg=color_events$down, xlab="Temporal order", ylab="Cluster",  main=paste("Clusters ordered by significant order of occurance of events ", title_testType, sep="", cex=0.5), xlim=c(0, max(as.numeric(mat_eventPoints[,cols_matFifty$col_x]))), ylim=c(-1, max(as.numeric(mat_eventPoints[,cols_matFifty$col_clus]))), xaxt="n", bty="n")
	# adjusted <- as.numeric(mat_downEvents[,cols_matFifty$col_y]) + 0.15
	# graphics::points(mat_downEvents[,cols_matFifty$col_x], adjusted, pch=15, col=color_events$down, bg=color_events$down)

	# gray lines in events points maps
	straightGrayLines = mat_grayLines_v2[mat_grayLines_v2[,cols_grayLines_v2$col_isSemiCirc] == FALSE,]
	# print(straightGrayLines)
	graphics::segments(x0=as.numeric(straightGrayLines[,cols_grayLines_v2$col_x0]), y0=rep(yPos_eventPts, length(straightGrayLines)), x1=as.numeric(straightGrayLines[,cols_grayLines_v2$col_x1]), y1=rep(yPos_eventPts, length(straightGrayLines)), lwd=1, col="#A2A2A2") #, lty="dotted")

	# gray curves in events points map
	mat_curvePoints = mat_grayLines_v2[mat_grayLines_v2[,cols_grayLines_v2$col_isSemiCirc] == TRUE,]
	# print(mat_curvePoints)
	# print(mat_curvePoints[1, cols_grayLines_v2$col_x0])
	# print(mat_curvePoints[1, cols_grayLines_v2$col_x1])
	print("the gray lines v2")
	print(mat_curvePoints)
	print(nrow(mat_curvePoints))

	if (nrow(mat_curvePoints) > 0){
		for (i in 1:nrow(mat_curvePoints)){
			list_pts = calcCurves(mat_curvePoints[i, cols_grayLines_v2$col_x0], mat_curvePoints[i, cols_grayLines_v2$col_x1], -1)
			graphics::lines(list_pts[[1]], list_pts[[2]], col="#A2A2A2")
		}
	}




	# printing up and down events (in arc events map)
	y_adjust = 0.15

	# Arrows(1,3,4,6,lwd=2, arr.type="triangle")

	shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=rep(-1.2, nrow(mat_upEvents)), y1=rep(-1, nrow(mat_upEvents)), arr.type="triangle", arr.length=0.01, arr.width=0.01, col=color_events$up, lwd=2, lend=2, arr.adj=0)

	shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=rep(-0.8, nrow(mat_downEvents)), y1=rep(-1, nrow(mat_downEvents)), arr.type="triangle", arr.length=0.01, arr.width=0.01, col=color_events$down, lwd=2, lend=2, arr.adj=0)


	shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=rep(-1.2, nrow(mat_upEvents)), y1=rep(-1, nrow(mat_upEvents)), arr.type="triangle", arr.length=0.13, arr.width=0.15, col=color_events$up, arr.adj=0, segment=FALSE)

	shape::Arrows(x0=as.numeric(mat_downEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_downEvents[,cols_matFifty$col_x]), y0=rep(-0.8, nrow(mat_downEvents)), y1=rep(-1, nrow(mat_downEvents)), arr.type="triangle", arr.length=0.13, arr.width=0.15, col=color_events$down, arr.adj=0, segment=FALSE)
#
#	shape::Arrows(x0=as.numeric(mat_upEvents[,cols_matFifty$col_x]), x1=as.numeric(mat_upEvents[,cols_matFifty$col_x]), y0=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.05-0.2,  y1=as.numeric(mat_upEvents[,cols_matFifty$col_y])-0.05, arr.type="triangle", arr.length=0.13, arr.width=0.15, col=color_events$up, arr.adj=0, segment=TRUE)

	# 	 rep(yPos_eventPts-0.08, nrow(mat_upEvents)), pch = 15, col=color_events$up, bg=color_events$up, cex=0.5)
	# graphics::points(mat_upEvents[,cols_matFifty$col_x], rep(yPos_eventPts-0.08, nrow(mat_upEvents)), pch = 15, col=color_events$up, bg=color_events$up, cex=0.5)
	# graphics::points(mat_downEvents[,cols_matFifty$col_x], rep(yPos_eventPts+y_adjust, nrow(mat_downEvents)), pch = 15, col=color_events$down, bg=color_events$down, cex=0.5)


	# graphics::points(mat_upEvents[,cols_matFifty$col_x], rep(yPos_eventPts+y_adjust, nrow(mat_upEvents)), pch=24, col=color_events$up, bg=color_events$up)
	# graphics::points(mat_downEvents[,cols_matFifty$col_x], rep(-1-0.06, nrow(mat_downEvents)), pch=25, col=color_events$down, bg=color_events$down)

	# print(mat_xLabelsClus)
	graphics::text(x=as.numeric(mat_xLabelsClus[,1]), y=as.numeric(mat_xLabelsClus[,2]), labels=mat_xLabelsClus[,3], offset=0, cex=0.5)


	# colors_orderedEvents$incr, bg=colors_orderedEvents$incr)


	# polygon(1:9, c(2,1,2,1,NA,2,1,2,1),
    #    col = c("red", "blue"),
    #    border = c("green", "yellow"),
    #    lwd = 3, lty = c("dashed", "solid"))

	# }

	# dev.off()
	return (mat_fiftyPts_withOrder)
}

calcCurves <- function(x0, x1, yCentPos){
	print("curves part")
	print(x0)
	print(x1)
	x_10pts <- seq(x0, x1, length.out=50)

	y_10pts <- topSemiCircleEq(x_10pts, yCentPos)
	# print(x_10pts)

	list_pts <- list(x_10pts, y_10pts)
	return (list_pts)

}

topSemiCircleEq <- function(xPts, yCentPos){
	yPts = c()


	r = abs(xPts[1] - xPts[length(xPts)])/2
	xCentPos = (xPts[1] + xPts[length(xPts)]) /2
	# print(r^2)
	for (xPt in xPts){
		yPt = sqrt(abs(r^2 - (xPt-xCentPos)^2)) + yCentPos + 0.4
		yPts = c(yPts, yPt)
	}
	# print(yPts)
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




performWilcoxonSignedRankTabular <- function(list_distributions){
	mat_pVal = matrix(nrow=length(list_distributions), ncol=length(list_distributions))
	mat_statistic = matrix(nrow=length(list_distributions), ncol=length(list_distributions))

	for (i in 1:length(list_distributions)){
		for (j in 1:length(list_distributions)){

			if (i != j){
				res <- stats::wilcox.test(list_distributions[[i]][,cols_matFifty$col_x], list_distributions[[j]][,cols_matFifty$col_x], paired=FALSE, conf.int=TRUE)
				# return (res)
				# print (res)
				mat_pVal[i, j] <- res$p.value
				mat_statistic[i, j] <- res$estimate
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
getXaxisLabels <- function(list_eventsOrder, mat_eventPoints, signifs){
	print (list_eventsOrder)
	mat_xLabelsClus = matrix(ncol=3) #x0, y0, labelNames
	x0 = -1; x0_incr = 1; y0_decr = 0.5
	for (i in 1:length(list_eventsOrder)){

		if (i > 1 && isAnyNonSignifFromPrevAll(list_eventsOrder, i, signifs) == TRUE){
			 print("closer ")
			 print(list_eventsOrder[[i]])
			x0 = x0 + 0.25;
		}
		else {
			x0 = x0 + x0_incr;
		}

		y0 = -1.8
		for (event in sort(list_eventsOrder[[i]])){
			rowVal = c(x0, y0, mat_eventPoints[event, cols_matFifty$col_clus])
			mat_xLabelsClus <- rbind(mat_xLabelsClus, rowVal)

			y0 = y0 - y0_decr;
		}


	}

	print(mat_xLabelsClus)
	mat_xLabelsClus <- mat_xLabelsClus[-1,]
	return (mat_xLabelsClus)
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

					# print(paste(j_layer, i_layer, j, i))
					# print(i_layer)

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

getRectPoints <- function(list_eventsOrder, mat_eventPoints, signifs){

	rectPoints = matrix(ncol=2) # c(rect_start, rect_end, col) col alternate between white and background gray

	# print("get rect points")
	# print(list_eventsOrder)

	rectRow = 1; rectStart = 1; rectEnd = 2; # rectCol;
	rectPoints[rectRow, rectStart] = -0.5

	for (i in 2:length(list_eventsOrder)){
		isNonSignif = FALSE


		if (isAnyNonSignifFromPrevAll(list_eventsOrder, i, signifs) == FALSE) {
			# print(as.numeric(mat_eventPoints[list_eventsOrder[[i]][1], cols_matFifty$col_x]))
			rectPoints[rectRow, rectEnd] = as.numeric(mat_eventPoints[list_eventsOrder[[i-1]][1], cols_matFifty$col_x]) + 0.5 # x of current event + 0.5;
			rectRow = rectRow + 1
			rectPoints <- rbind(rectPoints, c(0,0))
			rectPoints[rectRow, rectStart]<- rectPoints[rectRow-1, rectEnd]
		}
	}
	# rectPoints <- rectPoints[-1,]

	rectPoints[rectRow, rectEnd] <- max(as.numeric(mat_eventPoints[, cols_matFifty$col_x])) + 0.5

	if (rectPoints[rectRow, rectStart] == rectPoints[rectRow, rectEnd]){
		rectPoints = rectPoints[-rectRow,]
	}
	# print("Rect points")
	# print(rectPoints)

	return (rectPoints)

}


getClusConnLines <- function(mat_eventPoints){
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
		mat_clusConnLines[i, cols_clusPlotObjs$col_x1] <- getTheEndX(mat_eventPoints, i, nrow(mat_eventPoints))


		# adding color
		mat_clusConnLines[i, cols_clusPlotObjs$col_col] <- mat_eventPoints[i, cols_matFifty$col_dir]


		## handling the adding of the basal
		if (prevClus == -1 ){
			# just store everything
			mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_x0] = -1;
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
			mat_basalConnLines[basalCounter, cols_clusPlotObjs$col_x0] = -1;
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




getTheClusLines <- function(mat_fiftyPoints, list_eventsOrder, signifs){


	mat_eventPoints <- matrix(nrow=nrow(mat_fiftyPoints), ncol=4)

	x <- 0
	y <- max(mat_fiftyPoints[,1])

	for (i in 1:length(list_eventsOrder)){
		for (j in 1:length(list_eventsOrder[[i]])){
			eventNum <- list_eventsOrder[[i]][j]

			# adding y coordinate
			prevY <- getYIfClusPres(mat_fiftyPoints[eventNum, cols_matFifty$col_clus], mat_eventPoints)
			if (prevY == -1){
				mat_eventPoints[eventNum, cols_matFifty$col_y] <- y
				y <- y - 1
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
		}

		# adding x coordinate
		if (i == 1){

			mat_eventPoints <- addXForAll(list_eventsOrder[[i]], mat_eventPoints, x)
		}
		else{
			# if (isAnyNonSignifFromPrev(list_eventsOrder[[i]], list_eventsOrder[[i-1]], signifs) == TRUE){ # place events closer together
			if (isAnyNonSignifFromPrevAll(list_eventsOrder, i, signifs) == TRUE){
				x <- x + 0.25
			}
			else {
				x <- x + 1
			}

			mat_eventPoints <- addXForAll(list_eventsOrder[[i]], mat_eventPoints, x)

		}
	}
	return (mat_eventPoints)

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
