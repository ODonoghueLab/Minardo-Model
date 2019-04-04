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
#' @return Figure - containing the ordered clusters by occurance of first event.
#'
#'
#' @importFrom methods is
#' @importFrom stats p.adjust wilcox.test t.test
#' @importFrom nem transitive.reduction
#' @importFrom igraph graph_from_adjacency_matrix adjacent_vertices V are.connected vertex_connectivity
#' @importFrom graphics plot segments axis
#'
#' @seealso \code{\link[e1071]{cmeans}} for clustering time profiles, \code{\link{calc50crossing}} for identifying events.
#'
#' @export
orderTheEvents <- function(Tc, clustered, mat_fiftyPoints, test="wilcox", fdrSignif=0.05 ){
	test_param = "t-test"
	test_nonParam = "wilcox"

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

	}
	else {
		print("Running t-test.")
		list_pVal_stat <- performTtestsTabular(list_distributions)
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
	# print(list_eventsOrder)

	## Getting other things ready for plotting
	# if (FALSE){
	mat_eventPoints <- getTheClusLines(mat_fiftyPoints, list_eventsOrder, signifs)
	mat_grayLines <- getGrayLines(mat_eventPoints, signifs)
	mat_clusConnLines <- getClusConnLines(mat_eventPoints)
	vec_labels <- getYaxisLabels(mat_eventPoints)


	# plotting

	graphics::plot(mat_eventPoints[,cols_matFifty$col_x], mat_eventPoints[,cols_matFifty$col_y], asp=NA, yaxt="n", pch=19, col=mat_eventPoints[,cols_matFifty$col_dir], xlab="Time", ylab="Cluster",  main="Clusters ordered by significant order of occurance of events", xlim=c(0, max(as.numeric(mat_eventPoints[,cols_matFifty$col_x]))), ylim=c(0, max(as.numeric(mat_eventPoints[,cols_matFifty$col_clus]))), xaxt="n")
	graphics::axis(side=2, at=seq(1,length(vec_labels)), labels=rev(vec_labels))


	graphics::segments(x0=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_x0]), y0=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_y0]), x1=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_x1]), y1=as.numeric(mat_clusConnLines[,cols_clusPlotObjs$col_y1]), col=mat_clusConnLines[,cols_clusPlotObjs$col_col])

	graphics::segments(x0=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_x0]), y0=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_y0]), x1=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_x1]), y1=as.numeric(mat_grayLines[,cols_clusPlotObjs$col_y1]), col=mat_grayLines[,cols_clusPlotObjs$col_col], lty="dashed")
	# }

}

################################## AUXILLARY

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

getGrayLines <- function(mat_eventPoints, signifs){
	theGrayLines = matrix(ncol=5)


	for (i in 1:nrow(signifs)){
		for (j in i:ncol(signifs)){
			if (i != j){
				if (signifs[j, i] == FALSE){
					# draw a gray50 line between the i and j events
					rowVal <-  c(mat_eventPoints[i, cols_matFifty$col_x], mat_eventPoints[i, cols_matFifty$col_y], mat_eventPoints[j, cols_matFifty$col_x], mat_eventPoints[j, cols_matFifty$col_y], "gray50")
					theGrayLines <- rbind(theGrayLines, rowVal)
				}
			}
		}
	}
	theGrayLines <- theGrayLines[-1,]

	return (theGrayLines)
}

getClusConnLines <- function(mat_eventPoints){
	mat_clusConnLines = matrix(ncol=5,nrow=nrow(mat_eventPoints))

	# cols_clusPlotObjs$col_x0 <- 1
	# cols_clusPlotObjs$col_y0 <- 2
	# cols_clusPlotObjs$col_x1 <- 3
	# cols_clusPlotObjs$col_y1 <- 4
	# cols_clusPlotObjs$col_col <- 5


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

	}

	return (mat_clusConnLines)

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
				mat_eventPoints[eventNum, cols_matFifty$col_dir] = "red"
			}
			else{
				mat_eventPoints[eventNum, cols_matFifty$col_dir] = "blue"
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
				x <- x + 0.5
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
