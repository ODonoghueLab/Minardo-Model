#' @title
#' Plot clusters using a monochromatic low opacity gray colour.
#'
#' @description
#' Using this colour scheme enables the clear visualisation of the density of profiles within each cluster.
#'
#' @param Tc A matrix containing time course data.
#' @param clustered A fclust object containing the clustering details
#' @param plotNumCol Number of cluster-plots to display within a given row.
#'
#' @return None
#'
#'
#' @importFrom graphics axis lines par plot
#' @importFrom grDevices adjustcolor
#' @importFrom methods is
#'
#'
#' @seealso \code{\link[e1071]{cmeans}} for clustering time profiles and generating the \code{clustered} object.
#'
#' @export
plotClusters <- function (Tc, clustered, plotNumCol=5){

  stopifnot(is(clustered, "fclust"), is(Tc, "matrix"), plotNumCol > 0)


  plotNumRow = ceiling(max(clustered$cluster)/plotNumCol)
  graphics::par(mfrow=c(plotNumRow,plotNumCol))


  allProfMemberships = as.matrix(apply(clustered$membership, 1, max)) # Getting each of the profile's best cluster num.


  for(clusNum in c(1:max(clustered$cluster))){

    profilesInClus <- Tc[clustered$cluster == clusNum,]
    profMemberships <- as.matrix(allProfMemberships[clustered$cluster == clusNum,])

    # profMemberships.sorted <- sort(profMemberships, index=TRUE)
	# print(colnames(Tc))
    graphics::plot(NA, xlim=c(1,ncol(Tc)), xaxt="n", ylim=c(min(Tc)-0.2, max(Tc)+0.2), xlab="Timepoints", ylab="Standardised profiles",  main = paste("Cluster ", toString(clusNum), "; size=", nrow(profilesInClus), sep=""))
	graphics::axis(side=1, at=c(1:ncol(Tc)), labels=colnames(Tc))

    for (j in 1:nrow(profilesInClus)){
        # idx = profMemberships.sorted$ix[j]
		# print(profMemberships[j])
        graphics::lines(profilesInClus[j,], col=grDevices::adjustcolor("gray50", alpha.f=0.2))
    }

    graphics::lines(clustered$center[clusNum,], col="black", lwd=1.5)

  }

}



#### MAIN FUNCTION 3.
#' @title
#' Plots of clusters with events overlaid.
#'
#' @description
#' Shown overlaid on the cluster plots are events computed for the cluster centroids. The increasing and decreasing events are shown using red and blue, respectively.
#'
#' @param Tc A matrix containing time course data.
#' @param clustered A fclust object containing the clustering details
#' @param mat_events A matrix containing the event information generated by running the \code{calcEvents} function.
#' @param plotNumCol Number of cluster-plots to display within a given row.
#'
#' @return None
#'
#' @export
#'
#' @importFrom graphics axis lines par plot points
#' @importFrom grDevices adjustcolor
#' @importFrom methods is
#'
#'
#' @seealso \code{\link[e1071]{cmeans}} for clustering time profiles and generating the \code{clustered} object, \code{\link{calcEvents}}
plotClusters_withEvents <- function (Tc, clustered, mat_events, plotNumCol=5){

  stopifnot(is(clustered, "fclust"), is(Tc, "matrix"), is(mat_events, "matrix"), plotNumCol > 0)


  plotNumRow = ceiling(max(clustered$cluster)/plotNumCol)
  graphics::par(mfrow=c(plotNumRow,plotNumCol))


  allProfMemberships = as.matrix(apply(clustered$membership, 1, max)) # Getting each of the profile's best cluster num.


  for(clusNum in c(1:max(clustered$cluster))){

    profilesInClus <- Tc[clustered$cluster == clusNum,]
    profMemberships <- as.matrix(allProfMemberships[clustered$cluster == clusNum,])


    graphics::plot(NA, xlim=c(1,ncol(Tc)), xaxt="n", ylim=c(min(Tc)-0.2, max(Tc)+0.2), xlab="Timepoints", ylab="Standardised profiles",  main = paste("Cluster ", toString(clusNum), "; size=", nrow(profilesInClus), sep=""))
    graphics::axis(side=1, at=c(1:ncol(Tc)), labels=colnames(Tc))

    for (j in 1:nrow(profilesInClus)){

        graphics::lines(profilesInClus[j,], col=grDevices::adjustcolor("gray50", alpha.f=0.2))
    }

    graphics::lines(clustered$center[clusNum,], col="black", lwd=1.5)

	isYEncountered = FALSE

	for (rowNum in 1:nrow(mat_events)){
		if (mat_events[rowNum, cols_matFifty$col_clus] == clusNum){
			# plot y line
			# if (isYEncountered == FALSE){
			#  	isYEncountered = TRUE

			# 	graphics::lines(x=c(1,2,3,4,5,6,7,8,9), y=(rep(mat_fiftyPoints[rowNum, cols_matFifty$col_y], 9)), lty="dashed")

			#}

			if (mat_events[rowNum, cols_matFifty$col_dir] == 1){
				color="red"
			}
			else{
				color="#003EFF"
			}

			tpStart = mat_events[rowNum, 5]
			tpEnd = mat_events[rowNum, 6]
			graphics::lines(x=seq(tpStart, tpEnd), y=(rep(mat_events[rowNum, cols_matFifty$col_y], tpEnd-tpStart+1)), lty="dashed", col=color)


			graphics::points(x=mat_events[rowNum, cols_matFifty$col_x], y=mat_events[rowNum, cols_matFifty$col_y], col=color, pch=19)
		}
	}


  }

}
