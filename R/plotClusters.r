#' Prints a plot of clusters using a monochromatic low opacity gray (gray50) colour.
#'
#' Using this colour scheme enables you to clearly visualise the density of profiles within each cluster.
#'
#' @param Tc A matrix containing time course data.
#' @param clustered A fclust object containing the clustering details
#' @param plotNumCol Number of columns in the multi-cluster plot
#'
#' @return None
#'
#' @export
#'
#' @importFrom graphics axis lines par plot
#' @importFrom grDevices adjustcolor
plotClusters <- function (Tc, clustered, plotNumCol=5){

  stopifnot(is(clustered, "fclust"), is(Tc, "matrix"))


  plotNumRow = ceiling(max(clustered$cluster)/plotNumCol)
  graphics::par(mfrow=c(plotNumRow,plotNumCol))


  allProfMemberships = as.matrix(apply(clustered$membership, 1, max)) # Getting each of the profile's best cluster num.


  for(clusNum in c(1:max(clustered$cluster))){

    profilesInClus <- Tc[clustered$cluster == clusNum,]
    profMemberships <- as.matrix(allProfMemberships[clustered$cluster == clusNum,])

    # profMemberships.sorted <- sort(profMemberships, index=TRUE)

    graphics::plot(NA, xlim=c(1,ncol(Tc)), ylim=c(min(Tc)-0.2, max(Tc)+0.2), xlab="Timepoints", ylab="Standardised profiles",  main = paste("Cluster ", toString(clusNum), "; size=", nrow(profilesInClus), sep=""))
    graphics::axis(1, at=c(1:ncol(Tc)))

    for (j in 1:nrow(profilesInClus)){
        # idx = profMemberships.sorted$ix[j]

        graphics::lines(profilesInClus[j,], col=grDevices::adjustcolor("gray50", alpha.f=0.2))
    }

    graphics::lines(clustered$center[clusNum,], col="black", lwd=1.5)

  }

}



#' Prints a plot of clusters using a monochromatic low opacity gray50 colour.
#'
#' @param Tc A matrix containing time course data.
#' @param clustered A fclust object containing the clustering details
#' @param mat_fiftyPoints A matrix containing the event information generated by running the "calc50crossing" function.
#' @param plotNumCol Number of columns in the multi-cluster plot
#'
#' @return None
#'
#' @seealso \code{\link{plotClusters}}
#'
#' @export
#'
#' @importFrom graphics axis lines par plot points
#' @importFrom grDevices adjustcolor
plotClusters_fifty <- function (Tc, clustered, mat_fiftyPoints, plotNumCol=5){

  stopifnot(is(clustered, "fclust"), is(Tc, "matrix"))


  plotNumRow = ceiling(max(clustered$cluster)/plotNumCol)
  graphics::par(mfrow=c(plotNumRow,plotNumCol))


  allProfMemberships = as.matrix(apply(clustered$membership, 1, max)) # Getting each of the profile's best cluster num.


  for(clusNum in c(1:max(clustered$cluster))){

    profilesInClus <- Tc[clustered$cluster == clusNum,]
    profMemberships <- as.matrix(allProfMemberships[clustered$cluster == clusNum,])


    graphics::plot(NA, xlim=c(1,ncol(Tc)), ylim=c(min(Tc)-0.2, max(Tc)+0.2), xlab="Timepoints", ylab="Standardised profiles",  main = paste("Cluster ", toString(clusNum), "; size=", nrow(profilesInClus), sep=""))
    graphics::axis(1, at=c(1:ncol(Tc)))

    for (j in 1:nrow(profilesInClus)){

        graphics::lines(profilesInClus[j,], col=grDevices::adjustcolor("gray50", alpha.f=0.2))
    }

    graphics::lines(clustered$center[clusNum,], col="black", lwd=1.5)

	isYEncountered = FALSE

	for (rowNum in 1:nrow(mat_fiftyPoints)){
		if (mat_fiftyPoints[rowNum, cols_matFifty$col_clus] == clusNum){
			# plot y line
			if (isYEncountered == FALSE){
			 	isYEncountered = TRUE

				graphics::lines(x=c(1,2,3,4,5,6,7,8,9), y=(rep(mat_fiftyPoints[rowNum, cols_matFifty$col_y], 9)), lty="dashed")

			}

			if (mat_fiftyPoints[rowNum, cols_matFifty$col_dir] == 1){
				color="red"
			}
			else{
				color="#003EFF"
			}

			graphics::points(x=mat_fiftyPoints[rowNum, cols_matFifty$col_x], y=mat_fiftyPoints[rowNum, cols_matFifty$col_y], col=color, pch=19)
		}
	}


  }

}
