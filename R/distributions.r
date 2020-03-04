#' @title
#' Densities of standardized abundances at each time for each cluster.
#'
#' @description
#' A pdf is saved containing the density plots at each time point for each cluster. These distributions should be approximately normally distributed to allow the application of GLMs.
#'
#'
#' @param Tc A matrix containing time course data.
#' @param clusters A vector with same length as number of profiles. The cluster labels are assumed to be numbers, starting from 1.
#' @param outfile Name of the output pdf file.
#' @param pdfHeight Height of a pdf page.
#' @param pdfWidth Width of a pdf page.
#' @param bandwidth Bandwidth of the density function.
#'
#'
#' @importFrom methods is
#' @importFrom stats density dnorm sd
#' @importFrom graphics polygon plot axis
#' @importFrom grDevices pdf dev.off
#'
#'
#'
#' @seealso \code{\link[e1071]{cmeans}} for clustering time profiles.
#'
#' @export
clusTpDistributions <- function(Tc, clusters, outfile="abundanceDistributions.pdf", pdfHeight=10, pdfWidth=10, bandwidth=0.5){

	list_matricies <- splitIntoSubMatrices(Tc, clusters)


	pdf(outfile, height=pdfHeight, width=pdfWidth)
	# par(mfrow = c(1, 1))



	plot(0, 0, ylim=c(0, 1.5), xlim=c(-3,3), yaxs="i", ylab="Density", main="Standardized distributions of standardized abundances at each time point for every cluster", xaxt="n", bty="n", yaxt="n", pch='.', xlab="Standardized abundance")
	axis(side =1 ,  at = seq(floor(min(Tc)), ceiling(max(Tc)), 1)) # x-axis label


	for (clusNum in 1:length(list_matricies)){

		for(tp in 1:ncol(list_matricies[[clusNum]])){
			tpMean <- mean(list_matricies[[clusNum]][,tp])
			tpSd <- sd(list_matricies[[clusNum]][,tp])

			standTpVals <- (list_matricies[[clusNum]][,tp] - tpMean)/tpSd

			d <- stats::density(standTpVals, bw=bandwidth)


			lines(d, col=rgb(0.7, 0.7, 0.7), lwd=0.25, xaxt="n", yaxt="n")
		}
	}

	x <- seq(-3, 3, length=100)
	# x <- seq(-3, 3, by = 0.05)
	y <- dnorm(x, mean = 0, sd = 1)

	lines(x, y, xaxt="n", yaxt="n")

	dev.off()
}





#' @title
#' Distribution density of event times.
#'
#' @description
#' A pdf is saved containing the density plots of every identified event's time distribution (calculated by calculating the time of event occurance for each profile). If these distributions are approximately normally distributed, then the Student's t-test can be used to compare them, otherwise the Mann-Whitney U test should be utilized.
#'
#'
#' @param Tc A matrix containing time course data.
#' @param clusters A vector with same length as number of profiles. The cluster labels are assumed to be numbers, starting from 1.
#' @param mat_events A matrix containing the time regions (and events for the centroid) generated by running the \code{calcEvents} function.
#' @param phosEventTh Value between 0 and 1, to define the threshold at which a phosphorylation (or increasing) event occurs. Here, a threshold of 0 corresponds to the minimum value within an interval. (Such an event is designated as 1).
#' @param dephosEventTh Value between 0 and 1, to define the threshold at which a dephosphorylation (or decreasing) event occurs. A threshold of 0 corresponds to the maximum value within an interval. (Such an event is designated as -1).
#' @param outfile Name of the output pdf file.
#' @param pdfHeight Height of a pdf page.
#' @param pdfWidth Width of a pdf page.
#' @param bandwidth Bandwidth of the density function.
#'
#'
#' @importFrom methods is
#' @importFrom stats density
#' @importFrom graphics legend polygon plot axis title lines par
#' @importFrom shape Arrows
#' @importFrom grDevices pdf dev.off rgb
#'
#'
#'
#' @seealso \code{\link[e1071]{cmeans}} for clustering time profiles, \code{\link{calcEvents}} for identifying events.
#'
#' @export
eventDistributions <- function(Tc, clusters, mat_events, phosEventTh, dephosEventTh, outfile="eventDistributions.pdf", pdfHeight=35, pdfWidth=10, bandwidth=0.25) {

	stopifnot(is(Tc, "matrix"), (length(clusters) == nrow(Tc)), is(mat_events, "matrix"), (phosEventTh >= 0 && phosEventTh <= 1), (dephosEventTh >= 0 && dephosEventTh <= 1))

	list_matricies <- splitIntoSubMatrices(Tc, clusters)

	list_distributions <- getDistOfAllEvents_v2(mat_events, list_matricies, phosEventTh, dephosEventTh)


	pdf(outfile, height=pdfHeight, width=pdfWidth)

	par(mfrow = c(length(list_matricies), 1))

	plotsList = list()
	counter_plotsList = 1
	for (eventNum in 1:nrow(mat_events)){
		addIt = FALSE
		plotCol = color_eventDist$up # phos colour
		if (mat_events[eventNum, 4] == -1){ # if dephos event
			plotCol = color_eventDist$down # dephos colour
		}


		if ((eventNum > 1) && (mat_events[eventNum, cols_matFifty$col_clus] ==  mat_events[eventNum-1, cols_matFifty$col_clus])) {

			isAnyNa <- is.na(list_distributions[[eventNum]][,2])
			if (any(isAnyNa)){
				naRmMat <- list_distributions[[eventNum]][-which(isAnyNa),]
			}
			else{
				naRmMat <- list_distributions[[eventNum]]
			}

			d <- getDensityValues(naRmMat, bandwidth)

			polygon(d$density, col=plotCol, border=plotCol)
			lines(c(d$x_med, d$x_med), c(0,d$y_med), col=plotCol)
			lines(c(d$x_mean, d$x_mean), c(0, d$y_mean), col=plotCol, lty="dashed")

			eventCol = color_events$down
			adjustment = -0.2
			if (mat_events[eventNum, 4] == 1){
				eventCol = color_events$up
				adjustment = 0.2
			}

			plotArrows(d, mat_events, eventNum, eventCol, adjustment);
		}

		else {

			isAnyNa <- is.na(list_distributions[[eventNum]][,2])
			if (any(isAnyNa)){
				naRmMat <- list_distributions[[eventNum]][-which(isAnyNa),]
			}
			else{
				naRmMat <- list_distributions[[eventNum]]
			}
			d <- getDensityValues(naRmMat, bandwidth)

			plot(d$density, col=plotCol, xlim=c(0, ncol(Tc)), ylim=c(0,2), xlab="Time", bty='n', yaxt="n", xaxt="n",  yaxs="i", main="", ylab=paste("Cluster ", mat_events[eventNum, 1]))

			polygon(d$density, col=plotCol, border=plotCol)
			lines(c(d$x_med, d$x_med), c(0,d$y_med), col=plotCol)
			lines(c(d$x_mean, d$x_mean), c(0, d$y_mean), col=plotCol, lty="dashed")

			eventCol = color_events$down
			adjustment = -0.2
			if (mat_events[eventNum, 4] == 1){
				eventCol = color_events$up
				adjustment = 0.2
			}

			plotArrows(d, mat_events, eventNum, eventCol, adjustment);

			if (eventNum == 1){
				 graphics::legend("topright", c("Increasing event", "Decreasing event", "Median density", "Mean density"), col=c(color_eventDist$up, color_eventDist$down, 'black', 'black'), lwd=c(10,10, 1, 1), lty=c(1, 1, 1, 2))
			}
			if (eventNum == nrow(mat_events)){
				axis(side =1 ,  at = seq(1, ncol(Tc), 1), labels= colnames(Tc))
			}
			else{
				axis(1, at=seq(1,ncol(Tc),1), col.ticks = 1, labels=FALSE)
			}
		}

	}

	dev.off()
}

plotArrows <- function(d, mat_events, eventNum, eventCol, adjustment){
	if (mat_events[eventNum, 4] == -1){
		shape::Arrows(x0=d$x_med, x1=d$x_med, y0=d$y_med + 0.4 , y1=d$y_med+adjustment + 0.4 , arr.type="triangle", arr.length=0.01, arr.width=0.01, col=eventCol, lwd=2, lend=2, arr.adj=0)

		## sharp arrow heads
		shape::Arrows(x0=d$x_med, x1=d$x_med, y0=d$y_med + 0.4,  y1=d$y_med+adjustment + 0.4, arr.type="triangle", arr.length=0.13, arr.width=0.15, col=eventCol, arr.adj=0, segment=FALSE)
	}
	else{
		shape::Arrows(x0=d$x_med, x1=d$x_med, y0=d$y_med, y1=d$y_med+adjustment, arr.type="triangle", arr.length=0.01, arr.width=0.01, col=eventCol, lwd=2, lend=2, arr.adj=0)

		## sharp arrow heads
		shape::Arrows(x0=d$x_med, x1=d$x_med, y0=d$y_med,  y1=d$y_med+adjustment, arr.type="triangle", arr.length=0.13, arr.width=0.15, col=eventCol, arr.adj=0, segment=FALSE)

	}
}

getDensityValues <- function(naRmMat, bandwidth){
	d <- stats::density(naRmMat[,2], bw=bandwidth)

	dx <- mean(diff(d$x))
	n <- length(d$y)
	y.unit <- sum(d$y) * dx
	y.cs <- cumsum(d$y)

	x.med <- d$x[i.med <- length(y.cs[2*y.cs <= y.cs[n]])]
	y.med <- d$y[i.med]
	x.mean <- sum(d$y * d$x) * dx
	y.mean <- d$y[length(d$x[d$x < x.mean])]
	d$y[1] = 0
	d$y[length(d$y)] = 0

	return (list(density=d, x_med=x.med, y_med=y.med, x_mean=x.mean, y_mean=y.mean))
}