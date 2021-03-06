% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotClusters.r
\name{plotClusters_withEvents}
\alias{plotClusters_withEvents}
\title{Plots of clusters with events overlaid.}
\usage{
plotClusters_withEvents(
  Tc,
  clusters,
  mat_events,
  plotNumCol = 5,
  centroids = NA
)
}
\arguments{
\item{Tc}{A matrix containing time course data.}

\item{clusters}{A vector with same length as number of profiles. The cluster labels are assumed to be numbers, starting from 1.}

\item{mat_events}{A matrix containing the event information generated by running the \code{calcEvents} function.}

\item{plotNumCol}{Number of cluster-plots to display within a given row.}

\item{centroids}{A matrix containing centroids. If no matrix is provided, the centroid are then compuated, and are a simple average of all the profiles within each cluster. Alternatively, for example, centroids from FCM clustering can be provided where profiles are weighted by the membership score.}
}
\value{
None
}
\description{
Shown overlaid on the cluster plots are events computed for the cluster centroids. The increasing and decreasing events are shown using red and blue, respectively.
}
\seealso{
\code{\link[e1071]{cmeans}} for clustering, \code{\link{calcEvents}}
}
