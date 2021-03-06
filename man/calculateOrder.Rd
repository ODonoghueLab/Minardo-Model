% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orderedEvents.r, R/orderedEvents_forCombined.r,
%   R/orderedEvents_forFlowChart.r
\name{calculateOrder}
\alias{calculateOrder}
\title{Calculate the ordering of events}
\usage{
calculateOrder(
  Tc,
  clusters,
  mat_events,
  test = "wilcox",
  fdrSignif = 0.05,
  phosEventTh = 0.5,
  dephosEventTh = 0.5
)

calculateOrder(
  Tc,
  clusters,
  mat_events,
  test = "wilcox",
  fdrSignif = 0.05,
  phosEventTh = 0.5,
  dephosEventTh = 0.5
)

calculateOrder(
  Tc,
  clusters,
  mat_events,
  test = "wilcox",
  fdrSignif = 0.05,
  phosEventTh = 0.5,
  dephosEventTh = 0.5
)
}
\arguments{
\item{Tc}{A matrix containing time course data.}

\item{clusters}{A vector with same length as number of profiles. The cluster labels are assumed to be numbers, starting from 1.}

\item{mat_events}{A matrix containing the time regions (and events for the centroid) generated by running the \code{calcEvents} function.}

\item{test}{The test used to compare the event time distributions. The choice is between a parametric or a non-parametric test:
\describe{
        \item{\code{wilcox}}{Choosing this will result in event times being compared using the Mann Whitney U test}
        \item{\code{t-test}}{Choosing this will result in event times being compared using the Students t-test}
}}

\item{fdrSignif}{The significance cutoff for FDR correction, performed after pairwise comparisons of the event times.}

\item{phosEventTh}{Value between 0 and 1, to define the threshold at which a phosphorylation (or increasing) event occurs. Here, a threshold of 0 corresponds to the minimum value within an interval. (Such an event is designated as 1).}

\item{dephosEventTh}{Value between 0 and 1, to define the threshold at which a dephosphorylation (or decreasing) event occurs. A threshold of 0 corresponds to the maximum value within an interval. (Such an event is designated as -1).}
}
\value{
A list with four named objects is returned, which encodes the order of the events and the clusters.

A list with four named objects is returned, which encodes the order of the events and the clusters.

A list with four named objects is returned, which encodes the order of the events and the clusters.
}
\description{
For each event, we generate distributions of time of event occurance using profiles within clusters. These distributions are then pair wise compared and the results of all comparisons are represented in the graph space, where events are the graph nodes, and the directional edges represent that the source node event occurs significantly early than the target node. This graph is then utilized to generate the ordering.

For each event, we generate distributions of time of event occurance using profiles within clusters. These distributions are then pair wise compared and the results of all comparisons are represented in the graph space, where events are the graph nodes, and the directional edges represent that the source node event occurs significantly early than the target node. This graph is then utilized to generate the ordering.

For each event, we generate distributions of time of event occurance using profiles within clusters. These distributions are then pair wise compared and the results of all comparisons are represented in the graph space, where events are the graph nodes, and the directional edges represent that the source node event occurs significantly early than the target node. This graph is then utilized to generate the ordering.
}
\seealso{
\code{\link[e1071]{cmeans}} for clustering time profiles, \code{\link{calcEvents}} for identifying events.

\code{\link[e1071]{cmeans}} for clustering time profiles, \code{\link{calcEvents}} for identifying events.

\code{\link[e1071]{cmeans}} for clustering time profiles, \code{\link{calcEvents}} for identifying events.
}
