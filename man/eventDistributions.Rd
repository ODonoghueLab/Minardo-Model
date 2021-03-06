% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.r
\name{eventDistributions}
\alias{eventDistributions}
\title{Distribution density of event times.}
\usage{
eventDistributions(
  Tc,
  clusters,
  mat_events,
  phosEventTh,
  dephosEventTh,
  outfile = "eventDistributions.pdf",
  pdfHeight = 35,
  pdfWidth = 10,
  bandwidth = 0.25
)
}
\arguments{
\item{Tc}{A matrix containing time course data.}

\item{clusters}{A vector with same length as number of profiles. The cluster labels are assumed to be numbers, starting from 1.}

\item{mat_events}{A matrix containing the time regions (and events for the centroid) generated by running the \code{calcEvents} function.}

\item{phosEventTh}{Value between 0 and 1, to define the threshold at which a phosphorylation (or increasing) event occurs. Here, a threshold of 0 corresponds to the minimum value within an interval. (Such an event is designated as 1).}

\item{dephosEventTh}{Value between 0 and 1, to define the threshold at which a dephosphorylation (or decreasing) event occurs. A threshold of 0 corresponds to the maximum value within an interval. (Such an event is designated as -1).}

\item{outfile}{Name of the output pdf file.}

\item{pdfHeight}{Height of a pdf page.}

\item{pdfWidth}{Width of a pdf page.}

\item{bandwidth}{Bandwidth of the density function.}
}
\description{
A pdf is saved containing the density plots of every identified event's time distribution (calculated by calculating the time of event occurance for each profile). If these distributions are approximately normally distributed, then the Student's t-test can be used to compare them, otherwise the Mann-Whitney U test should be utilized.
}
\seealso{
\code{\link[e1071]{cmeans}} for clustering time profiles, \code{\link{calcEvents}} for identifying events.
}
