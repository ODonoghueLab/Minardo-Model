#' @title
#' Multi-omics data sets and clusters.
#'
#' @description
#' This dataset was taken from Yang \emph{et al.} 2019 Cell Syst ( 8, 427-445.e10). It is a multiomics data set, comprising transciptomics, proteomics and phosphoproteomics measurements at multiple time points. The clustering of these time series data is also made available. The phosphoproteomics and proteomics data were clustered using Mfuzz, and the transcriptomics data were clustered using STEM (The STEM clustering details can be found in Kaur \emph{et al.} 2020 npj Sys. Bio.).
#'
#'
#' @name multiomics
#' @docType data
#' @references Yang \emph{et al.} 2019. Available at: \url{http://doi.org/10.1016/j.cels.2019.03.012}
#' @keywords data
NULL

#'
#' @rdname multiomics
"phosTsData.stand"

#' @rdname multiomics
"protTsData.stand"

#' @rdname multiomics
"rnaSeqTsData.stand"

#' @rdname multiomics
"phosClusVec"

#' @rdname multiomics
"protClusVec"

#' @rdname multiomics
"rnaSeqClusVec"
