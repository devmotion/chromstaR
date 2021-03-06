#' chromstaR plotting functions
#' 
#' This page provides an overview of all \pkg{\link{chromstaR}} plotting functions.
#'
#' Plotting functions that work on \code{\link{uniHMM}} objects:
#' \describe{
#'   \item{\code{\link{plotHistogram}}}{Read count histogram with fitted mixture distributions.}
#'   \item{\code{\link{plotKaryogram}}}{Karyogram with read counts and peak calls.}
#' }
#' Plotting functions that work on \code{\link{multiHMM}} objects:
#' \describe{
#'   \item{\code{\link{heatmapCountCorrelation}}}{Heatmap of read count correlations.}
#'   \item{\code{\link{heatmapTransitionProbs}}}{Heatmap of transition probabilities of the Hidden Markov Model.}
#'   \item{\code{\link{heatmapCombinations}}}{Binary presence/absence pattern of combinatorial states.}
#'   \item{\code{\link{plotExpression}}}{Boxplot of expression values that overlap combinatorial states.}
#' }
#' Plotting functions that work on \code{\link{multiHMM}} and \code{\link{combinedMultiHMM}} objects:
#' \describe{
#'   \item{\code{\link{heatmapCountCorrelation}}}{Heatmap of read count correlations.}
#'   \item{\code{\link{plotEnrichCountHeatmap}}}{Heatmap of read counts around annotation.}
#'   \item{\code{\link{plotEnrichment}}}{Enrichment of combinatorial states around annotation.}
#'   \item{\code{\link{plotFoldEnrichHeatmap}}}{Enrichment of combinatorial states at multiple annotations.}
#'   \item{\code{\link{plotExpression}}}{Boxplot of expression values that overlap combinatorial states.}
#' }
#' Other plotting functions:
#' \describe{
#'   \item{\code{\link{heatmapCombinations}}}{Binary presence/absence pattern of combinatorial states.}
#' }
#' @name plotting
NULL
