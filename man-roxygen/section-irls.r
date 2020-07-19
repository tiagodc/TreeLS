#' @section Iterative Reweighted Least Squares (IRLS):
#'
#' \emph{irls} \code{circle} or \code{cylinder} estimation methods 
#' perform automatic outlier assigning through iterative reweighting
#' with M-estimators, followed by least squares optimization as a mean
#' to determine the best circle/cylinder parameters for a given point
#' cloud. The reweighting strategy used in \emph{TreeLS} is based on 
#' Liang et al. (2012)
#' @template reference-liang
