#' @section Iterative Reweighted Least Squares (IRLS) Algorithm:
#'
#' \emph{irls} \code{circle} or \code{cylinder} estimation methods 
#' perform automatic outlier assigning through iterative reweighting
#' with M-estimators, followed by a Nelder-Mead optimization of squared distance sums
#' to determine the best circle/cylinder parameters for a given point
#' cloud. The reweighting strategy used in \emph{TreeLS} is based on 
#' Liang et al. (2012).The Nelder-Mead algorithm implemented in Rcpp was provided by 
#' \href{https://github.com/kthohr/optim}{kthohr/optim}.
#' @template reference-liang
