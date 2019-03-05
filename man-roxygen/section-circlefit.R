#' @section Least Squares Circle Fit:
#'
#' The circle fit method applied in \emph{TreeLS} estimates the circle parameters from a 
#' pre-selected (denoised) set of points 
#' by \href{https://en.wikipedia.org/wiki/QR_decomposition}{QR decompostion}.
#' The optimization criterion for selecting the best circle parameters among several possible candidates
#' is the least squares method, that selects a set of circle parameters that minimize the sum of squared 
#' distances between the model circle and its originating points.
