#' @section Least Squares Circle Fit:
#'
#' The circle fit methods applied in \emph{TreeLS} estimate the circle parameters (its center's XY coordinates and radius) 
#' from a pre-selected (denoised) set of points in a least squares fashion 
#' by applying either \href{https://en.wikipedia.org/wiki/QR_decomposition}{QR decompostion}, used in combination
#' with the RANSAC algorithm, or \href{https://en.wikipedia.org/wiki/Nelder-Mead_method}{Nelder-Mead simplex} 
#' optimization combined the IRLS approach.