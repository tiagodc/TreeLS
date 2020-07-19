#' @section Least Squares Cylinder Fit:
#'
#' The cylinder fit methods implemented in \emph{TreeLS} estimate parameters of a 3D
#' cylinder that describe its axis orientation and radius.The algorithm used internally 
#' to optimize those parameters is the
#' \href{https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method}{Nelser-Mead simplex}, 
#' which takes as objective function the model describing the point to model surface distance of a 
#' regular 3D cylinder point cloud. The objective function and its parameters are described in Liang et al. (2012)
#' @template reference-liang
