#' @section Adapted Hough Transform:
#'
#' The Hough Transform circle search algorithm used in
#' TreeLS applies a constrained circle search on discretized 
#' point cloud layers. Tree-wise, the circle search is  
#' recursive, in which the search for circle parameters 
#' of a stem section is constrained to the 
#' \emph{feature space} of the stem section underneath it.
#' Initial estimates of the stem's \emph{feature space} 
#' are performed on a \emph{baselise} stem segment - i.e.
#' a low height interval where a tree's bole is expected  
#' to be clearly visible in the point cloud.
#' The algorithm is described in detail by Conto et al. (2017).
#'
#' This adapted version of the algorithm is very robust against outliers, 
#' but not against forked or leaning stems.
