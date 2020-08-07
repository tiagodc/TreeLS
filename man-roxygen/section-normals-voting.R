#' @section Radius Estimation Through Normal Vectors Intersection:
#'
#' \code{stemPoints} methods that filter points based on eigen decomposition
#' metrics (knn or voxel) provide a rough estimation of stem segments radii
#' by splitting every stem segment into a local voxel space and counting
#' the number of times that point normals intersect on every voxel (votes).
#' Every stem point then has a radius assigned to it, corresponding to the distance 
#' between the point and the voxel with most votes its normal intersects. The average of
#' all points' radii in a stem segment is the segment's radius. For approximately straight 
#' vertical stem segments, the voting can be done in 2D (pixels). 
#' 
#' The point normals of this method are extracted from the eigen vectors calculated by 
#' \code{\link{fastPointMetrics}}. On top of the point metrics used for stem point filtering, the following 
#' fields are also added to the \code{LAS} object:
#' 
#' \itemize{
#'    \item \code{Votes}: number of normal vector intersections crossing the point's normal at its estimated center 
#'    \item \code{VotesWeight}: ratio of (votes count) / (highest votes count) per \emph{TreeID}
#'    \item \code{Radius}: estimated stem segment radius
#' }
#' 
#' This method was inspired by the denoising algorithm developed by Olofsson & Holmgren (2016), 
#' but it is not an exact reproduction of their work. 