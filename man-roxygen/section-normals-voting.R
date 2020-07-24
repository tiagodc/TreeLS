#' @section Radius Estimation Through Normal Vectors Intersection:
#'
#' \code{stemPoints} methods that filter points based on eigen decomposition
#' metrics (knn or voxel) provide a rough estimation of stem segments radii
#' by splitting every stem segment into a local voxel space and counting
#' the number of times that point normals intersect on every voxel (votes).
#' Every stem point then has a radius assigned to it, corresponding to the distance 
#' between the point and the voxel with most votes its normal intersects. The average of
#' all points' radii in a stem segment is the segment's radius.