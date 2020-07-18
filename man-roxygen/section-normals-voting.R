#' @section Normal Vectors Intersection Voting for Radius Estimation:
#'
#' \code{stemPoints} methods that filter points based on eigen decomposition
#' metrics (knn or voxel) provide a rough estimation of stem segments radii
#' by splitting every stem segment into a local pixel or voxel space, and counting
#' the number of times that point/voxel normals intersect on every pixel/voxel (voting).
#' The radius is then estimated by the the average distance of the pixels or voxels
#' with most votes to their connected points.