#' @section Brute Force Cylinder Fit:
#'
#' The brute force cylinder fit approach estimates the axis rotation
#' angles by brute force combined with 2D ransac circle fit. The coordinates
#' of a point cloud representing a single cylinder are iteratively rotated up
#' to a pre defined threshold, and for every iteration a circle is estimated after 
#' rotation is performed. The rotation that minimizes the circle parameters the most
#' is used to describe the axis direction of the cylinder with the circle's radius. 