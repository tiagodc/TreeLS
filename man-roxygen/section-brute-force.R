#' @section Brute Force Cylinder Fit:
#'
#' The brute force cylinder fit approach estimates the axis rotation
#' angles by brute force combined with 2D ransac circle fit. The coordinates
#' of a point cloud representing a single cylinder are iteratively rotated up
#' to a pre defined threshold, and for every iteration a circle is estimated after
#' rotation is performed. The rotation that minimizes the circle parameters the most
#' is used to describe the axis direction of the cylinder with the circle's radius.
#'
#' The parameters returned by the brute force cylinder fit method are:
#' \itemize{
#'    \item \code{X,Y}: 2D circle center coordinates after rotation
#'    \item \code{Radius}: 3D circle radius, in point cloud units
#'    \item \code{Error}: model circle error from the RANSAC least squares fit, after rotation
#'    \item \code{DX,DY}: absolute rotation angles (in degrees) applied to the X and Y axes, respectively
#'    \item \code{AvgHeight}: average height of the stem segment's points
#'    \item \code{N}: number of points belonging to the stem segment
#'  }
