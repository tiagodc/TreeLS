#' @section Least Squares Cylinder Fit:
#'
#' \loadmathjax
#'
#' The cylinder fit methods implemented in \emph{TreeLS} estimate a 3D
#' cylinder`s axis direction and radius. The algorithm used internally
#' to optimize the cylinder parameters is the
#' \href{https://en.wikipedia.org/wiki/Nelder-Mead_method}{Nelder-Mead simplex},
#' which takes as objective function the model describing the distance from any point
#' to a modelled cylinder`s surface on a regular 3D cylinder point cloud:
#'
#' \mjdeqn{D_{p} = |(p - q) \times a| - r}{Dp = abs((p - q) x a) - r}
#'
#' where:
#'
#' \itemize{
#'    \item \emph{Dp}: distance from a point to the model cylinder`s surface
#'    \item \emph{p}: a point on the cylinder`s surface
#'    \item \emph{q}: a point on the cylinder`s axis
#'    \item \emph{a}: unit vector of cylinder`s direction
#'    \item \emph{r}: cylinder`s radius
#' }
#'
#' The Nelder-Mead algorithm minimizes the sum of squared \emph{Dp} from
#' a set of points belonging to a stem segment - in the context of \emph{TreeLS}.
#'
#' The parameters returned by the cylinder fit methods are:
#' \itemize{
#'    \item \code{rho,theta,phi,alpha}: 3D cylinder estimated axis parameters (Liang et al. 2012)
#'    \item \code{Radius}: 3D cylinder radius, in point cloud units
#'    \item \code{Error}: model cylinder error from the least squares fit
#'    \item \code{AvgHeight}: average height of the stem segment's points
#'    \item \code{N}: number of points belonging to the stem segment
#'    \item \code{PX,PY,PZ}: absolute center positions of the stem segment points, in point cloud units (used for plotting)
#'  }
