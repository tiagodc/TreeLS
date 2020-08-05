#' @section RANSAC Algorithm:
#'
#' \loadmathjax
#' 
#' The \strong{RAN}dom \strong{SA}mple \strong{C}onsensus algorithm is a method that relies on resampling 
#' a data set as many times as necessary to find a subset comprised of only inliers - e.g. observations
#' belonging to a desired model. The RANSAC algorithm provides a way of estimating the necessary number of
#' iterations necessary to fit a model using inliers only, at least once, as shown in the equation:
#' \mjdeqn{k = log(1 - p) / log(1 - w^{n})}{}
#' where:
#' \itemize{
#' \item \emph{k}: number of iterations
#' \item \emph{p}: confidence level, i.e. desired probability of success 
#' \item \emph{w}: proportion of inliers expected in the \emph{full} dataset
#' \item \emph{n}: number of observations sampled on every iteration
#' }
#'
#' The models reiterated in \emph{TreeLS} usually relate to circle or cylinder 
#' fitting over a set of 3D coordinates, selecting the best possible model through the RANSAC algorithm
#'
#' For more information, checkout \href{https://en.wikipedia.org/wiki/Random_sample_consensus}{this wikipedia page}.
