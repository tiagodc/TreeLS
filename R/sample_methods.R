#' Point sampling algorithm: systematic voxel grid
#' @description This function is meant to be used inside \code{\link{tlsSample}}. It selects one random point per voxel at a given spatial resolution.
#' @param spacing \code{numeric} - voxel side length.
#' @examples
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#' print(tls)
#'
#' ### thin point cloud - note the point count
#' tls = tlsSample(tls, voxelize(0.05))
#' print(tls)
#'
#' @export
voxelize = function(spacing = 0.05){

  if(spacing <= 0)
    stop('spacing must be a positive number')

  func = function(las){
    las %>% las2xyz %>% thinCloud(spacing) %>% return()
  }

  func %<>% setAttribute('tls_sample_mtd')

  return(func)
}

#' Point sampling algorithm: random sample
#' @description This function is meant to be used inside \code{\link{tlsSample}}. It selects points randomly, returning a fraction of the input point cloud.
#' @param p \code{numeric} - between 0 and 1 - proportion of points to keep.
#' @examples
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#' print(tls)
#'
#' ### thin point cloud - note the point count
#' tls = tlsSample(tls, randomize(0.33))
#' print(tls)
#'
#' @importFrom stats rbinom
#' @export
randomize = function(p = 0.5){

  if(p <= 0)
    stop('p must be a positive number')

  if(p >= 1)
    stop('p must be a number between 0 and 1 for random sampling')

  func = function(las){
    n = nrow(las@data)
    keep = rbinom(n, 1, p) == 1
    return(keep)
  }

  func %<>% setAttribute('tls_sample_mtd')

  return(func)
}
