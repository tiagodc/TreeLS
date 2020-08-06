# ===============================================================================
#
# Developers:
#
# Tiago de Conto - tdc.florestal@gmail.com -  https://github.com/tiagodc/
#
# COPYRIGHT: Tiago de Conto, 2020
#
# This piece of software is open and free to use, redistribution and modifications
# should be done in accordance to the GNU General Public License >= 3
# Use this software as you wish, but no warranty is provided whatsoever. For any
# comments or questions on TreeLS, please contact the developer (prefereably through my github account)
#
# If publishing any work/study/research that used the tools in TreeLS,
# please don't forget to cite the proper sources!
#
# Enjoy!
#
# ===============================================================================

#' Point sampling algorithm: systematic voxel grid
#' @description This function is meant to be used inside \code{\link{tlsSample}}. It selects one random point per voxel at a given spatial resolution.
#' @param spacing \code{numeric} - voxel side length, in point cloud units.
#' @export
smp.voxelize = function(spacing = 0.05){

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
#' @param p \code{numeric} - sampling probability (from 0 to 1).
#' @importFrom stats rbinom
#' @export
smp.randomize = function(p = 0.5){

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
