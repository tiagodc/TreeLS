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

#' Stem segmentation algorithm: RANSAC circle fit
#' @description This function is meant to be used inside \code{\link{stemSegmentation}}. It applies a least squares circle fit algorithm in a RANSAC fashion over stem segments. More details are given in the sections below.
#' @template param-tol
#' @template param-n-ransac
#' @template param-conf
#' @template param-inliers
#' @template section-ransac
#' @template section-circlefit
#' @template reference-olofsson
#' @template reference-thesis
#' @export
sgt.ransac.circle = function(tol=0.1, n = 10, conf = 0.99, inliers = 0.8){

  params = list(
    tol = tol,
    n = n,
    conf = conf,
    inliers = inliers
  )

  for(i in names(params)){
    val = params[[i]]

    if(!is.numeric(val))
      stop( i %>% paste('must be Numeric') )

    if(length(val) > 1)
      stop( i %>% paste('must be of length 1') )

    if(val <= 0)
      stop( i %>% paste('must be positive') )
  }

  if(n < 3)
    stop('n must be at least 3')

  if(conf >= 1)
    stop('conf must be between 0 and 1')

  if(inliers >= 1)
    stop('inliers must be between 0 and 1')

  if(n > 20)
    message('beware that a large n value increases processing time exponentially')

  func = function(las){

    if(!hasField(las, 'TreeID')){

      message('performing single stem segmentation')

      las = filter_poi(las, Stem == TRUE)

      estimates = ransacStemCircle(las %>% las2xyz, las@data$Segment, las@data$Radius, n, conf, inliers, tol) %>% do.call(what = rbind) %>% as.data.table
      names(estimates) = c('X', 'Y', 'Radius', 'Error', 'Segment')

      z = las@data[, .(AvgHeight = mean(Z), N = .N), Segment]

      estimates %<>% merge(z, by='Segment')

      # Segment = NULL
      estimates = estimates[order(Segment)]

      estimates %<>% setAttribute("single_stem_dt")

    }else{

      message('performing multiple stems segmentation')

      las = filter_poi(las, Stem == TRUE)

      estimates = ransacPlotCircles(las %>% las2xyz, las$TreeID, las$Segment, las$Radius, n, conf, inliers, tol) %>% sapply(do.call, what=rbind) %>% do.call(what = rbind) %>% as.data.table
      names(estimates) = c('X', 'Y', 'Radius', 'Error', 'Segment', 'TreeID')

      z = las@data[, .(AvgHeight = mean(Z), N = .N), .(TreeID, Segment)]

      estimates %<>% merge(z, by=c('TreeID', 'Segment'))

      # TreeID = Segment = NULL
      estimates = estimates[order(TreeID, Segment)]

      estimates %<>% setAttribute("multiple_stems_dt")
    }

    return(estimates)
  }

  func %>% setAttribute('stem_sgmt_mtd') %>% return()

}


#' Stem segmentation algorithm: RANSAC cylinder fit
#' @description This function is meant to be used inside \code{\link{stemSegmentation}}. It applies a least squares cylinder fit algorithm in a RANSAC fashion over stem segments. More details are given in the sections below.
#' @template param-tol
#' @template param-n-ransac
#' @template param-conf
#' @template param-inliers
#' @template section-ransac
#' @template section-cylinderfit
#' @template reference-liang
#' @template reference-olofsson
#' @template reference-thesis
#' @export
sgt.ransac.cylinder = function(tol=0.1, n = 10, conf = 0.95, inliers = 0.9){

  params = list(
    tol = tol,
    n = n,
    conf = conf,
    inliers = inliers
  )

  for(i in names(params)){
    val = params[[i]]

    if(!is.numeric(val))
      stop( i %>% paste('must be Numeric') )

    if(length(val) > 1)
      stop( i %>% paste('must be of length 1') )

    if(val <= 0)
      stop( i %>% paste('must be positive') )
  }

  if(n < 5)
    stop('n must be at least 5')

  if(conf >= 1)
    stop('conf must be between 0 and 1')

  if(inliers >= 1)
    stop('inliers must be between 0 and 1')

  if(n > 15)
    message('beware that a large n value increases processing time exponentially')

  func = function(las){

    if(!hasField(las, 'TreeID')){

      message('performing single stem segmentation')

      las = filter_poi(las, Stem == TRUE)

      estimates = ransacStemCylinder(las %>% las2xyz, las@data$Segment, las@data$Radius, n, conf, inliers, tol) %>% do.call(what = rbind) %>% as.data.table
      names(estimates) = c('rho', 'theta', 'phi', 'alpha', 'Radius', 'Error', 'Segment')

      z = las@data[, .(AvgHeight = mean(Z), N = .N), Segment]
      poses = las@data[,rangeMeans(X,Y,Z),by='Segment']

      estimates %<>% merge(z, by='Segment') %>% merge(poses, by='Segment')

      # Segment = NULL
      estimates = estimates[order(Segment)]

      estimates %<>% setAttribute("single_stem_dt")

    }else{

      message('performing multiple stems segmentation')

      las = filter_poi(las, Stem == TRUE)

      estimates = ransacPlotCylinders(las %>% las2xyz, las$TreeID, las$Segment, las$Radius, n, conf, inliers, tol) %>% lapply(do.call, what=rbind) %>% do.call(what = rbind) %>% as.data.table
      names(estimates) = c('rho', 'theta', 'phi', 'alpha', 'Radius', 'Error', 'Segment', 'TreeID')

      z = las@data[, .(AvgHeight = mean(Z), N = .N), .(TreeID, Segment)]
      poses = las@data[,rangeMeans(X,Y,Z),by=c('TreeID','Segment')]

      estimates %<>% merge(z, by=c('TreeID', 'Segment')) %>% merge(poses, by=c('TreeID', 'Segment'))

      # TreeID = Segment = NULL
      estimates = estimates[order(TreeID, Segment)]

      estimates %<>% setAttribute("multiple_stems_dt")
    }

    return(estimates)
  }

  func %>% setAttribute('stem_sgmt_mtd') %>% return()

}


#' Stem segmentation algorithm: Iterated Reweighted Least Squares circle fit
#' @description This function is meant to be used inside \code{\link{stemSegmentation}}. It applies a reweighted least squares circle fit algorithm using M-estimators in order to remove outlier effects.
#' @template param-tol
#' @param n \code{numeric} - maximum number of points to sample for fitting stem segments.
#' @template section-irls
#' @template section-circlefit
#' @template reference-liang
#' @template reference-thesis
#' @export
sgt.irls.circle = function(tol=0.1, n = 500){

  params = list(
    tol = tol,
    n = n
  )

  for(i in names(params)){
    val = params[[i]]

    if(!is.numeric(val))
      stop( i %>% paste('must be Numeric') )

    if(length(val) > 1)
      stop( i %>% paste('must be of length 1') )

    if(val <= 0)
      stop( i %>% paste('must be positive') )
  }

  if(n < 3)
    stop('n must be at least 3')

  if(tol > 1)
    stop('tol must be smaller than 1')

  if(tol <= 0)
    stop('tol must be larger than 0')

  func = function(las){

    if(!hasField(las, 'TreeID')){

      message('performing single stem segmentation')

      las = filter_poi(las, Stem == TRUE)

      estimates = irlsStemCircle(las %>% las2xyz, las@data$Segment, las@data$Radius, n, tol) %>% do.call(what = rbind) %>% as.data.table
      names(estimates) = c('X', 'Y', 'Radius', 'SSQ', 'Error', 'Segment')

      z = las@data[, .(AvgHeight = mean(Z), N = .N), Segment]

      estimates %<>% merge(z, by='Segment')

      # Segment = NULL
      estimates = estimates[order(Segment)]

      estimates %<>% setAttribute("single_stem_dt")

    }else{

      message('performing multiple stems segmentation')

      las = filter_poi(las, Stem == TRUE)

      estimates = irlsPlotCircles(las %>% las2xyz, las$TreeID, las$Segment, las$Radius, n, tol) %>% sapply(do.call, what=rbind) %>% do.call(what = rbind) %>% as.data.table
      names(estimates) = c('X', 'Y', 'Radius', 'SSQ', 'Error', 'Segment', 'TreeID')

      z = las@data[, .(AvgHeight = mean(Z), N = .N), .(TreeID, Segment)]

      estimates %<>% merge(z, by=c('TreeID', 'Segment'))

      # TreeID = Segment = NULL
      estimates = estimates[order(TreeID, Segment)]

      estimates %<>% setAttribute("multiple_stems_dt")
    }

    estimates$SSQ = NULL
    return(estimates)
  }

  func %>% setAttribute('stem_sgmt_mtd') %>% return()

}


#' Stem segmentation algorithm: Iterated Reweighted Least Squares cylinder fit
#' @description This function is meant to be used inside \code{\link{stemSegmentation}}. It applies a reweighted least squares cylinder fit algorithm using M-estimators and Nelder-Mead optimization in order to remove outlier effects.
#' @template param-tol
#' @param n \code{numeric} - maximum number of points to sample for fitting stem segments.
#' @template section-irls
#' @template section-cylinderfit
#' @template reference-liang
#' @template reference-thesis
#' @export
sgt.irls.cylinder = function(tol=0.1, n = 100){

  params = list(
    tol = tol,
    n = n
  )

  for(i in names(params)){
    val = params[[i]]

    if(!is.numeric(val))
      stop( i %>% paste('must be Numeric') )

    if(length(val) > 1)
      stop( i %>% paste('must be of length 1') )

    if(val <= 0)
      stop( i %>% paste('must be positive') )
  }

  if(n < 5)
    stop('n must be at least 5')

  if(tol > 1)
    stop('tol must be smaller than 1')

  if(tol <= 0)
    stop('tol must be larger than 0')

  func = function(las){

    if(!hasField(las, 'TreeID')){

      message('performing single stem segmentation')

      las = filter_poi(las, Stem == TRUE)

      estimates = irlsStemCylinder(las %>% las2xyz, las@data$Segment, las@data$Radius, n, tol) %>% do.call(what = rbind) %>% as.data.table
      names(estimates) = c('rho', 'theta', 'phi', 'alpha', 'Radius', 'Error', 'Segment')

      z = las@data[, .(AvgHeight = mean(Z), N = .N), Segment]
      poses = las@data[,rangeMeans(X,Y,Z),by='Segment']

      estimates %<>% merge(z, by='Segment') %>% merge(poses, by='Segment')

      # Segment = NULL
      estimates = estimates[order(Segment)]

      estimates %<>% setAttribute("single_stem_dt")

    }else{

      message('performing multiple stems segmentation')

      las = filter_poi(las, Stem == TRUE)

      estimates = irlsPlotCylinders(las %>% las2xyz, las$TreeID, las$Segment, las$Radius, n, tol) %>% sapply(do.call, what=rbind) %>% do.call(what = rbind) %>% as.data.table
      names(estimates) = c('rho', 'theta', 'phi', 'alpha', 'Radius', 'Error', 'Segment', 'TreeID')

      z = las@data[, .(AvgHeight = mean(Z), N = .N), .(TreeID, Segment)]
      poses = las@data[,rangeMeans(X,Y,Z),by=c('TreeID','Segment')]

      estimates %<>% merge(z, by=c('TreeID', 'Segment')) %>% merge(poses, by=c('TreeID', 'Segment'))

      # TreeID = Segment = NULL
      estimates = estimates[order(TreeID, Segment)]

      estimates %<>% setAttribute("multiple_stems_dt")
    }

    return(estimates)
  }

  func %>% setAttribute('stem_sgmt_mtd') %>% return()

}


#' Stem segmentation algorithm: Brute Force cylinder fit
#' @description This function is meant to be used inside \code{\link{stemSegmentation}}. It applies a least squares cylinder fit algorithm in a RANSAC fashion over stem segments. More details are given in the sections below.
#' @template param-tol
#' @template param-n-ransac
#' @template param-conf
#' @template param-inliers
#' @template param-z-dev
#' @template section-brute-force
#' @export
sgt.bf.cylinder = function(tol=0.1, n = 10, conf = 0.95, inliers = 0.9, z_dev = 30){

  params = list(
    tol = tol,
    n = n,
    conf = conf,
    inliers = inliers,
    z_dev = z_dev
  )

  for(i in names(params)){
    val = params[[i]]

    if(!is.numeric(val))
      stop( i %>% paste('must be Numeric') )

    if(length(val) > 1)
      stop( i %>% paste('must be of length 1') )

    if(val <= 0)
      stop( i %>% paste('must be positive') )
  }

  if(n < 5)
    stop('n must be at least 5')

  if(conf >= 1)
    stop('conf must be between 0 and 1')

  if(inliers >= 1)
    stop('inliers must be between 0 and 1')

  if(n > 15)
    message('beware that a large n value increases processing time exponentially')

  func = function(las){

    if(!hasField(las, 'TreeID')){

      message('performing single stem segmentation')

      las = filter_poi(las, Stem == TRUE)

      estimates = bfStemCylinder(las %>% las2xyz, las@data$Segment, las@data$Radius, n, conf, inliers, z_dev, tol) %>%
        do.call(what = rbind) %>% as.data.table
      names(estimates) = c('X', 'Y', 'Radius', 'Error', 'DX', 'DY', 'Segment')

      z = las@data[, .(AvgHeight = mean(Z), N = .N), Segment]

      estimates %<>% merge(z, by='Segment')

      # Segment = NULL
      estimates = estimates[order(Segment)]

      estimates %<>% setAttribute("single_stem_dt")

    }else{

      message('performing multiple stems segmentation')

      las = filter_poi(las, Stem == TRUE)

      estimates = bfPlotCylinders(las %>% las2xyz, las$TreeID, las$Segment, las$Radius, n, conf, inliers, z_dev, tol) %>%
        lapply(do.call, what=rbind) %>%
        do.call(what = rbind) %>% as.data.table
      names(estimates) = c('X', 'Y', 'Radius', 'Error', 'DX', 'DY', 'Segment', 'TreeID')

      z = las@data[, .(AvgHeight = mean(Z), N = .N), .(TreeID, Segment)]

      estimates %<>% merge(z, by=c('TreeID', 'Segment'))

      # TreeID = Segment = NULL
      estimates = estimates[order(TreeID, Segment)]

      estimates %<>% setAttribute("multiple_stems_dt")
    }

    return(estimates)
  }

  func %>% setAttribute('stem_sgmt_mtd') %>% return()

}
