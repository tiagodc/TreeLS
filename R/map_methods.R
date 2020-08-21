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

#' Tree mapping algorithm: Hough Transform
#' @description This function is meant to be used inside \code{\link{treeMap}}. It applies an adapted version of the Hough Transform for circle search. Mode details are given in the sections below.
#' @template param-min_h-max_h
#' @template param-h_step
#' @template param-pixel-size
#' @template param-max-d
#' @template param-min-density
#' @template param-min-votes
#' @section \code{LAS@data} Special Fields:
#'
#' Each point in the \code{LAS} object output represents a pixel center that is
#' \emph{possibly} also a stem cross-section center.
#'
#' The variables describing each point in the output are:
#'
#' \itemize{
#' \item \code{Intensity}: number of votes received by that point
#' \item \code{PointSourceID}: unique stem segment ID (among all trees)
#' \item \code{Keypoint_flag}: if \code{TRUE}, the point is the most likely circle center
#' of its stem segment (\code{PointSourceID})
#' \item \code{Radii}: approximate radius estimated by that point - always a multiple of the \code{pixel_size}
#' \item \code{TreeID}: unique tree ID of the point
#' \item \code{TreePosition}: if \code{TRUE}, the point represents the tree's position coordinate
#' }
#'
#' @template section-hough-transform
#' @section Tree Selection:
#'
#' An initial tree filter is used to select \emph{probable} trees in the input point
#' cloud. Parallel stacked layers, each one as thick as \code{hstep}, undergo the
#' circle search within the \code{hmin}/\code{hmax} limits. On every layer, pixels
#' above the \code{min_votes} criterion are clustered, forming
#' \emph{probability zones}. \emph{Probability zones} vertically aligned
#' on at least 3/4 of the stacked layers are assigned as \emph{tree occurrence regions}
#' and exported in the output map.
#'
#' @template reference-thesis
#' @template example-tree-map
#' @export
map.hough = function(min_h = 1, max_h = 3, h_step = 0.5, pixel_size = 0.025, max_d = 0.5, min_density = 0.1, min_votes = 3){

  if(max_h <= min_h)
    stop('max_h must be larger than min_h')

  params = list(
    h_step = h_step,
    pixel_size = pixel_size,
    max_d = max_d,
    min_density = min_density,
    min_votes = min_votes
  )

  for(i in names(params)){
    val = params[[i]]

    if(length(val) != 1)
      stop( i %>% paste('must be of length 1') )

    if(!is.numeric(val))
      stop( i %>% paste('must be Numeric') )

    if(val <= 0)
      stop( i %>% paste('must be positive') )
  }

  if(min_density > 1)
    stop('min_density must be between 0 and 1')

  func = function(las){

    rgz = las$Z %>% range

    if(max_h < rgz[1])
      stop('hmax is too low - below the point cloud')

    if(min_h > rgz[2])
      stop('hmin is too high - above the point cloud')

    map = stackMap(las %>% las2xyz, min_h, max_h, h_step, pixel_size, max_d/2, min_density, min_votes) %>%
      do.call(what=cbind) %>% as.data.table

    map$Intensity %<>% as.integer
    map$Keypoint_flag %<>% as.logical
    map$PointSourceID %<>% as.integer
    map$TreePosition %<>% as.logical
    map = suppressMessages(map %>% LAS %>% setHeaderTLS %>% setAttribute('tree_map'))

    return(map)
  }

  func %<>% setAttribute('tls_map_mtd')
  return(func)
}


#' Tree mapping algorithm: KNN point geometry
#' @description This function is meant to be used inside \code{\link{treeMap}}. It applies a KNN filter to select points with specific neighborhood features. For more details on geometry features, check out \code{\link{fastPointMetrics}}.
#' @template param-max-curvature
#' @template param-max-verticality
#' @param max_mean_dist \code{numeric} - maximum mean distance tolerated from a point to its nearest neighbors.
#' @template param-max-d
#' @template param-min_h-max_h
#' @details
#' Point metrics are calculated for every point. Points are then removed depending on their point
#' metrics parameters and clustered to represent individual tree regions.
#' Clusters are defined as a function of the expected maximum diameter. Any fields added to the
#' point cloud are described in \code{\link{fastPointMetrics}}.
#' @template section-eigen-decomposition
#' @export
map.eigen.knn = function(max_curvature = .1, max_verticality = 10, max_mean_dist = .1, max_d = .5, min_h = 1.5, max_h = 3){

  params = list(
    max_curvature = max_curvature,
    max_verticality = max_verticality,
    max_d = max_d,
    min_h = min_h,
    max_h = max_h
  )

  for(i in names(params)){
    val = params[[i]]

    if(length(val) != 1)
      stop( i %>% paste('must be of length 1') )

    if(!is.numeric(val))
      stop( i %>% paste('must be Numeric') )

    if(val <= 0)
      stop( i %>% paste('must be positive') )
  }

  if(max_curvature > 1){
    stop('max_curvature must be a number between 0 and 1')
  }

  if(max_verticality > 180){
    stop('max_verticality must be a number between 0 and 180')
  }

  func = function(las){
    las = filter_poi(las, Classification != 2 & Z > min_h & Z < max_h)

    if(lidR::is.empty(las)){
      stop('no points found in the specified min_h/max_h range')
    }

    mtrlst = c('Curvature', 'Verticality', 'MeanDist')

    check_pt_metrics = mtrlst %>% sapply(function(x) hasField(las, x)) %>% as.logical %>% all
    if(!check_pt_metrics){
      message('Calculating knn fastPointMetrics')
      las = fastPointMetrics(las, ptm.knn(), mtrlst)
    }

    f3d = nrow(las@data) / (area(las) * abs(diff(range(las$Z))))
    n1 = ceiling(f3d * (.1^3) * 3) + 1
    n2 = ceiling(f3d * (.25^3) * 3) + 1

    md = 20/f3d
    if(md > max_mean_dist) max_mean_dist = md

    las = filter_poi(las, Curvature < max_curvature & abs(Verticality - 90) < max_verticality & MeanDist < max_mean_dist)
    if(is.empty(las)) stop('map.eigen.knn parameters too restrictive, try increasing some of them.')
    las %<>% nnFilter(.1, n1)
    if(is.empty(las)) stop('map.eigen.knn parameters too restrictive, try increasing some of them.')
    las %<>% nnFilter(.25, n2)
    if(is.empty(las)) stop('map.eigen.knn parameters too restrictive, try increasing some of them.')

    las@data$TreeID = 0
    maxdst = max_d*2
    i = 1
    while(any(las@data$TreeID == 0)){
      xy = las@data[TreeID == 0,.(X,Y)][1,] %>% as.double
      dst = las@data[,sqrt( (X-xy[1])^2 + (Y-xy[2])^2 )] %>% as.double
      las@data[TreeID == 0 & dst < maxdst]$TreeID = i
      i=i+1
    }

    hn = las@data[,.(H=max(Z) - min(Z), .N), by=TreeID]
    if(any(hn$N > f3d)){
      las = filter_poi(las, hn$N[TreeID] > f3d)
    }

    las %<>% setAttribute('tree_map')
    return(las)
  }

  func %<>% setAttribute('tls_map_mtd')
  return(func)
}


#' Tree mapping algorithm: Voxel geometry
#' @description This function is meant to be used inside \code{\link{treeMap}}. It applies a filter to select points belonging to voxels with specific features. For more details on geometry features, check out \code{\link{fastPointMetrics}}.
#' @template param-max-curvature
#' @template param-max-verticality
#' @param voxel_spacing \code{numeric} - voxel side length, in point cloud units.
#' @template param-max-d
#' @template param-min_h-max_h
#' @details
#' Point metrics are calculated for every voxel. Points are then removed depending on their voxel's metrics
#' metrics parameters and clustered to represent individual tree regions.
#' Clusters are defined as a function of the expected maximum diameter. Any fields added to the
#' point cloud are described in \code{\link{fastPointMetrics}}.
#' @template section-eigen-decomposition
#' @export
map.eigen.voxel = function(max_curvature = .15, max_verticality = 15, voxel_spacing = .1, max_d = .5, min_h = 1.5, max_h = 3){

  params = list(
    max_curvature = max_curvature,
    max_verticality = max_verticality,
    voxel_spacing = voxel_spacing,
    max_d = max_d,
    min_h = min_h,
    max_h = max_h
  )

  for(i in names(params)){
    val = params[[i]]

    if(length(val) != 1)
      stop( i %>% paste('must be of length 1') )

    if(!is.numeric(val))
      stop( i %>% paste('must be Numeric') )

    if(val <= 0)
      stop( i %>% paste('must be positive') )
  }

  if(max_curvature > 1){
    stop('max_curvature must be a number between 0 and 1')
  }

  if(max_verticality > 180){
    stop('max_verticality must be a number between 0 and 180')
  }

  func = function(las){
    las = filter_poi(las, Classification != 2 & Z > min_h & Z < max_h)

    if(lidR::is.empty(las)){
      stop('no points found in the specified min_h/max_h range')
    }

    mtrlst = c('N', 'Curvature', 'Verticality')
    check_pt_metrics = c(mtrlst, 'VoxelID') %>% sapply(function(x) hasField(las, x)) %>% as.logical %>% all
    if(!check_pt_metrics){
      message('Calculating voxel fastPointMetrics')
      las = fastPointMetrics(las, ptm.voxel(voxel_spacing), mtrlst)
    }

    f3d = nrow(las@data) / (area(las) * abs(diff(range(las$Z))))
    n1 = ceiling(f3d * (.1^3) * 3) + 1
    n2 = ceiling(f3d * (.25^3) * 3) + 1

    las = filter_poi(las, N > 3 & Curvature < max_curvature & abs(Verticality - 90) < max_verticality)
    if(is.empty(las)) stop('map.eigen.voxel parameters too restrictive, try increasing some of them.')
    las %<>% nnFilter(.1, n1)
    if(is.empty(las)) stop('map.eigen.voxel parameters too restrictive, try increasing some of them.')
    las %<>% nnFilter(.25, n2)
    if(is.empty(las)) stop('map.eigen.voxel parameters too restrictive, try increasing some of them.')

    las@data$TreeID = 0
    maxdst = max_d*2
    i = 1
    while(any(las@data$TreeID == 0)){
      xy = las@data[TreeID == 0,.(X,Y)][1,] %>% as.double
      dst = las@data[,sqrt( (X-xy[1])^2 + (Y-xy[2])^2 )] %>% as.double
      las@data[TreeID == 0 & dst < maxdst]$TreeID = i
      i=i+1
    }

    hn = las@data[,.(H=max(Z) - min(Z), .N), by=TreeID]

    if(any(hn$N > f3d)){
      las = filter_poi(las, hn$N[TreeID] > f3d)
    }

    las %<>% setAttribute('tree_map')
    return(las)
  }

  func %<>% setAttribute('tls_map_mtd')
  return(func)
}


#' Tree mapping algorithm: pick trees manually
#' @description This function is meant to be used inside \code{\link{treeMap}}. It opens an interactive \code{rgl} plot where the user can specify tree locations by clicking.
#' @param map optional tree map to be manually updated.
#' @template param-min_h-max_h
#' @export
map.pick = function(map = NULL, min_h=1, max_h=5){

  if(min_h >= max_h){
    stop('max_h must be larger than min_h')
  }

  if(!is.null(map)){
    if(hasAttribute(map, 'tree_map_dt')){
      # map = map
    }else if(hasAttribute(map, 'tree_map')){
      map %<>% treeMap.positions(F)
    }else{
      stop('map is not a tree_map object: check ?treeMap')
    }
  }

  func = function(las){

    rgz = las$Z %>% range

    if(max_h < rgz[1])
      stop('hmax is too low - below the point cloud')

    if(min_h > rgz[2])
      stop('hmin is too high - above the point cloud')

    if(!is.null(min_h)){
      las = filter_poi(las, Z > min_h)
    }

    if(!is.null(max_h)){
      las = filter_poi(las, Z < max_h)
    }

    if(lidR::is.empty(las)){
      stop('no points found in the specified min_h/max_h range')
    }

    plot(las, size = 1.5, clear_artifacts=F)

    if(!is.null(map)){
      spheres3d(map$X, map$Y, median(las$Z), .33, color='white')
      text3d(map$X, map$Y, min(las$Z)-.5, map$TreeID, cex=1.5, col='yellow')
    }

    axes3d(col='white')
    pts = las@data %$% identify3d(X, Y, Z, tolerance = 50)

    if(length(pts) == 0) return(map)

    ids = 1:length(pts)
    if(!is.null(map)){
      ids = ids + max(map$TreeID)
    }

    tmap = data.table(ids, las@data$X[pts], las@data$Y[pts]) # %>% toLAS(c('X','Y','TreeID'))
    names(tmap) = c('TreeID', 'X','Y')

    if(!is.null(map)){
      tmap = rbind(map, tmap)
    }

    tmap %<>% setAttribute('tree_map_dt')
    return(tmap)
  }

  func %<>% setAttribute('tls_map_mtd')
  return(func)
}
