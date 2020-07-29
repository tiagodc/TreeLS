# ===============================================================================
#
# Developers:
#
# Tiago de Conto - tdc.florestal@gmail.com -  https://github.com/tiagodc/
#
# COPYRIGHT: Tiago de Conto, 2019
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

ptmMetricsLog = function(which_metrics){
  metrics_log = AVAILABLE_POINT_METRICS %in% which_metrics

  if(all(!metrics_log)) stop('Please provide at least one valid metric. See ?fastPointMetrics.available')

  metrics_names = POINT_METRICS_NAMES[1:9][metrics_log[1:9]]
  if(metrics_log[10]) metrics_names %<>% c(POINT_METRICS_NAMES[10:12])
  if(metrics_log[11]) metrics_names %<>% c(POINT_METRICS_NAMES[13:21])

  return(list(log = metrics_log, names = metrics_names))

}

ptmStatistics = function(las, knn, which_metrics = ENABLED_POINT_METRICS){

  pick_metrics = ptmMetricsLog(which_metrics)

  kid = knn$nn.idx
  kds = knn$nn.dists
  kds[kid == 0] = 0

  ptm = data.table()

  if(any(pick_metrics$log[1:11])){
    ptm =  fastPointMetricsCpp(las %>% las2xyz, kid, pick_metrics$log) %>% do.call(what = rbind) %>% as.data.table
    colnames(ptm) = pick_metrics$names
  }

  dist_metrics = AVAILABLE_POINT_METRICS[ 12:17 ][ pick_metrics$log[12:17] ]
  if(length(dist_metrics) > 0){
    dtm = cppFastApply(kds[,-1], dist_metrics) %>% do.call(what=rbind) %>% as.data.table
    colnames(dtm) = dist_metrics
    ptm = cbind(ptm, dtm)
  }

  ptm = as.matrix(ptm)
  ptm[is.na(ptm)] = ptm[is.nan(ptm)] = ptm[is.null(ptm)] = ptm[is.infinite(ptm)] = 0
  ptm = as.data.table(ptm)

  return(ptm)
}


#' Point metrics algorithm: Voxel-wise metrics
#' @description This function is meant to be used inside \code{\link{fastPointMetrics}}. It calculates metrics voxel-wise.
#' @param d \code{numeric} - voxel spacing, in point cloud units.
#' @param exact \code{logical} - use exact voxel search? If \code{FALSE}, applies approximate voxel search unisg integer index hashing, much faster on large point clouds (several million points).
#' @export
ptm.voxel = function(d = .1, exact=FALSE){

  if(!is.numeric(d)) stop('d must be a number')
  if(d <= 0) stop('d must be a positive number')
  if(!is.logical(exact)) stop('exact must be logical')

  func = function(las, which_metrics){

    pick_metrics = ptmMetricsLog(which_metrics)

    if(exact){

      df = data.table()
      offset = las@data[1,1:3]
      for(var in c('X', 'Y', 'Z')){
        dst = floor( (las[[var]] - offset[[var]]) / d )
        df %<>% cbind(dst)
      }

      vx = paste(df[[1]], df[[2]], df[[3]], sep='_') %>% as.factor %>% as.integer

    }else{

      vx = voxelIndex(las2xyz(las), d)

      uid = unique(vx)
      vxid = data.table(hash = uid, id = 1:length(uid))

      vx = data.table(hash = vx)
      vx = merge(vx, vxid, by='hash', sort=F)
      vx = vx$id %>% as.double
    }

    idx = split(0:(length(vx)-1), vx)
    vtm = voxelMetrics(las2xyz(las), idx, pick_metrics$log) %>% do.call(what = rbind) %>% as.data.table
    colnames(vtm) = pick_metrics$names

    keep_names = colnames(las@data)[ !( colnames(las@data) %in% colnames(vtm) ) ]
    las@data = las@data[, ..keep_names]
    vtm$VoxelID = names(idx) %>% as.double
    las@data$VoxelID = vx
    las@data = merge(las@data, vtm, by='VoxelID', sort=F)
    las@data = las@data[,-c('VoxelID')]
    las@data$VoxelID = vx

    return(las)
  }

  func %<>% setAttribute('ptm_mtd')
  return(func)
}


#' Point metrics algorithm: Nearest Neighborhood metrics
#' @description This function is meant to be used inside \code{\link{fastPointMetrics}}. It calculates metrics from a point's nearest neighborhood (KNN).
#' @param k \code{numeric} - number of closest points to search per neighborhood.
#' @param r \code{numeric} - limit radius for nearest neighbor search. If \code{r == 0}, no distance limit is applied.
#' @export
ptm.knn = function(k = 20, r = 0){

  if(!is.numeric(k) || !is.numeric(r)) stop('k and r must be numbers')
  if(k < 3) stop('k must be a number larger than 3')
  if(r < 0) stop('r must be 0 or higher')

  func = function(las, which_metrics){

    k = nabor::knn(las %>% las2xyz, k = k+1, radius = r)
    ptm = ptmStatistics(las, k, which_metrics)
    las@data[,colnames(ptm)] = ptm
    return(las)
  }

  func %<>% setAttribute('ptm_mtd')
  return(func)
}
