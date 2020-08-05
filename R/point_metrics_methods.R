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

POINT_METRICS_NAMES = c(
  'N',
  'MinDist',
  'MaxDist',
  'MeanDist',
  'SdDist',
  'Linearity',
  'Planarity',
  'Scattering',
  'Omnivariance',
  'Anisotropy',
  'Eigentropy',
  'EigenSum',
  'Curvature',
  'KnnRadius',
  'KnnDensity',
  'Verticality',
  'ZRange',
  'ZSd',
  'KnnRadius2d',
  'KnnDensity2d',
  'EigenSum2d',
  'EigenRatio2d',
  'EigenValue1',
  'EigenValue2',
  'EigenValue3',
  'EigenVector11',
  'EigenVector21',
  'EigenVector31',
  'EigenVector12',
  'EigenVector22',
  'EigenVector32',
  'EigenVector13',
  'EigenVector23',
  'EigenVector33'
)

ENABLED_POINT_METRICS = new.env()
ENABLED_POINT_METRICS$names = POINT_METRICS_NAMES[c(1,13,16,4,32:34)]

ptmMetricsLog = function(which_metrics){
  metrics_log = POINT_METRICS_NAMES %in% which_metrics
  if(all(!metrics_log)) stop('Please provide at least one valid metric. See ?fastPointMetrics.available')
  metrics_names = POINT_METRICS_NAMES[metrics_log]
  return(list(log = metrics_log, names = metrics_names))
}

ptmStatistics = function(las, knn, which_metrics = ENABLED_POINT_METRICS$names){

  pick_metrics = ptmMetricsLog(which_metrics)

  kid = knn$nn.idx
  # kds = knn$nn.dists
  # kds[kid == 0] = 0

  ptm = pointMetricsCpp(las %>% las2xyz, kid, pick_metrics$log) %>% do.call(what = rbind) %>% as.data.table
  colnames(ptm) = pick_metrics$names
  ptm = cleanFields(ptm, colnames(ptm))

  return(ptm)
}


#' Point metrics algorithm: Voxel metrics
#' @description This function is meant to be used inside \code{\link{fastPointMetrics}}. It calculates metrics per voxel.
#' @param d \code{numeric} - voxel spacing, in point cloud units.
#' @param exact \code{logical} - use exact voxel search? If \code{FALSE}, applies approximate voxel search using integer index hashing, much faster on large point clouds (several million points).
#' @export
ptm.voxel = function(d = .1, exact=FALSE){

  if(!is.numeric(d)) stop('d must be a number')
  if(d <= 0) stop('d must be a positive number')
  if(!is.logical(exact)) stop('exact must be logical')

  func = function(las, which_metrics){

    special_case = any(c('KnnDensity', 'KnnDensity2d') %in% which_metrics) && all(which_metrics != 'N')

    if(special_case){
      which_metrics = c('N', which_metrics)
    }

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
    las = cleanFields(las, pick_metrics$names)

    if(hasField(las,'KnnDensity')){
      las@data$KnnDensity = las$N / (d^3)
    }

    if(hasField(las,'KnnDensity2d')){
      las@data$KnnDensity2d = las$N / (d^2)
    }

    if(special_case){
      las@data$N = NULL
    }

    return(las)
  }

  func %<>% setAttribute('ptm_mtd')
  return(func)
}


#' Point metrics algorithm: K Nearest Neighbors metrics
#' @description This function is meant to be used inside \code{\link{fastPointMetrics}}. It calculates metrics for every point using its nearest neighbors (KNN).
#' @param k \code{numeric} - number of nearest points to search per neighborhood.
#' @param r \code{numeric} - search radius limit. If \code{r == 0}, no distance limit is applied.
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
