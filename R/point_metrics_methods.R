ptmMetricsLog = function(metrics_list){
  metrics_log = point.metrics.check %in% metrics_list

  if(all(!metrics_log)) stop('Please provide at least one known metric. See ?availablePointMetrics')

  metrics_names = point.metrics.names[1:9][metrics_log[1:9]]
  if(metrics_log[10]) metrics_names %<>% c(point.metrics.names[10:12])
  if(metrics_log[11]) metrics_names %<>% c(point.metrics.names[13:21])

  return(list(log = metrics_log, names = metrics_names))

}

ptmStatistics = function(las, knn, metrics_list = point.metrics.check){

  pickMetrics = ptmMetricsLog(metrics_list)

  kid = knn$nn.idx
  kds = knn$nn.dists
  kds[kid == 0] = 0

  ptm = data.table()

  if(any(pickMetrics$log[1:11])){
    ptm =  pointMetricsCpp(las %>% las2xyz, kid, pickMetrics$log) %>% do.call(what = rbind) %>% as.data.table
    colnames(ptm) = pickMetrics$names
  }

  distMetrics = point.metrics.check[ 12:17 ][ pickMetrics$log[12:17] ]
  if(length(distMetrics) > 0){
    dtm = cppFastApply(kds[,-1], distMetrics) %>% do.call(what=rbind) %>% as.data.table
    colnames(dtm) = distMetrics
    ptm = cbind(ptm, dtm)
  }

  ptm = as.matrix(ptm)
  ptm[is.na(ptm)] = ptm[is.nan(ptm)] = ptm[is.null(ptm)] = ptm[is.infinite(ptm)] = 0
  ptm = as.data.table(ptm)

  return(ptm)
}

ptm.voxel = function(d = .1, exact=F){

  func = function(las, metrics_list){

    pickMetrics = ptmMetricsLog(metrics_list)

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
    vtm = voxelMetrics(las2xyz(las), idx, pickMetrics$log) %>% do.call(what = rbind) %>% as.data.table
    colnames(vtm) = pickMetrics$names

    keepNames = colnames(las@data)[ !( colnames(las@data) %in% colnames(vtm) ) ]
    las@data = las@data[, ..keepNames]
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

ptm.knn = function(k = 30, r = 0){

  func = function(las, metrics_list){

    k = nabor::knn(las %>% las2xyz, k = k+1, radius = r)
    ptm = ptmStatistics(las, k, metrics_list)
    las@data[,colnames(ptm)] = ptm
    # zclass = splitByIndex(las)
    # zuq = unique(zclass)

    # df = data.table()
    # for(i in zuq){
    #   temp = lasfilter(las, zclass == i)
    #   knn = nabor::knn(temp %>% las2xyz, k = k+1, radius = r)
    #   ptm = ptmStatistics(temp, knn, metrics_list)
    #   temp@data[,colnames(ptm)] = ptm
    #   df = rbind(df, temp@data)
    # }

    # las@data = df
    # las = resetLAS(las)
    return(las)
  }

  func %<>% setAttribute('ptm_mtd')
  return(func)
}
