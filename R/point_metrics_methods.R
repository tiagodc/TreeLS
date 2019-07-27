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
  kds[kid == 0] = NA

  ptm = data.table()

  if(any(pickMetrics$log[1:11])){
    ptm =  pointMetricsCpp(las %>% las2xyz, kid, pickMetrics$log) %>% do.call(what = rbind) %>% as.data.table
    colnames(ptm) = pickMetrics$names
  }

  if(pickMetrics$log[12]) ptm$MeanDistance = kds[,-1] %>% rowMeans(na.rm=T)
  if(pickMetrics$log[13]) ptm$MedianDistance = suppressWarnings(kds[,-1] %>% apply(1, median, na.rm=T))
  if(pickMetrics$log[14]) ptm$MinDistance = suppressWarnings(kds[,-1] %>% apply(1, min, na.rm=T))
  if(pickMetrics$log[15]) ptm$MaxDistance = suppressWarnings(kds[,-1] %>% apply(1, max, na.rm=T))

  ptm = as.matrix(ptm)
  ptm[is.na(ptm)] = ptm[is.nan(ptm)] = ptm[is.null(ptm)] = ptm[is.infinite(ptm)] = 0
  ptm = as.data.table(ptm)

  return(ptm)
}

ptm.voxels = function(d = .1, exact=F){

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

ptm.knn = function(k = 30){

  func = function(las, metrics_list){
    zclass = splitByIndex(las)
    zuq = unique(zclass)

    df = data.table()
    for(i in zuq){
      temp = lasfilter(las, zclass == i)
      knn = RANN::nn2(temp %>% las2xyz, k = k, treetype = 'kd', searchtype = 'standard')
      ptm = ptmStatistics(temp, knn, metrics_list)
      temp@data[,colnames(ptm)] = ptm
      df = rbind(df, temp@data)
    }

    las@data = df
    las = resetLAS(las)
    return(las)
  }

  func %<>% setAttribute('ptm_mtd')
  return(func)
}

ptm.radius = function(r = 0.1, max_k = 30){

  func = function(las, metrics_list){
    zclass = splitByIndex(las)
    zuq = unique(zclass)

    df = data.table()
    for(i in zuq){
      temp = lasfilter(las, zclass == i)
      knn = RANN::nn2(las %>% las2xyz, k = max_k, treetype = 'kd', searchtype = 'radius', radius = r)
      ptm = ptmStatistics(temp, knn, metrics_list)
      temp@data[,colnames(ptm)] = ptm
      df = rbind(df, temp@data)
    }

    las@data = df
    las = resetLAS(las)
    return(las)
  }

  func %<>% setAttribute('ptm_mtd')
  return(func)
}
