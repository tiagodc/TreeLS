require(magrittr)
require(lidR)

setHeaderTLS = function(las, xfac = 0.0001, yfac = 0.0001, zfac = 0.0001){

  if(class(las) != "LAS")
    stop("las must be a LAS object")

  if(las@header@PHB$`X scale factor` > xfac)
    las@header@PHB$`X scale factor` = xfac

  if(las@header@PHB$`Y scale factor` > yfac)
    las@header@PHB$`Y scale factor` = yfac

  if(las@header@PHB$`Z scale factor` > xfac)
    las@header@PHB$`Z scale factor` = zfac

  return(las)
}

resetLAS = function(las){

  if(class(las) != "LAS")
    stop("las must be a LAS object")

  las = las@data %>% LAS %>% setHeaderTLS

  return(las)
}

toLAS = function(dataMatrix, namesVector=NULL){

  if(ncol(dataMatrix) < 3)
    stop('dataMatrix must have at least 3 columns')

  dataMatrix %<>% as.data.frame

  if(!is.null(namesVector)){

    if(length(namesVector) != ncol(dataMatrix))
      stop('dataMatrix must have the same number of columns as there are names in namesVector')

    checkXYZ = c('X', 'Y', 'Z') %in% namesVector

    if(!all(checkXYZ))
      stop('X, Y and Z must be declared explicitly in namesVector (in uppercase)')

    names(dataMatrix) = namesVector

  } else {

    if(ncol(dataMatrix) > 3)
      warning('only the first three columns are being converted (assumed XYZ coordinates)')

    dataMatrix = dataMatrix[,1:3]
    names(dataMatrix) = c('X', 'Y', 'Z')

  }

  dataMatrix %<>% LAS %>% setHeaderTLS
  return(dataMatrix)
}

readTLS = function(file, colNames = NULL, ...){

  format = sub('.+\\.(.+$)', '\\1', file) %>% tolower

  if(format %in% c('laz', 'las')){

    las = readLAS(file, ...)

  }else{

    las = read.table(file, ...) %>% toLAS(colNames)

  }

  return(las)
}

las2xyz = function(las){

  if(class(las) != "LAS")
    stop("las must be a LAS object")

  las = las@data[,1:3] %>% as.matrix
  return(las)
}

treeMap = function(las, hmin = 1, hmax = 3, hstep = 0.5, pixel = 0.025, rad_max = 0.25, min_den = 0.1, min_votes = 3){

  if(hmax <= hmin)
    stop('hmax must be larger than hmin')

  if("Classification" %in% names(las@data))
    las %<>% lasfilter(Classification != 2)

  map = stackMap(las@data %>% as.matrix, hmin, hmax, hstep, pixel, rad_max, min_den, min_votes) %>%
    do.call(what=cbind) %>% as.data.frame

  map$Intensity %<>% as.integer
  map$Keypoint_flag %<>% as.logical
  map$PointSourceID %<>% as.integer
  map %<>% LAS %>% setHeaderTLS

  return(map)
}

