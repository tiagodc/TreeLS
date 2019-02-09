require(magrittr)
require(lidR)

#' @import magrittr
#' @import lidR
#' @useDynLib TreeLS, .registration = TRUE

preCheck = function(las){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  hasClass = "Classification" %in% names(las@data)

  if( !hasClass ){
    warning('no Classification field found in the dataset')

  }else{

    hasGround = any(las$Classification == 2)

    if(!hasGround){
      warning('ground points not classified')

    }else{

      meanGround = las$Z[ las$Classification == 2 ] %>% mean(na.rm=T) %>% abs

      if(meanGround > 0.2)
        warning("point cloud apparently not normalized")
    }
  }

}

#' Set LAS Header to TLS applications
#' @description Alters the header of an input \code{LAS} object by increasing the precision of XYZ scale factors when suitable
#' @param las \code{LAS} object
#' @param xfac,yfac,zfac scale factor of XYZ axes
#' @return copy of \code{las} with the altered header
#' @export
setHeaderTLS = function(las, xfac = 0.0001, yfac = 0.0001, zfac = 0.0001){

  if(class(las)[1] != "LAS")
    stop("las must be a LAS object")

  if(las@header@PHB$`X scale factor` > xfac)
    las@header@PHB$`X scale factor` = xfac

  if(las@header@PHB$`Y scale factor` > yfac)
    las@header@PHB$`Y scale factor` = yfac

  if(las@header@PHB$`Z scale factor` > xfac)
    las@header@PHB$`Z scale factor` = zfac

  return(las)
}

#' Reset the properties of a \code{LAS} object
#' @description Rebuilds the \code{LAS} object based on its data table only - all header params are recalculated
#' @param las \code{LAS} object
#' @return newly assigned \code{LAS} object
#' @export
resetLAS = function(las){

  if(class(las)[1] != "LAS")
    stop("las must be a LAS object")

  las = las@data %>% LAS %>% setHeaderTLS

  return(las)
}

#' Convert a table-like object to \code{LAS}
#' @description Generates a \code{LAS} object from a \code{matrix}, \code{data.frame} or similar structure
#' @param dataMatrix object to be converted to \code{LAS}
#' @param namesVector \code{colnames} of the variables in the \code{data} compartment of the \code{LAS} otput.
#' If \code{NULL} it ignores any names, considering the first 3 columns as XYZ data and nothing else.
#' @return \code{LAS} object
#' @export
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

#' Resets or creates a \code{LAS} object depending on the input's type
#' @description Resets the input whenever it is a \code{LAS} object, or generates a new \code{LAS} from a table-like input
#' @param cloud object to be converted or reset
#' @param ... further arguments passed to \code{\link{toLAS}} for table-like inputs
#' @return \code{LAS} object
#' @export
setLAS = function(cloud, ...){

  if(class(cloud)[1] == 'LAS'){
    cloud %<>% resetLAS
  }else{
    cloud %<>% toLAS(vars, ...)
  }

  return(cloud)
}

#' Wrapper to read point clouds straight to LAS objects suitable for TLS applications
#' @description Reads \emph{las} or \emph{laz} files with \code{\link{readLAS}}, or tries to read other file formats with \code{\link{read.table}}
#' @param file object to be converted or reset
#' @param colNames parameter passed to \code{\link{toLAS}} whenever reading table-like files - default = \code{NULL}
#' @param ... further arguments passed to either \code{readLAS} or \code{read.table}
#' @return \code{LAS} object
#' @export
readTLS = function(file, colNames=NULL, ...){

  format = sub('.+\\.(.+$)', '\\1', file) %>% tolower

  if(format %in% c('laz', 'las')){

    las = readLAS(file, ...) %>% setHeaderTLS

  }else{

    las = read.table(file, ...) %>% toLAS(colNames)

  }

  return(las)
}

#' Extracts XYZ data from \code{LAS}
#' @description Extracts the XYZ slots from a \code{LAS} object into a \code{matrix}
#' @param las \code{LAS} object
#' @return XYZ \code{matrix} object
#' @export
las2xyz = function(las){

  if(class(las)[1] != "LAS")
    stop("las must be a LAS object")

  las = las@data[,1:3] %>% as.matrix
  return(las)
}

#' Samples a point cloud randomly or systematically
#' @description Applies a random sample or voxel thinning algorithm te keep a fraction of the point cloud.
#' @param las \code{LAS} object
#' @param by sampling method: \emph{voxel} for systematic 3D sampling or \emph{random} for random sampling
#' @param val Sampling parameter value. For \code{by = 'voxel'}, \code{val} must be the voxel side length. For \code{bu = 'random'}, it must be the proportion of points to be kept - between 0 and 1.
#' @return \code{LAS} object
#' @export
tlsSample = function(las, by='voxel', val=0.05){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  methods = c('voxel', 'random')

  if(!(by %in% methods))
    stop('choose a valid method: voxel or random')

  if(val <= 0)
    stop('val must be a positive number')

  if(by == methods[1]){

    keep = las@data[,1:3] %>% as.matrix %>% thinCloud(val)

  }else if(by == methods[2]){

    if(val >= 1)
      stop('val must be a number between 0 and 1 for random sampling')

    n = nrow(las@data)
    keep = rbinom(n, 1, val) == 1

  }

  las %<>% lasfilter(keep)
  return(las)

}

#' Normalize a TLS point cloud
#' @description Normalizes a TLS point based on a Digital Terrain Model of the ground points. If the input's ground points are not classified, the \code{\link{csf}} algorithm is applied internally.
#' @param las \code{LAS} object
#' @param res resolution of the DTM used for normalization
#' @param keepGround default = \code{TRUE} - if \code{FALSE}, returns a point cloud with ground points removed
#' @return \code{LAS} object
#' @export
tlsNormalize = function(las, res=.5, keepGround=T){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  if(res <= 0)
    stop('res must be a positive number')

  if(!any(las$Classification == 2)){
    warning('no ground points found, performing ground segmentation')
    las %<>% lasground(csf(class_threshold = 0.2, cloth_resolution = 0.1), last_returns = F)
  }

  grid = las %>% extent %>% raster
  res(grid) = res

  dtm = grid_terrain(las, res = grid, algorithm = knnidw())

  las %<>% lasnormalize(dtm)

  if(!keepGround) las %<>% lasfilter(Classification != 2)

  return(las)

}

#' Map tree occurrences from TLS data
#' @description Estimates tree probability regions from a point cloud based on a Hough Transform circle search
#' @param las \code{LAS} object
#' @param hmin ...
#' @param hmax ...
#' @param hstep ...
#' @param pixel ...
#' @param rad_max ...
#' @param min_den ...
#' @param min_votes ...
#' @return \code{LAS} object
#' @export
treeMap = function(las, hmin = 1, hmax = 3, hstep = 0.5, pixel = 0.025, rad_max = 0.25, min_den = 0.1, min_votes = 3){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  if(hmax <= hmin)
    stop('hmax must be larger than hmin')

  params = list(
    hstep = hstep,
    pixel = pixel,
    rad_max = rad_max,
    min_den = min_den,
    min_votes = min_votes
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

  if(min_den > 1)
    stop('min_den must be between 0 and 1')

  rgz = las$Z %>% range

  if(hmax < rgz[1])
    stop('hmax is too low - below the point cloud')

  if(hmin > rgz[2])
    stop('hmin is too high - above the point cloud')

  preCheck(las)

  if("Classification" %in% names(las@data))
    las %<>% lasfilter(Classification != 2)

  map = stackMap(las %>% las2xyz, hmin, hmax, hstep, pixel, rad_max, min_den, min_votes) %>%
    do.call(what=cbind) %>% as.data.frame

  map$Intensity %<>% as.integer
  map$Keypoint_flag %<>% as.logical
  map$PointSourceID %<>% as.integer
  map$TreePosition %<>% as.logical
  map %<>% LAS %>% setHeaderTLS

  return(map)
}

#' Get tree XY positions from a TLS-tree map
#' @description Estimates tree probability regions from a point cloud based on a Hough Transform circle search
#' @param las \code{LAS} object - \code{treeMap}'s output
#' @return \code{data.frame} of tree IDs and XY coordinates
#' @export
treePositions = function(las){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  if(!('TreePosition' %in% names(las@data)))
    stop('no TreePosition field found in input: check ?treeMap')

  las %<>% lasfilter(TreePosition)

  pos = las@data[,c('TreeID', 'X', 'Y')] %>% as.data.frame

  return(pos)
}

#' Single tree stem point classification
#' @description Classify stem points through the Hough Transform algorithm - it searches for only \emph{one} stem
#' @param las \code{LAS} object
#' @param hstep ...
#' @param max_radius ...
#' @param hbase ...
#' @param pixel_size ...
#' @param min_density ...
#' @param min_votes ...
#' @return \code{LAS} object
#' @export
stemPoints = function(las, hstep=0.5, max_radius=0.25, hbase = c(1,2.5), pixel_size=0.025, min_density=0.1, min_votes=3){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  if(length(hbase) != 2)
    stop('hbase must be a numeric vector of length 2')

  if(diff(hbase) <= 0)
    stop('hbase[2] must be larger than hbase[1]')

  params = list(
    hstep = hstep,
    max_radius = max_radius,
    pixel_size = pixel_size,
    min_density = min_density,
    min_votes = min_votes
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

  if(min_density > 1)
    stop('min_density must be between 0 and 1')

  if(max(las$Z) < 0)
    stop('input Z coordinates are all negative')

  rg = apply(las@data[,1:2], 2, function(x) max(x) - min(x)) %>% as.double

  if(any(rg > 10))
    warning("point cloud doesn't look like a single tree - XY extents are too large")

  if(min(las$Z) < 0)
    warning("points with Z below 0 will be ignored")

  if(min(las$Z) > 5)
    warning("point cloud doesn't look normalized - Z values too high")

  results = houghStemPoints(las %>% las2xyz, hbase[1], hbase[2], hstep, max_radius, pixel_size, min_density, min_votes)

  las@data$Stem = results$Stem
  las@data$Radius = results$Radius
  las@data$Votes = results$Votes

  las %<>% resetLAS

  return(las)

}


