# ===============================================================================
#
# Developers:
#
# Tiago de Conto - ti@forlidar.com.br -  https://github.com/tiagodc/
#
# COPYRIGHT: Tiago de Conto, 2019
#
# This piece of software is open and free to use, redistribution and modifications
# should be done in accordance to the GNU General Public License >= 3
#
# Use this software as you wish, but no warranty is provided whatsoever. For any
# comments or questions on TreeLS, please contact the developer (prefereably through my github account)
#
# If publishing any work/study/research that used the tools in TreeLS,
# please don't forget to cite the proper sources!
#
# Enjoy!
#
# ===============================================================================

#' @import Rcpp
#' @import magrittr
#' @import lidR
#' @importFrom utils read.table
#' @importFrom stats rbinom
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

resetLAS = function(las){

  if(class(las)[1] != "LAS")
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

    colnames(dataMatrix) = namesVector

  } else {

    if(ncol(dataMatrix) > 3)
      warning('only the first three columns were converted (assumed XYZ coordinates)')

    dataMatrix = dataMatrix[,1:3]
    colnames(dataMatrix) = c('X', 'Y', 'Z')

  }

  dataMatrix %<>% LAS %>% setHeaderTLS
  return(dataMatrix)
}

las2xyz = function(las){

  if(class(las)[1] != "LAS")
    stop("las must be a LAS object")

  las = las@data[,1:3] %>% as.matrix
  return(las)
}

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

#' Resets or creates a \code{LAS} object depending on the input's type
#' @description Resets the input whenever it is a \code{LAS} object, or generates a new \code{LAS} from a table-like input
#' @param cloud object to be converted or reset
#' @param colNames optional - only used for non-LAS objects. It states the \code{file}'s column names - if not set, only the 3 first columns will be converted to \code{LAS}
#' @return \code{LAS} object
#' @export
setTLS = function(cloud, colNames=NULL){

  if(class(cloud)[1] == 'LAS'){
    cloud %<>% resetLAS
  }else{
    cloud %<>% toLAS(colNames)
  }

  return(cloud)
}

#' Wrapper to read point clouds straight to LAS objects suitable for TLS applications
#' @description Reads \emph{las} or \emph{laz} files with \code{\link{readLAS}} and alters the header defaults, or tries to read other file formats with \code{\link{read.table}}
#' @param file object to be converted or reset
#' @param colNames optional - only used for table-like files. It states the \code{file}'s column names - if not set, only the 3 first columns will be imported as XYZ
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

#' Crops a point cloud using a circle or square
#' @description Returns a cropped point cloud of all points inside or outside specified boundaries.
#' @param las \code{LAS} object
#' @param x,y X and Y center coordinates of the area to be cropped
#' @param len if \code{circle = TRUE}, \code{len} is the circle's radius, otherwise it is the side length of a square
#' @param circle if \code{TRUE}, crops a circle, otherwise a square
#' @param negative if \code{TRUE}, returns all points outside the specified circle/square, otherwise returns all points inside the circle/square
#' @return \code{LAS} object
#' @export
tlsCrop = function(las, x, y, len, circle=T, negative=F){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  if(!is.numeric(x))
    stop('x must be Numeric')

  if(!is.numeric(y))
    stop('y must be Numeric')

  if(length(x) != 1)
    stop('x must of length 1')

  if(length(y) != 1)
    stop('y must of length 1')

  if(!is.numeric(len))
    stop('len must be Numeric')

  if(len <= 0)
    stop('len must be a positive number')

  if(length(len) != 1)
    stop('len must of length 1')

  if(!is.logical(circle))
    stop('circle must be Logical')

  if(length(circle) != 1)
    stop('circle must of length 1')

  if(!is.logical(negative))
    stop('negative must be Logical')

  if(length(negative) != 1)
    stop('negative must of length 1')

  bool = RCropCloud(las %>% las2xyz, x, y, len, circle, negative)
  las %<>% lasfilter(bool)

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
    warning('no ground points found initially, ground segmentation was performed internally')
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
    warning("point cloud unlikely a single tree - XY extents too large")

  if(min(las$Z) < 0)
    warning("points with Z below 0 were be ignored")

  if(min(las$Z) > 5)
    warning("point cloud didn't look normalized - Z values too high")

  results = houghStemPoints(las %>% las2xyz, hbase[1], hbase[2], hstep, max_radius, pixel_size, min_density, min_votes)

  las@data$Stem = results$Stem
  las@data$Radius = results$Radius
  las@data$Votes = results$Votes

  las %<>% resetLAS

  return(las)

}

#' Plot-wise stem point classification
#' @description Classify stem points through the Hough Transform algorithm - it searches for one stem per tree in a point cloud
#' @param las \code{LAS} object
#' @param map map of tree positions - output from \code{\link{treeMap}} or \code{\link{treePositions}}
#' @param hstep ...
#' @param max_radius ...
#' @param hbase ...
#' @param pixel_size ...
#' @param min_density ...
#' @param min_votes ...
#' @return \code{LAS} object
#' @export
stemPoints_plot = function(las, map, hstep=0.5, max_radius=0.25, hbase = c(1,2.5), pixel_size=0.025, min_density=0.1, min_votes=3){

  if(class(las)[1] != 'LAS')
    stop('input las must be a LAS object')

  mapNames = c('X', 'Y', 'TreeID')
  if(class(map)[1] == 'LAS'){

    if( !all(mapNames %in% names(map@data)) ){
      stop('input map must be the output from treeMap or treePositions')
    }

    map %<>% treePositions()

  }else if( !all(mapNames %in% names(map)) ){

    stop('input map must be the output from treeMap or treePositions')

  }

  if( any(map$TreeID %>% duplicated) )
    stop('input map must have unique TreeIDs')

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

  if(min(las$Z) < 0)
    warning("points with Z below 0 were be ignored")

  if(min(las$Z) > 5)
    warning("point cloud didn't look normalized - Z values too high")

  hasGround = F
  if(any(names(las@data) == 'Classification')){
    if(any(las$Classification == 2)){
      hasGround = T
      ground = lasfilterground(las)
      las %<>% lasfilter(Classification != 2)
    }
  }

  results = houghStemPlot(las %>% las2xyz, map %>% as.matrix, hbase[1], hbase[2], hstep, max_radius, pixel_size, min_density, min_votes)

  las@data$TreeID = results$TreeID
  las@data$Stem = results$Stem
  las@data$Radius = results$Radius
  las@data$Votes = results$Votes

  if(hasGround){
    ground@data$TreeID = 0
    ground@data$Stem = F
    ground@data$Radius = 0
    ground@data$Votes = 0

    las = las@data %>% rbind(ground@data) %>% toLAS(las@data %>% names)
  }

  las %<>% resetLAS

  return(las)

}

