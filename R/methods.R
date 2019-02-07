require(magrittr)
require(lidR)

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
      meanGround = las@data[ las$Classification == 2 , 'Z'] %>% mean %>% abs
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
#' @description Reads \emph{las} or \emph{laz} files with \code{\link{lidR::readLAS}}, or tries to read other file formats with \code{\link{read.table}}
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

