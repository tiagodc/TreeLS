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

#' @import data.table
#' @import lidR
#' @import magrittr
#' @import rgl
#' @useDynLib TreeLS, .registration = TRUE

. = Z = Classification = TreePosition = TreeID = Stem = Segment = NULL

tls.marker = 'tlsAttribute'

setAttribute = function(obj, atnm){
  attr(obj, tls.marker) = atnm
  return(obj)
}

hasAttribute = function(obj, atnm){
  tlsatt = attr(obj, tls.marker)
  bool = is.null(tlsatt) || tlsatt != atnm
  return(!bool)
}

plot.cylinder = function(xCenter = 0, yCenter = 0, hBase = 0, hTop = 1, radius = 0.5, col = 'yellow'){

  axis = matrix(c(
    rep(xCenter, 2),
    rep(yCenter, 2),
    seq(hBase, hTop, length.out = 2)
  ), ncol = 3, byrow = F)

  cyl = cylinder3d(axis, radius = radius)

  mesh = shade3d(addNormals(subdivision3d(cyl, depth = 0)), col = col)
  # mesh = shade3d(cyl, col=col)
}

preCheck = function(las){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  hasClass = "Classification" %in% names(las@data)

  if( !hasClass ){
    message('no Classification field found in the dataset')

  }else{

    hasGround = any(las$Classification == 2)

    if(hasGround){

      meanGround = las$Z[ las$Classification == 2 ] %>% mean(na.rm=T) %>% abs

      if(meanGround > 0.2)
        message("point cloud apparently not normalized")
    }
  }

}

resetLAS = function(las){

  if(class(las)[1] != "LAS")
    stop("las must be a LAS object")

  prj4 = las@proj4string
  vlr = las@header@VLR

  las = las@data %>% LAS %>% setHeaderTLS

  las@proj4string = prj4
  las@header@VLR = vlr

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
      message('converting first three columns only (assumed XYZ coordinates)')

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

hasField = function(las, field){
  any(names(las@data) == field) %>% return()
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

#' @importFrom stats runif
tlsCylinder = function(n=10000, h=100, rad=30, dev=0){

  rad = runif(n, rad-dev, rad+dev)

  z=runif(n = n, min = 0, max = h)

  angs = runif(n, 0, 2*pi)
  x = sin(angs)*rad
  y = cos(angs)*rad

  return(cbind(x,y,z) %>% toLAS)
}


#' Reset or create a \code{LAS} object depending on the input's type
#' @description Reset the input's header if it is a \code{LAS} object, or generate a new \code{LAS} from a table-like input. For more information, checkout the \code{\link[lidR:LAS]{lidR::LAS}} description page.
#' @param cloud \code{LAS}, \code{data.frame}, \code{matrix} or similar object to be converted or reset.
#' @template param-colnames
#' @template return-las
#' @examples
#' cld = matrix(runif(300, 0, 10), ncol=3)
#' cld = setTLS(cld)
#' summary(cld)
#' @export
setTLS = function(cloud, colNames=NULL){

  if(class(cloud)[1] == 'LAS'){
    cloud %<>% resetLAS
  }else{
    cloud %<>% toLAS(colNames)
  }

  return(cloud)
}


#' Import a point cloud file into a LAS object
#' @description Wrapper to read point clouds straight to LAS objects suitable for TLS applications. Reads \emph{las} or \emph{laz} files with \code{\link[lidR:readLAS]{readLAS}} and alters the header defaults. Other file formats are (or try to be) read using \code{\link[utils:read.table]{read.table}}.
#' @param file file path.
#' @template param-colnames
#' @param ... further arguments passed to either \code{readLAS} or \code{read.table}.
#' @template return-las
#' @examples
#' file = system.file("extdata", "model_boles.laz", package="TreeLS")
#' tls = readTLS(file)
#' summary(tls)
#' @importFrom utils read.table
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


#' Resample a point cloud
#' @description Applies an algorithm that returns a thinned point cloud.
#' @template param-las
#' @param method point sampling algorithm - currently available: \code{\link{voxelize}} or \code{\link{randomize}}
#' @template return-las
#' @examples
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#' summary(tls)
#'
#' vx = tlsSample(tls, voxelize(0.05))
#' summary(vx)
#'
#' rd = tlsSample(tls, randomize(0.5))
#' summary(rd)
#' @export
tlsSample = function(las, method = voxelize()){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  if(!hasAttribute(method, 'tls_sample_mtd'))
    stop('invalid method: check ?tlsSample')

  las %<>% lasfilter(method(las))

  return(las)
}

#' Point cloud cropping
#' @description Returns a cropped point cloud of all points inside or outside specified boundaries of circle or square shapes.
#' @template param-las
#' @param x,y \code{numeric} -  X and Y center coordinates of the area to be cropped.
#' @param len \code{numeric} -  if \code{circle = TRUE}, \code{len} is the circle's radius, otherwise it is the side length of a square.
#' @param circle \code{logical} -  if \code{TRUE} (default), crops a circle, otherwise a square.
#' @param negative \code{logical} - if \code{TRUE}, returns all points outside the specified circle/square boundaries, otherwise returns all points inside the circle/square (default).
#' @template return-las
#' @examples
#' file = system.file("extdata", "model_boles.laz", package="TreeLS")
#' tls = readTLS(file)
#' plot(tls)
#'
#' tls = tlsCrop(tls, 2, 3, 1.5, TRUE, TRUE)
#' plot(tls)
#'
#' tls = tlsCrop(tls, 15, 10, 3, FALSE, FALSE)
#' plot(tls)
#' @export
tlsCrop = function(las, x, y, len, circle=TRUE, negative=FALSE){

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
#' @description Fast normalization of TLS point clouds based on a Digital Terrain Model (DTM) of the ground points. If the input's ground points are not classified, the \code{\link[lidR:csf]{csf}} algorithm is applied internally.
#' @template param-las
#' @param res \code{numeric} - resolution of the DTM used for normalization.
#' @param keepGround \code{logical} - if \code{TRUE} (default), returns a normalized point cloud with classified ground, otherwise removes the ground points.
#' @template return-las
#' @examples
#' file = system.file("extdata", "model_boles.laz", package="TreeLS")
#' tls = readTLS(file)
#' plot(tls)
#'
#' tls = tlsNormalize(tls, 0.5, FALSE)
#' plot(tls)
#' @importFrom raster raster extent res<-
#' @export
tlsNormalize = function(las, res=.5, keepGround=TRUE){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  if(res <= 0)
    stop('res must be a positive number')

  if(!any(las$Classification == 2)){
    message('no ground points found, performing ground segmentation')
    las %<>% lasground(csf(class_threshold = 0.2, cloth_resolution = 0.1), last_returns = F)
  }

  grid = las %>% extent %>% raster
  res(grid) = res

  dtm = grid_terrain(las, res = grid, algorithm = knnidw())

  las %<>% lasnormalize(dtm, TRUE)

  if(!keepGround) las %<>% lasfilter(Classification != 2)

  return(las)

}


#' Map tree occurrences from TLS data
#' @description Estimates tree occurrence regions from a \strong{normalized} point cloud.
#' @template param-las
#' @param method tree mapping algorithm - currently available: \code{\link{map.hough}}.
#' @return \code{LAS} object with \emph{tree_map} signature.
#' @template example-tree-map
#' @section Output:
#' The output is a \code{LAS} object with extra fields in the \code{data} slot. For more details on the output fields checkout \code{\link{map.hough}}'s help page.
#' @export
treeMap = function(las, method = map.hough()){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  if(!hasAttribute(method, 'tls_map_mtd'))
    stop('invalid method: check ?treeMap')

  preCheck(las)

  if("Classification" %in% names(las@data))
    las %<>% lasfilter(Classification != 2)

  map = method(las) %>% setAttribute('tree_map')

  return(map)
}


#' Get unique tree positions from a \emph{tree_map}
#' @description Extracts the tree XY positions from a \emph{tree_map} \code{LAS} object
#' @param las \code{LAS} object - output from \code{\link{treeMap}}.
#' @param plot \code{logical} - plot the tree map?
#' @return \code{data.table} of tree IDs and XY coordinates with \emph{tree_map_dt} signature.
#' @template example-tree-map
#' @export
treePositions = function(las, plot=T){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  if(!hasAttribute(las, 'tree_map'))
    stop('las is not a tree_map object: check ?treeMap')

  las %<>% lasfilter(TreePosition)

  pos = las@data[,c('TreeID', 'X', 'Y')]

  # TreeID = NULL
  pos = pos[order(TreeID)]

  pos %<>% setAttribute('tree_map_dt')

  if(plot){
    pos %$% plot(Y ~ X, cex=3, pch=20, main='tree map', xlab='X', ylab='Y')
  }

  return(pos)
}


#' Stem points classification
#' @description Classify stem points on a \strong{normalized} point cloud.
#' @template param-las
#' @param map optional - map of tree positions (output from \code{\link{treeMap}} or \code{\link{treePositions}}). If omitted, the algorithm assumes \code{las} is a single tree.
#' @param method stem denoising algorithm - currently available: \code{\link{stem.hough}}.
#' @return \code{LAS} object with \emph{stem_points} signature.
#' @examples
#' ### single tree
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#' tls = stemPoints(tls)
#' plot(tls, color='Stem')
#'
#' ### forest plot
#' file = system.file("extdata", "square.laz", package="TreeLS")
#' tls = readTLS(file)
#'
#' # map the trees on a resampled point cloud so all trees have approximately the same point density
#' thin = tlsSample(tls, voxelize(0.02))
#' map = treeMap(thin, map.hough(min_density = 0.05))
#'
#' tls = stemPoints(tls, map)
#' plot(tls, color='Stem')
#'
#' @export
stemPoints = function(las, map = NULL, method = stem.hough()){

  if(class(las)[1] != 'LAS')
    stop('input las must be a LAS object')

  if(!is.null(map)){

    if(hasAttribute(map, 'tree_map_dt')){
      # map = map
    }else if(hasAttribute(map, 'tree_map')){
      map %<>% treePositions(F)
    }else{
      stop('las is not a tree_map object: check ?treeMap')
    }

    if( map$TreeID %>% duplicated %>% any )
      stop('input map must have unique TreeIDs')

  }else{
    rg = apply(las@data[,1:2], 2, function(x) max(x) - min(x)) %>% as.double

    if(any(rg > 10))
      message("point cloud unlikely a single tree (XY extents too large) - check ?tlsCrop")
  }

  if(max(las$Z) < 0)
    stop('input Z coordinates are all negative')

  if(!hasAttribute(method, 'stem_pts_mtd'))
    stop('invalid method: check ?stemPoints')

  las %<>% method(map)

  return(las)

}


#' Stem segmentation
#' @description Measure stem segments from a classified point cloud.
#' @template param-las
#' @param method stem segmentation algorithm - currently available: \code{\link{sgmt.ransac.circle}}.
#' @return \code{data.table} of stem segments with \emph{stem_dt} signature.
#' @examples
#' ### single tree
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#' tls = stemPoints(tls)
#' df = stemSegmentation(tls)
#' head(df)
#'
#' ### forest plot
#' file = system.file("extdata", "square.laz", package="TreeLS")
#' tls = readTLS(file)
#'
#' # map the trees on a resampled point cloud so all trees have approximately the same point density
#' thin = tlsSample(tls, voxelize(0.02))
#' map = treeMap(thin, map.hough(min_density = 0.05))
#'
#' tls = stemPoints(tls, map)
#' df = stemSegmentation(tls)
#' head(df)
#'
#' @export
stemSegmentation = function(las, method=sgmt.ransac.circle()){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  if(!hasAttribute(method, 'stem_sgmt_mtd'))
    stop('invalid method: check ?stemSegmentation')

  las %>% method %>% return()

}
