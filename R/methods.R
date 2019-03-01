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
#' @import rgl

tls.marker = 'tlsAttribute'

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

tlsCylinder = function(n=10000, h=100, rad=30, dev=0){

  rad = runif(n, rad-dev, rad+dev)

  z=runif(n = n, min = 0, max = h)

  angs = runif(n, 0, 2*pi)
  x = sin(angs)*rad
  y = cos(angs)*rad

  return(cbind(x,y,z) %>% toLAS)
}


#' Resets or creates a \code{LAS} object depending on the input's type
#' @description Resets the input's header if it is a \code{LAS} object, or generates a new \code{LAS} from a table-like input.
#' @param cloud \code{LAS}, \code{data.frame}, \code{matrix} or similar object to be converted or reset.
#' @template param-colnames
#' @template return-las
#' @examples
#' cld = runif(300, 0, 10) %>% matrix(ncol=3)
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


#' Wrapper to read point clouds straight to LAS objects suitable for TLS applications
#' @description Reads \emph{las} or \emph{laz} files with \code{\link[lidR:readLAS]{readLAS}} and alters the header defaults, or tries to read other file formats with \code{\link[utils:read.table]{read.table}}.
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


#' Samples a point cloud randomly or systematically
#' @description Applies a random sample or voxel thinning algorithm te keep a fraction of the point cloud.
#' @template param-las
#' @param by \code{character} - sampling method: \emph{"voxel"} for systematic 3D sampling or \emph{"random"} for random sampling.
#' @param val \code{numeric} - sampling parameter value. For \code{by = 'voxel'}, \code{val} must be the voxel side length. For \code{by = 'random'}, it must be the proportion of points to be kept (between 0 and 1).
#' @template return-las
#' @examples
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#' summary(tls)
#'
#' tls = tlsSample(tls, by='voxel', val=0.05)
#' summary(tls)
#'
#' tls = tlsSample(tls, by='random', val=0.5)
#' summary(tls)
#' @importFrom stats rbinom
#' @useDynLib TreeLS, .registration = TRUE
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
#' @template param-las
#' @param x,y \code{numeric} -  X and Y center coordinates of the area to be cropped.
#' @param len \code{numeric} -  if \code{circle = TRUE}, \code{len} is the circle's radius, otherwise it is the side length of a square.
#' @param circle \code{logical} -  if \code{TRUE} (default), crops a circle, otherwise a square.
#' @param negative \code{logical} - if \code{TRUE}, returns all points outside the specified circle/square, otherwise returns all points inside the circle/square (default).
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
#' @useDynLib TreeLS, .registration = TRUE
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
#' @description Normalizes a TLS point based on a Digital Terrain Model of the ground points. If the input's ground points are not classified, the \code{\link[lidR:csf]{csf}} algorithm is applied internally.
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
#' @export
tlsNormalize = function(las, res=.5, keepGround=T){

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

  las %<>% lasnormalize(dtm)

  if(!keepGround) las %<>% lasfilter(Classification != 2)

  return(las)

}


#' Map tree occurrences from TLS data
#' @description Estimates tree occurrence regions from a \strong{normalized} point cloud based on the Hough Transform circle search.
#' @template param-las
#' @template param-hmin-hmax
#' @template param-hstep
#' @template param-pixel-size
#' @template param-max-radius
#' @template param-min-density
#' @template param-min-votes
#' @template return-las
#' @section Output Header:
#'
#' Each point in the \code{LAS} object output represents a pixel center that is
#' \emph{possibly} also stem cross-section center.
#'
#' The variables describing each point in the output are:
#'
#' \itemize{
#' \item \code{Intensity}: number of votes received by that point
#' \item \code{PointSourceID}: stem segment ID (among all trees)
#' \item \code{Keypoint_flag}: if \code{TRUE}, the point is the most likely circle center
#' of its stem segment (\code{PointSourceID})
#' \item \code{Radii}: approximate radius estimated by that point - always a multiple of the \code{pixel_size}
#' \item \code{TreeID}: (possible) unique tree ID of the point
#' \item \code{TreePosition}: if \code{TRUE}, the point represents its tree's approximate coordinate
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
#' @useDynLib TreeLS, .registration = TRUE
#' @export
treeMap = function(las, hmin = 1, hmax = 3, hstep = 0.5, pixel_size = 0.025, max_radius = 0.25, min_density = 0.1, min_votes = 3){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  if(hmax <= hmin)
    stop('hmax must be larger than hmin')

  params = list(
    hstep = hstep,
    pixel_size = pixel_size,
    max_radius = max_radius,
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
    stop('min_den must be between 0 and 1')

  rgz = las$Z %>% range

  if(hmax < rgz[1])
    stop('hmax is too low - below the point cloud')

  if(hmin > rgz[2])
    stop('hmin is too high - above the point cloud')

  preCheck(las)

  if("Classification" %in% names(las@data))
    las %<>% lasfilter(Classification != 2)

  map = stackMap(las %>% las2xyz, hmin, hmax, hstep, pixel_size, max_radius, min_density, min_votes) %>%
    do.call(what=cbind) %>% as.data.frame

  map$Intensity %<>% as.integer
  map$Keypoint_flag %<>% as.logical
  map$PointSourceID %<>% as.integer
  map$TreePosition %<>% as.logical
  map %<>% LAS %>% setHeaderTLS

  attributes(map)[[tls.marker]] = 'tree_map'

  return(map)
}


#' Get tree XY positions from a TLS-tree map
#' @description Extracts the tree positions from a \code{LAS} object tree map
#' @param las \code{LAS} object - output from \code{\link{treeMap}}.
#' @param plot \code{logical} - plot the tree map?
#' @return \code{data.frame} of tree IDs and XY coordinates.
#' @template example-tree-map
#' @export
treePositions = function(las, plot=T){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  tlsatt = attributes(las)[[tls.marker]]
  if(is.null(tlsatt) || tlsatt != 'tree_map')
    stop('las is not a tree_map object: check ?treeMap')

  las %<>% lasfilter(TreePosition)

  pos = las@data[,c('TreeID', 'X', 'Y')] %>% as.data.frame

  attributes(pos)[[tls.marker]] = 'tree_map_df'

  if(plot){
    pos %$% plot(Y ~ X, cex=3, pch=20, main='tree map', xlab='X', ylab='Y')
  }

  return(pos)
}


#' Stem points classification
#' @description Classify stem points from a \strong{normalized} point cloud through the Hough Transform algorithm.
#' @template param-las
#' @param map \emph{optional} - map of tree positions (output from \code{\link{treeMap}} or \code{\link{treePositions}}). If omitted, the algorithm assumes that \code{las} is a single tree.
#' @template param-hstep
#' @template param-max-radius
#' @template param-hbase
#' @template param-pixel-size
#' @template param-min-density
#' @template param-min-votes
#' @template return-las
#' @section Output Header:
#'
#' Meaninful fields in the output:
#'
#' \itemize{
#' \item \code{TreeID}:  (possible) unique tree ID of the point - available when a tree \code{map} is provided
#' \item \code{Stem}: \code{TRUE} for stem points
#' \item \code{Segment}: stem segment number (from bottom to top)
#' \item \code{Radius}: approximate radius of the point's stem segment estimated by the Hough Transform - always a multiple of the \code{pixel_size}
#' \item \code{Votes}: votes received by the stem segment's center during the Hough Transform
#' }
#'
#' @template section-hough-transform
#' @template reference-thesis
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
#' thin = tlsSample(tls, 'voxel', 0.02)
#' map = treeMap(thin, min_density = 0.05)
#'
#' tls = stemPoints(tls, map)
#' plot(tls, color='Stem')
#'
#' @useDynLib TreeLS, .registration = TRUE
#' @export
stemPoints = function(las, map = NULL, hstep=0.5, max_radius=0.25, hbase = c(1,2.5), pixel_size=0.025, min_density=0.1, min_votes=3){

  if(class(las)[1] != 'LAS')
    stop('input las must be a LAS object')

  if(!is.null(map)){
    tlsatt = attributes(map)[[tls.marker]]

    if(is.null(tlsatt)){
      stop('las is not a tree_map object: check ?treeMap')
    }else if(tlsatt == 'tree_map_df'){
      map %<>% as.data.frame
    }else if(tlsatt == 'tree_map'){
      map %<>% treePositions(F)
    }else{
      stop('las is not a tree_map object: check ?treeMap')
    }

    if( map$TreeID %>% duplicated %>% any )
      stop('input map must have unique TreeIDs')

  }else{
    rg = apply(las@data[,1:2], 2, function(x) max(x) - min(x)) %>% as.double

    if(any(rg > 10))
      message("point cloud unlikely a single tree - XY extents too large")
  }

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

    if(length(val) != 1)
      stop( i %>% paste('must be of length 1') )

    if(!is.numeric(val))
      stop( i %>% paste('must be Numeric') )

    if(val <= 0)
      stop( i %>% paste('must be positive') )
  }

  if(min_density > 1)
    stop('min_density must be between 0 and 1')

  if(max(las$Z) < 0)
    stop('input Z coordinates are all negative')

  if(min(las$Z) < 0)
    message("points with Z below 0 will be ignored")

  if(min(las$Z) > 5)
    message("point cloud didn't look normalized - Z values too high")

  groundPts = if(las %>% hasField('Classification')){
    las$Classification == 2
  }else{
    rep(F, las@data %>% nrow)
  }

  if(map %>% is.null){
    message('no tree_map provided: performing single stem point classification')
    results = houghStemPoints(las2xyz(las)[!groundPts,], hbase[1], hbase[2], hstep, max_radius, pixel_size, min_density, min_votes)
  }else{
    message('performing point classification on multiple stems')
    results = houghStemPlot(las2xyz(las)[!groundPts,], map %>% as.matrix, hbase[1], hbase[2], hstep, max_radius, pixel_size, min_density, min_votes)
    las@data$TreeID = 0
    las@data$TreeID[!groundPts] = results$TreeID
  }

  las@data$Stem = F
  las@data$Stem[!groundPts] = results$Stem

  las@data$Segment = 0
  las@data$Segment[!groundPts] = results$Segment

  las@data$Radius = 0
  las@data$Radius[!groundPts] = results$Radius

  las@data$Votes = 0
  las@data$Votes[!groundPts] = results$Votes

  las %<>% resetLAS

  attributes(las)[[tls.marker]] = ifelse(map %>% is.null, "single_stem_points", "multiple_stem_points")

  return(las)

}


#' Stem segmentation
#' @description Measure stem segments from a classified point cloud through least squares circle fit and the RANSAC algorithm.
#' @template param-las
#' @param tol \code{numeric} - tolerance offset between absolute radii estimates and hough transform estimates.
#' @param n \code{integer} - number of points selected on every RANSAC iteration.
#' @param conf \code{numeric} - confidence level.
#' @param inliers \code{numeric} - expected proportion of inliers among stem segments' point cloud chunks.
#' @return \code{data.frame} of stem segmetns.
#' @template reference-thesis
#' @section Output Fields:
#'
#' \itemize{
#' \item \code{TreeID}:  unique tree IDs - available only for multiple stems
#' \item \code{Segment}: stem segment number (from bottom to top)
#' \item \code{X}, \code{Y}: circle center coordinate
#' \item \code{Radius}: estimated circle radius
#' \item \code{Error}: least squares circle fit error
#' \item \code{AvgHeight}: average height of stem segment
#' \item \code{N}: number of points in the segment
#' }
#'
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
#' thin = tlsSample(tls, 'voxel', 0.02)
#' map = treeMap(thin, min_density = 0.05)
#'
#' tls = stemPoints(tls, map)
#' df = stemSegmentation(tls)
#' head(df)
#'
#' @useDynLib TreeLS, .registration = TRUE
#' @export
stemSegmentation = function(las, tol=0.025, n = 10, conf = 0.99, inliers = 0.8){

  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object')

  params = list(
    tol = tol,
    n = n,
    conf = conf,
    inliers = inliers
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

  if(n < 3)
    stop('n must be at least 3')

  if(conf >= 1)
    stop('conf must be between 0 and 1')

  if(inliers >= 1)
    stop('inliers must be between 0 and 1')

  tlsatt = attributes(las)[[tls.marker]]
  if(is.null(tlsatt)){
    stop('stem points identifier missing - check ?stemPoints')
  }

  las %<>% lasfilter(Stem)

  if(tlsatt == 'single_stem_points'){

    message('performing single stem segmentation')

    estimates = ransacStem(las %>% las2xyz, las@data$Segment, las@data$Radius, n, conf, inliers, tol) %>% do.call(what = rbind) %>% as.data.frame
    names(estimates) = c('X', 'Y', 'Radius', 'Error', 'Segment')

    z = tapply(las@data$Z, las@data$Segment, mean)
    z = cbind(AvgHeight = z, Segment = z %>% names %>% as.double)

    np = table(las@data$Segment)
    np = cbind(N = np, Segment = np %>% names %>% as.double)

    estimates %<>% base::merge(z, by='Segment') %>% base::merge(np, by='Segment')
    attributes(estimates)[[tls.marker]] = "single_stem_df"

  }else if(tlsatt == 'multiple_stem_points'){

    message('performing multiple stems segmentation')

    estimates = ransacPlot(las %>% las2xyz, las$TreeID, las$Segment, las$Radius, n, conf, inliers, tol) %>% sapply(do.call, what=rbind) %>% do.call(what = rbind) %>% as.data.frame

    names(estimates) = c('X', 'Y', 'Radius', 'Error', 'Segment', 'TreeID')
    estimates = estimates[,c('TreeID', 'Segment', 'X', 'Y', 'Radius', 'Error')]

    z = apply(estimates, 1, function(x){
      temp = las@data$Z[las@data$Segment == x[2] & las@data$TreeID == x[1]]
      return(data.frame(AvgHeight = mean(temp), N = length(temp)))
    }) %>% do.call(what = rbind)

    estimates %<>% base::cbind(z)
    attributes(estimates)[[tls.marker]] = "multiple_stems_df"

  }else{
    stop('stem points identifier missing - check ?stemPoints')
  }

  estimates = estimates[ estimates$N > n , ]

  return(estimates)

}
