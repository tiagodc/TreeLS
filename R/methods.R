# ===============================================================================
#
# Developers:
#
# Tiago de Conto - ti@forlidar.com.br -  https://github.com/tiagodc/
#
# COPYRIGHT: Tiago de Conto, 2019
#
# This piece of software is open and free to use, redistribution and modifications
# should be done in accordance to the GNU General Public License >= 3#
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

. = X = Y = Z = Classification = TreePosition = TreeID = Stem = Segment = gpstime = AvgHeight = Radius = NULL

point.metrics.names = c('Planarity', 'Verticality', 'LinearSaliency', 'PlanarSaliency', 'Scattering', 'Anisotropy', 'Zrange', 'Zsd', 'N', 'EigenValue1', 'EigenValue2', 'EigenValue3', 'EigenVector11', 'EigenVector21', 'EigenVector31', 'EigenVector12', 'EigenVector22', 'EigenVector32', 'EigenVector13', 'EigenVector23', 'EigenVector33')

point.metrics.check = c('Planarity', 'Verticality', 'LinearSaliency', 'PlanarSaliency', 'Scattering', 'Anisotropy', 'Zrange', 'Zsd', 'N', 'EigenValues', 'EigenVectors', 'MeanDistance', 'MedianDistance', 'MinDistance', 'MaxDistance')

tls.marker = 'tlsAttribute'

isLAS = function(las){
  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object - checkout ?setTLS')
}

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

  isLAS(las)

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

  las = las@data[,c('X','Y','Z')] %>% as.matrix
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

#' @importFrom stats cov
planeAngle = function(xyz, axis='z'){

  e = eigen(cov(xyz))
  if(axis != 'z') e$vectors[3,3] = 0

  global_axis = if(axis == 'z') c(0,0,1) else if(axis=='x') c(1,0,0) else c(0,1,0)

  ang = (( e$vectors[,3] %*% global_axis ) / ( sqrt(sum(e$vectors[,3]^2)) * sqrt(sum(global_axis^2)) )) %>%
    as.double %>% acos

  return(ang)
}

rotationMatrix = function (ax, az, ax2){

  Rx = matrix(c(1, 0, 0, 0, cos(ax), sin(ax), 0, -sin(ax), cos(ax)), ncol = 3, byrow = T)
  Rz = matrix(c(cos(az), 0, -sin(az), 0, 1, 0, sin(az), 0, cos(az)), ncol = 3, byrow = T)
  Rx2 = matrix(c(cos(ax2), sin(ax2), 0, -sin(ax2), cos(ax2), 0, 0, 0, 1), ncol = 3, byrow = T)

  mat = Rx2 %*% Rz %*% Rx

  return(mat)
}

tfMatrix = function (ax, az, ax2, x, y, z){

  mat = rotationMatrix(ax, az, ax2) %>%
    rbind(0) %>% cbind(c(x,y,z,1))

  return(mat)
}

splitByIndex = function(las, var='Z', max_size = 1E6){

  npts = nrow(las@data)
  zclass = 0

  if(npts > max_size){
    nparts = ceiling(npts/max_size)
    probs = seq(0, 1, by = 1/nparts)
    probs = probs[probs > 0 & probs < 1]
    zqts = quantile(las[[var]], probs) %>% as.double
    hs = c(min(las[[var]])-1, zqts, max(las$Z)+1)
    zclass = cut(las[[var]], hs) %>% as.integer
  }

  return(zclass)
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
#' file = system.file("extdata", "pine.laz", package="TreeLS")
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
#' ## sample points systematically in 3D
#' vx = tlsSample(tls, voxelize(0.05))
#' summary(vx)
#'
#' ## sample points randomly
#' rd = tlsSample(tls, randomize(0.5))
#' summary(rd)
#'
#' @export
tlsSample = function(las, method = voxelize()){

  isLAS(las)

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

  isLAS(las)

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
#' file = system.file("extdata", "pine_plot.laz", package="TreeLS")
#' tls = readTLS(file)
#' plot(tls)
#' rgl::axes3d(col='white')
#'
#' ### remove topography effect
#' tls = tlsNormalize(tls, 0.5, FALSE)
#' plot(tls)
#' rgl::axes3d(col='white')
#'
#' @importFrom raster raster extent res<-
#' @export
tlsNormalize = function(las, res=.5, keepGround=TRUE){

  isLAS(las)

  if(res <= 0)
    stop('res must be a positive number')

  if(!any(las$Classification == 2)){
    message('no ground points found, performing ground segmentation')
    las %<>% lasground(csf(class_threshold = 0.2, cloth_resolution = 0.1), last_returns = F)
  }

  grid = las %>% extent %>% raster
  res(grid) = res

  dtm = grid_terrain(las, res = grid, algorithm = knnidw())

  las %<>% lasnormalize(dtm, na.rm=TRUE)

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

  isLAS(las)

  if(!hasAttribute(method, 'tls_map_mtd'))
    stop('invalid method: check ?treeMap')

  preCheck(las)

  if(hasField(las, 'Classification'))
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

  isLAS(las)

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
#' @param method stem denoising algorithm - currently available: \code{\link{stem.hough}}.
#' @return \code{LAS} object with \emph{stem_points} signature.
#' @examples
#' ### single tree
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#' tls = stemPoints(tls)
#' plot(tls, color='Stem')
#'
#' ### forest plot - check the example given for the stem segmentation function
#' ?stemSegmentation
#'
#' @export
stemPoints = function(las, method = stem.hough()){

  isLAS(las)

  if(!hasField(las, 'TreeID')){
    rg = apply(las@data[,1:2], 2, function(x) max(x) - min(x)) %>% as.double

    if(any(rg > 10))
      message("point cloud unlikely a single tree (XY extents too large) - check ?tlsCrop or ?treePoints")
  }

  if(max(las$Z) < 0)
    stop('input Z coordinates are all negative')

  if(!hasAttribute(method, 'stem_pts_mtd'))
    stop('invalid method: check ?stemPoints')

  if(abs(min(las$Z)) > 0.5)
    warning('point cloud apparently not normalized')

  las %<>% method()

  return(las)

}


#' Stem segmentation
#' @description Measure stem segments from a classified point cloud.
#' @template param-las
#' @param method stem segmentation algorithm - currently available: \code{\link{sgmt.ransac.circle}}.
#' @return \code{data.table} of stem segments with \emph{stem_dt} signature.
#' @template example-segmentation
#' @export
stemSegmentation = function(las, method=sgmt.ransac.circle()){

  isLAS(las)

  if(!hasAttribute(method, 'stem_sgmt_mtd'))
    stop('invalid method: check ?stemSegmentation')

  las %>% method %>% return()

}


#' Filter points based on gpstime
#' @description This is a simple wrapper to \code{\link[lidR:lasfilter]{lasfilter}} that takes as inputs proportional values instead of absolute time stamp values for filtering a point cloud object based on the gpstime field. This function is particularly useful to check narrow intervals of point cloud frames from mobile scanning data.
#' @template param-las
#' @param from,to \code{numeric} - between 0 and 1 - gpstime quantile limits to filter by
#' @template return-las
#' @importFrom stats quantile
#' @examples
#' file = system.file("extdata", "model_boles.laz", package="TreeLS")
#' tls = readTLS(file)
#'
#' ### color points according to its chronological time stamp
#' plot(tls, color='gpstime')
#'
#' ### keep points registered in the 70% to 95% time interval
#' tls = gpsTimeFilter(tls, .7, .95)
#' plot(tls, color='gpstime')
#'
#' @export
gpsTimeFilter = function(las, from=0, to=1){

  isLAS(las)

  if(!hasField(las, 'gpstime'))
    stop('LAS object has no gpstime field')

  if(length(from) != 1 || from < 0 || from > 1)
    stop('from must be a single value between 0 and 1')

  if(length(to) != 1 || to < 0 || to > 1)
    stop('to must be a single value between 0 and 1')

  if(to <= from)
    stop('to must be larger than from')

  if(all(las@data$gpstime == las@data$gpstime[1])){
    message('all observations have the same gpstime stamp - no filtering performed')
    return(las)
  }

  qts = las@data$gpstime %>% quantile(c(from, to), na.rm=T)
  las %<>% lasfilter(gpstime >= qts[1] & gpstime <= qts[2])

  return(las)

}


#' Rotate point cloud towards a horizontal plane
#' @description Check for ground points and rotates the point cloud to align its ground surface to a horizontal plane (XY).
#' This function is especially useful for point clouds not georreferenced or generated through mobile scanning,
#' which might present a tilted global reference system. Since the coordinates are altered in this procedure, any
#' geographical information is erased from the LAS' header after rotation.
#' @template param-las
#' @template return-las
#' @examples
#' file = system.file("extdata", "pine_plot.laz", package="TreeLS")
#' tls = readTLS(file)
#'
#' ### note the tilted ground
#' plot(tls)
#' rgl::axes3d(col='white')
#'
#' ### after rotation
#' tls = tlsRotate(tls)
#' plot(tls)
#' rgl::axes3d(col='white')
#'
#' @export
tlsRotate = function(las){

  isLAS(las)

  . = NULL

  ground = las@data[,1:3] %>%
    toLAS %>%
    lasground(csf(class_threshold = .2), F) %>%
    lasfilter(Classification == 2)

  center = ground@data[,.(mean(X), mean(Y))] %>% as.double
  ground = tlsCrop(ground, center[1], center[2], 5) %>% las2xyz

  az = planeAngle(ground, 'z')
  ax = planeAngle(ground, 'x')
  ay = planeAngle(ground, 'y')

  rz = ifelse(az > pi/2, pi-az, -az)
  rx = ifelse(ay < pi/2, -ax, ax)

  rot = rotationMatrix(0, rz, rx) %>% as.matrix
  xyBack = rotationMatrix(0,0,-rx) %>% as.matrix

  minXYZ = apply(las@data[,1:3], 2, min) %>% as.double

  las@data$X = las@data$X - minXYZ[1]
  las@data$Y = las@data$Y - minXYZ[2]
  las@data$Z = las@data$Z - minXYZ[3]

  las@data[,1:3] = (las2xyz(las) %*% rot) %*% xyBack %>% as.data.table

  las@data$X = las@data$X + minXYZ[1]
  las@data$Y = las@data$Y + minXYZ[2]
  las@data$Z = las@data$Z + minXYZ[3]

  las %<>% resetLAS

  return(las)
}


#' Alter point cloud's coordinates
#' @description Apply transformations to the XYZ axes of a point cloud.
#' @template param-las
#' @param xyz \code{character} vector of length 3 - \code{las}' columns to be reassigned as XYZ, respectively.
#' Use minus signs to mirror the axes' coordinates - more details in the sections below.
#' @param bring_to_origin \code{logical} - force output to start at \code{c(0,0,0)}? If \code{TRUE},
#' removes any geographical information from the output.
#' @param rotate  \code{logical} - rotate the point cloud to align the ground points horizontally? If \code{TRUE},
#' removes any geographical information from the output. Checkout \code{\link{tlsRotate}} for more information.
#' @template return-las
#' @section XYZ Manipulation:
#'
#' The \code{xyz} argument can take a few different forms, it is useful to shift axes positions in a point cloud or
#' to mirror an axis' coordinates. All axes characters can be entered in lower or uppercase and also be preceeded
#' by a minus sign ('-'), indicating to invert (mirror) the axis' coordinates in the output.
#' Check the \emph{examples} section for a practical overview.
#'
#' @examples
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#'
#' ### swap the Y and Z axes
#' zy = tlsAlter(tls, c('x', 'z', 'y'))
#' bbox(zy)
#' range(zy$Z)
#' plot(zy, clear_artifacts=FALSE)
#' rgl::axes3d(col='white')
#'
#' ### return an upside down point cloud
#' ud = tlsAlter(tls, c('x', 'y', '-z'))
#' bbox(ud)
#' range(ud$Z)
#' plot(zy, clear_artifacts=FALSE)
#' rgl::axes3d(col='white')
#'
#' ### mirror all axes, then set the point cloud's starting point as the origin
#' rv = tlsAlter(tls, c('-x', '-y', '-z'), bring_to_origin=TRUE)
#' bbox(rv)
#' range(rv$Z)
#' plot(rv, clear_artifacts=FALSE)
#' rgl::axes3d(col='white')
#'
#' @export
tlsAlter = function(las, xyz = c('X', 'Y', 'Z'), bring_to_origin = FALSE, rotate = FALSE){

  isLAS(las)

  xyz %<>% toupper

  if(typeof(xyz) != 'character' || length(xyz) != 3)
    stop('xyz must be of type character and have length of exactly 3')

  if(!is.logical(bring_to_origin))
    stop('bring_to_origin must be logical')

  if(!is.logical(rotate))
    stop('rotate must be logical')

  XYZ = lapply(xyz, function(i){
    temp =
      if(i == 'X'){
        las$X
      }else if(i == '-X'){
        -las$X
      }else if(i == 'Y'){
        las$Y
      }else if(i == '-Y'){
        -las$Y
      }else if(i == 'Z'){
        las$Z
      }else if(i == '-Z'){
        -las$Z
      }else{
        NULL
      }

    if(temp %>% is.null)
      stop('invalid input in xyz:' %>% paste(i))

    return(temp)
  }) %>% do.call(what = cbind) %>% as.data.table

  . = NULL

  las@data[,1:3] = XYZ
  las %<>% resetLAS

  if(rotate){
    las %<>% tlsRotate()
  }

  if(bring_to_origin){
    las = LAS(las@data)

    mincoords = las@data[,.(min(X), min(Y), min(Z))] %>% as.double

    las@data$X = las@data$X - mincoords[1]
    las@data$Y = las@data$Y - mincoords[2]
    las@data$Z = las@data$Z - mincoords[3]

    las %<>% resetLAS
  }

  return(las)

}


#' Plot TLS outputs
#' @description Plot the \code{LAS} outputs of tls functions on the same scene using \code{rgl}. Check ?stemSegmentation
#' for usage examples.
#' @param las \code{LAS} object - ideally an output from \code{\link{stemPoints}}.
#' @param sgmt optional \code{data.table} - output from \code{\link{stemSegmentation}}.
#' @param map optional \code{LAS} object - output from \code{\link{treeMap}}.
#' @param treeID optional \code{numeric} - single \emph{TreeID} to extract from \code{las}.
#' @param sgmtColor optional - color of the plotted stem segment representations.
#' @examples
#' ### single tree
#' file = system.file("extdata", "spruce.laz", package="TreeLS")
#' tls = readTLS(file)
#' tls = stemPoints(tls)
#' df = stemSegmentation(tls)
#'
#' tlsPlot(tls, df)
#'
#' ### For further examples check:
#' ?stemSegmentation
#' @export
tlsPlot = function(las, sgmt = NULL, map = NULL, treeID = NULL, sgmtColor = 'yellow'){

  isLAS(las)

  if(!is.null(treeID) && (length(treeID) != 1 || !is.numeric(treeID)))
      stop('treeID must be numeric of length 1')

  . = NULL

  cp = pastel.colors(100)

  if(hasAttribute(las, 'multiple_stem_points')){

    if(is.null(treeID)){

      temp = lasfilter(las, Stem)@data[,.(X = mean(X), Y = mean(Y), Z = min(Z)-0.5), by=TreeID]

      corner = plot(las %>% lasfilter(Stem), color='TreeID', colorPalette=cp, size=1.5)
      las %<>% lasfilter(!Stem & Classification != 2)
      rgl.points(las$X - corner[1], las$Y - corner[2], las$Z, color='white', size=.5)

      temp %$% text3d(X - corner[1], Y - corner[2], Z, TreeID, cex=1.5, col=sgmtColor)

    }else{

      xy = las@data[TreeID == treeID, .(mean(X),mean(Y))]

      if(nrow(xy) == 0)
        stop('no match for treeID ==' %>% paste(treeID))

      xy %<>% as.double

      las %<>% tlsCrop(xy[1], xy[2], 1.5)

      temp = lasfilter(las, Stem)@data[,.(X = mean(X), Y = mean(Y), Z = min(Z)-0.5), by=TreeID]

      corner = plot(las %>% lasfilter(Stem), size=1.5)
      las %<>% lasfilter(!Stem & Classification != 2)
      rgl.points(las$X - corner[1], las$Y - corner[2], las$Z, color='white', size=.5)

      temp %$% text3d(X - corner[1], Y - corner[2], Z, TreeID, cex=1.5, col=sgmtColor)

    }

  }else if(hasAttribute(las, 'single_stem_points')){

    corner = plot(las %>% lasfilter(Stem), size=1.5)
    las %<>% lasfilter(!Stem & Classification != 2)
    rgl.points(las$X - corner[1], las$Y - corner[2], las$Z, color='white', size=.5)

  }else{
    corner = plot(las, size=1.5)
  }

  if(!is.null(sgmt)){
    if(hasAttribute(sgmt, 'multiple_stems_dt') || hasAttribute(sgmt, 'single_stem_dt')){

      if(!is.null(treeID) && hasAttribute(sgmt, 'multiple_stems_dt')){
        sgmt = sgmt[TreeID == treeID]
      }

      sgmt %$% spheres3d(X-corner[1], Y-corner[2], AvgHeight, Radius, color=sgmtColor)
    }else{
      message('sgmt does not have a stem_dt signature')
    }
  }

  if(!is.null(map)){
    if(hasAttribute(map, 'tree_map')){
      if(!is.null(treeID)){
        map %<>% lasfilter(TreeID == treeID)
      }
      map@data %$% rgl.points(X - corner[1], Y - corner[2], Z, color='green', size=.5)
    }else{
      message('map does not have a tree_map signature')
    }
  }

  return(corner %>% invisible)

}


####### NEW IMPLEMENTATIONS

treePoints = function(las, map, method=trees.voronoi()){

  isLAS(las)

  if(hasAttribute(map, 'tree_map_dt')){
    # map = map
  }else if(hasAttribute(map, 'tree_map')){
    map %<>% treePositions(F)
  }else{
    stop('las is not a tree_map object: check ?treeMap')
  }

  if( map$TreeID %>% duplicated %>% any )
    stop('input map must have unique TreeIDs')

  las %<>% method(map)
  return(las)

}

nnFilter = function(las, d = 0.05, n = 2, max_points = 1E6){

  zclass = splitByIndex(las, max_size = max_points)
  keep = rep(T, nrow(las@data))

  for(i in unique(zclass)){
    bool = zclass == i
    xyz = las@data[bool,.(X,Y,Z)]
    rnn = RANN::nn2(data = xyz, k = n+1, treetype = 'kd')$nn.dists[,-1] %>% as.data.frame

    knn = rep(0, nrow(rnn))
    for(j in 1:ncol(rnn)){
      temp = ifelse(rnn[,j] > d, 1, 0) %>% as.double
      knn = knn + temp
    }

    keep[bool] = knn < n
  }

  las %<>% lasfilter(keep)
  return(las)
}

ptm.voxels = function(las, d = .1, exact=F, metrics_list = point.metrics.check){

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

  vtm$VoxelID = names(idx) %>% as.double
  las@data$VoxelID = vx
  las@data = merge(las@data, vtm, by='VoxelID', sort=F)
  las@data = las@data[,-c('VoxelID')]
  las@data$VoxelID = vx

  return(las)

}

ptm.knn = function(las, k = 30, metrics_list = point.metrics.check){
  zclass = splitByIndex(las)
  zuq = unique(zclass)

  df = data.table()
  for(i in zuq){
    temp = lasfilter(las, zclass == i)
    knn = RANN::nn2(temp %>% las2xyz, k = k, treetype = 'kd', searchtype = 'standard')
    ptm = ptmStatistics(las, knn, metrics_list)
    temp@data[,colnames(ptm)] = ptm
    df = rbind(df, temp@data)
  }

  las@data = df
  las = resetLAS(las)
  return(las)
}

ptm.radius = function(las, r = 0.1, max_k = 30, metrics_list = point.metrics.check){
  zclass = splitByIndex(las)
  zuq = unique(zclass)

  df = data.table()
  for(i in zuq){
    temp = lasfilter(las, zclass == i)
    knn = RANN::nn2(las %>% las2xyz, k = max_k, treetype = 'kd', searchtype = 'radius', radius = r)
    ptm = ptmStatistics(las, knn, metrics_list)
    temp@data[,colnames(ptm)] = ptm
    df = rbind(df, temp@data)
  }

  las@data = df
  las = resetLAS(las)
  return(las)
}

ptmStatistics = function(las, knn, metrics_list = point.metrics.check){

  pickMetrics = ptmMetricsLog(metrics_list)

  kid = knn$nn.idx
  kds = knn$nn.dists
  kds[kid == 0] = NA

  ptm = data.table()

  if(any(pickMetrics$log[1:11])){
    ptm =  pointMetrics(las %>% las2xyz, kid, pickMetrics$log) %>% do.call(what = rbind) %>% as.data.table
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

ptmMetricsLog = function(metrics_list){
  metrics_log = point.metrics.check %in% metrics_list

  if(all(!metrics_log)) stop('Please provide at least one known metric. See ?availablePointMetrics')

  metrics_names = point.metrics.names[1:9][metrics_log[1:9]]
  if(metrics_log[10]) metrics_names %<>% c(point.metrics.names[10:12])
  if(metrics_log[11]) metrics_names %<>% c(point.metrics.names[13:21])

  return(list(log = metrics_log, names = metrics_names))

}

availablePointMetrics = function(){
  point.metrics.check %>% as.data.frame %>% print
  cat('\n')
  return(point.metrics.check)
}
