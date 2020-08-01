# ===============================================================================
#
# Developers:
#
# Tiago de Conto - tdc.florestal@gmail.com -  https://github.com/tiagodc/
#
# COPYRIGHT: Tiago de Conto, 2019
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

#' @import data.table
#' @import lidR
#' @import magrittr
#' @import rgl
#' @import dismo
#' @import deldir
#' @import nabor
#' @useDynLib TreeLS, .registration = TRUE

. = X = Y = Z = Classification = TreePosition = TreeID = Stem = Segment = gpstime = AvgHeight = Radius = NULL

TLS_MARKER = 'TLS_MARKER'

isLAS = function(las){
  if(class(las)[1] != 'LAS')
    stop('input data must be a LAS object - checkout ?setTLS')
}

setAttribute = function(obj, attribute_name){
  attr(obj, TLS_MARKER) = attribute_name
  return(obj)
}

hasAttribute = function(obj, attribute_name){
  tlsatt = attr(obj, TLS_MARKER)
  bool = is.null(tlsatt) || tlsatt != attribute_name
  return(!bool)
}

plot.cylinder = function(x_center = 0, y_center = 0, h_bottom = 0, h_top = 1, radius = 0.5, color = 'yellow'){

  axis = matrix(c(
    rep(x_center, 2),
    rep(y_center, 2),
    seq(h_bottom, h_top, length.out = 2)
  ), ncol = 3, byrow = F)

  cyl = cylinder3d(axis, radius = radius)

  mesh = shade3d(addNormals(subdivision3d(cyl, depth = 0)), col = color)
  # mesh = shade3d(cyl, col=col)
}

preCheck = function(las){

  isLAS(las)

  has_class = "Classification" %in% names(las@data)

  if(!has_class){
    message('no Classification field found in the dataset')

  }else{

    has_ground = any(las$Classification == 2)

    if(has_ground){

      mean_ground = las$Z[ las$Classification == 2 ] %>% mean(na.rm=T) %>% abs

      if(mean_ground > 0.2)
        message("point cloud apparently not normalized")
    } else {
      message("ground points not classified")
    }
  }

}

#' @importFrom benchmarkme get_ram
sizeCheck = function(las, n_new_fields, bytes=8){
  n = nrow(las@data)
  mem = sum(gc(full=TRUE)[,2])
  ram = as.double(get_ram()) / 1000000
  malloc = n * n_new_fields * bytes / 1000000
  can_malloc = ram - mem - malloc
  if(can_malloc < 0){
    paste('adding', n_new_fields, 'columns to the las object is not possible - not enough RAM available') %>% stop
  }else if((can_malloc / ram) < 0.1){
    paste('adding', n_new_fields, 'columns to the las object will take more than 90% of the available RAM') %>% warning
  }
}

toLAS = function(data_matrix, column_names=NULL){

  if(ncol(data_matrix) < 3)
    stop('data_matrix must have at least 3 columns')

  data_matrix %<>% as.data.table

  if(!is.null(column_names)){

    if(length(column_names) != ncol(data_matrix))
      stop('data_matrix must have the same number of columns as in column_names')

    checkXYZ = c('X', 'Y', 'Z') %in% column_names

    if(!all(checkXYZ))
      stop('X, Y and Z must be declared explicitly in column_names (in uppercase)')

    colnames(data_matrix) = column_names

  } else {

    if(ncol(data_matrix) > 3)
      message('converting first three columns only (assumed XYZ coordinates)')

    data_matrix = data_matrix[,1:3]
    colnames(data_matrix) = c('X', 'Y', 'Z')

  }

  data_matrix = suppressMessages(data_matrix %>% LAS %>% setHeaderTLS)
  return(data_matrix)
}

las2xyz = function(las){

  if(class(las)[1] != "LAS")
    stop("las must be a LAS object")

  las = las@data[,c('X','Y','Z')] %>% as.matrix
  return(las)
}

hasField = function(las, field_name){
  if(class(las)[1] == 'LAS') las = las@data
  any(colnames(las) == field_name) %>% return()
}

cleanFields = function(las, field_names){
  is_las = class(las)[1] == 'LAS'
  for(i in field_names){
    temp = if(is_las) las@data[,..i] else las[,..i]
    temp = unlist(temp)
    temp[is.na(temp) | is.nan(temp) | is.infinite(temp) | is.null(temp)] = ifelse(is.logical(temp), F, 0)
    if(is_las) las@data[,i] = temp else las[,i] = temp
  }
  return(las)
}

setHeaderTLS = function(las, x_scale = 0.0001, y_scale = 0.0001, z_scale = 0.0001){

  if(class(las)[1] != "LAS")
    stop("las must be a LAS object")

  if(las@header@PHB$`X scale factor` < x_scale)
    x_scale = las@header@PHB$`X scale factor`

  if(las@header@PHB$`Y scale factor` < y_scale)
    y_scale = las@header@PHB$`Y scale factor`

  if(las@header@PHB$`Z scale factor` < z_scale)
    z_scale = las@header@PHB$`Z scale factor`

  return(suppressMessages(las_rescale(las, x_scale, y_scale, z_scale)))
}

#' @importFrom stats runif
tlsCylinder = function(n=10000, h=100, radius=30, deviation=0){

  radius = runif(n, radius-deviation, radius+deviation)

  z=runif(n = n, min = 0, max = h)

  angs = runif(n, 0, 2*pi)
  x = sin(angs)*radius
  y = cos(angs)*radius

  return(cbind(x,y,z) %>% toLAS)
}

#' @importFrom stats cov
planeAngle = function(xyz, axis='z'){

  e = eigen(cov(xyz))
  if(axis != 'z') e$vectors[3,3] = 0

  global_axis = if(axis == 'z') c(0,0,1) else if(axis=='x') c(1,0,0) else c(0,1,0)

  angle = (( e$vectors[,3] %*% global_axis ) / ( sqrt(sum(e$vectors[,3]^2)) * sqrt(sum(global_axis^2)) )) %>%
    as.double %>% acos

  return(angle)
}

rotationMatrix = function (ax, ay, az){

  Rx = matrix(c(1, 0, 0, 0, cos(ax), sin(ax), 0, -sin(ax), cos(ax)), ncol = 3, byrow = T)
  Ry = matrix(c(cos(ay), 0, -sin(ay), 0, 1, 0, sin(ay), 0, cos(ay)), ncol = 3, byrow = T)
  Rz = matrix(c(cos(az), sin(az), 0, -sin(az), cos(az), 0, 0, 0, 1), ncol = 3, byrow = T)

  mat = Rz %*% Ry %*% Rx

  return(mat)
}

tfMatrix = function(ax, ay, az, x, y, z){

  mat = rotationMatrix(ax, ay, az) %>%
    rbind(0) %>% cbind(c(x,y,z,1))

  return(mat)
}

rangeMeans = function(X,Y,Z){
  meds = apply(cbind(X,Y,Z), 2, function(x) sum(range(x))/2) %>% t %>% as.data.table
  names(meds) = c('PX','PY','PZ')
  return(meds)
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
#' cloud = matrix(runif(300, 0, 10), ncol=3)
#' cloud = setTLS(cld)
#' summary(cld)
#' @export
setTLS = function(cloud, col_names=NULL){

  if(class(cloud)[1] == 'LAS'){
    cloud = setHeaderTLS(cloud)
  }else{
    cloud = toLAS(cloud, col_names)
  }

  return(cloud)
}


#' Import a point cloud file into a LAS object
#' @description Wrapper to read point clouds straight to LAS objects suitable for TLS applications. Reads \emph{las} or \emph{laz} files with \code{\link[lidR:readLAS]{readLAS}} and \emph{ply} files with \code{\link[rlas:read.las]{read.las}}. Other file formats are (or try to be) read using \code{\link[data.table:fread]{fread}}.
#' @param file file path.
#' @template param-colnames
#' @param ... further arguments passed to \code{readLAS}, \code{read.las} or \code{fread}.
#' @template return-las
#' @examples
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#' summary(tls)
#' @importFrom utils read.table
#' @importFrom rlas read.las
#' @importFrom data.table fread
#' @export
readTLS = function(file, col_names=NULL, ...){

  format = sub('.+\\.(.+$)', '\\1', file) %>% tolower

  if(format %in% c('laz', 'las')){

    las = readLAS(file, ...) %>% setHeaderTLS

  }else if(format == 'ply'){

    las = LAS(read.las(file, ...)) %>% setHeaderTLS

  }else{

    las = fread(file, ...) %>% toLAS(col_names)

  }

  return(las)
}


#' Nearest neighborhood point filter
#' @description Remove points from a \code{LAS} point cloud based on nearest neighborhood metrics.
#' @template param-las
#' @param d \code{numeric} - search radius.
#' @param n \code{integer} - number of nearest points to search within \code{d} radius.
#' @template return-las
#' @examples
#' file = system.file("extdata", "spruce.laz", package="TreeLS")
#' tls = readTLS(file)
#' summary(tls)
#'
#' nn_tls = nnFilter(tls, 0.05, 3)
#' summary(nn_tls)
#' @importFrom nabor knn
#' @export
nnFilter = function(las, d = 0.05, n = 2){

  isLAS(las)

  if(d <= 0)
    stop('d must be a number larger than 0')

  if(n <= 0 )
    stop('d must be an integer larger than 0')

  sizeCheck(las, n)

  rnn = knn(las %>% las2xyz, k = n+1)$nn.dists[,-1]

  keep = rep(T, nrow(las@data))
  for(i in 1:ncol(rnn)){
    keep = keep & rnn[,i] < d
  }

  las = filter_poi(las, keep)
  return(las)
}


#' Calculate metrics on point neighborhoods
#' @description Get statistics for every point's neighborhood. Check out \code{\link{fastPointMetrics.available}} for information about available metrics. Neighborhood search methods are prefixed by \code{ptm}.
#' @template param-las
#' @param method neighborhood search algorithm - currently available: \code{\link{ptm.voxels}} and \code{\link{ptm.knn}}.
#' @param which_metrics optional \code{character} vector - list of metrics to be calculated. Overwrites the gobally enabled metrics set with \code{fastPointMetrics.available()}.
#' @return \code{LAS} object with extra fields - one for each calculated metric.
#' @examples
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#' @export
fastPointMetrics = function(las, method = ptm.voxels(), which_metrics = ENABLED_POINT_METRICS){

  isLAS(las)

  if(!hasAttribute(method, 'ptm_mtd'))
    stop('invalid method: check ?fastPointMetrics')

  return(method(las, which_metrics))
}


#' Print available point metrics
#' @description print available point metrics - used in \code{\link{fastPointMetrics}}.
#' @param enable optional \code{integer} or \code{character} vector. Indices or names of the metrics you want to enable globally. Enabled metrics are calculated everytime you run \code{\link{fastPointMetrics}}.
#' @return \code{character} vector of available metrics.
#' @examples
#' m = fastPointMetrics.available()
#' print(m)
#' @export
fastPointMetrics.available = function(enable = ENABLED_POINT_METRICS){
  if(typeof(enable) != 'character') enable = POINT_METRICS_NAMES[enable]
  enable_metrics = ptmMetricsLog(enable)
  ENABLED_POINT_METRICS <<- enable_metrics$names
  temp = data.frame(INDEX = 1:length(POINT_METRICS_NAMES), METRIC = POINT_METRICS_NAMES, ENABLED = enable_metrics$log)
  print(temp)
  cat('\n')
  return(invisible(POINT_METRICS_NAMES))
}


#' Resample a point cloud
#' @description Applies an algorithm that reduces the point cloud's density. Sampling methods are prefixed by \code{smp}.
#' @template param-las
#' @param method point sampling algorithm - currently available: \code{\link{smp.voxelize}} or \code{\link{smp.randomize}}
#' @template return-las
#' @examples
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#' summary(tls)
#'
#' ## sample points systematically from a 3D voxel grid
#' vx = tlsSample(tls, smp.voxelize(0.05))
#' summary(vx)
#'
#' ## sample half of the points randomly
#' rd = tlsSample(tls, smp.randomize(0.5))
#' summary(rd)
#'
#' @export
tlsSample = function(las, method = smp.voxelize()){

  isLAS(las)

  if(!hasAttribute(method, 'tls_sample_mtd'))
    stop('invalid method: check ?tlsSample')

  las %<>% filter_poi(method(las))

  return(las)
}


#' Point cloud cropping
#' @description Returns a cropped point cloud of all points inside or outside specified boundaries of circle or square shapes.
#' @template param-las
#' @param x,y \code{numeric} -  X and Y center coordinates of the crop region.
#' @param len \code{numeric} -  if \code{circle = TRUE}, \code{len} is the circle's radius, otherwise it is the side length of a square.
#' @param circle \code{logical} -  if \code{TRUE} (default), crops a circle, otherwise a square.
#' @param negative \code{logical} - if \code{TRUE}, returns all points outside the specified circle/square perimeter, otherwise returns all points inside the circle/square (default).
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
  las %<>% filter_poi(bool)

  return(las)

}


#' Normalize a TLS point cloud
#' @description Fast normalization of TLS point clouds based on a Digital Terrain Model (DTM) of the ground points. If the input's ground points are not yet classified, the \code{\link[lidR:csf]{csf}} algorithm is applied internally.
#' @template param-las
#' @param min_res \code{numeric} - minimum resolution of the DTM used for normalization.
#' @param keep_ground \code{logical} - if \code{FALSE} removes the ground points from the output.
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
tlsNormalize = function(las, min_res=.25, keep_ground=TRUE){

  isLAS(las)

  if(min_res <= 0)
    stop('res must be a positive number')

  if(!any(las$Classification == 2)){
    message('no ground points found, performing ground segmentation')
    las = classify_ground(las, csf(class_threshold = 0.05, cloth_resolution = 0.05), last_returns = F)
  }

  res = area(las) / nrow(las@data[Classification == 2])
  res = ifelse(res < min_res, min_res, res)

  grid = las %>% extent %>% raster
  res(grid) = res

  dtm = grid_terrain(las, res = grid, algorithm = knnidw(), full_raster=TRUE)
  las = normalize_height(las, dtm, na.rm=TRUE, Wdegenerated = TRUE)

  if(!keep_ground) las = filter_poi(las, Classification != 2)

  return(las)

}


#' Map tree occurrences from TLS data
#' @description Estimates tree occurrence regions from a \strong{normalized} point cloud. Tree mapping methods are prefixed by \code{map}.
#' @template param-las
#' @param method tree mapping algorithm - currently available: \code{\link{map.hough}}, \code{\link{map.eigen.knn}}, \code{\link{map.eigen.voxel}} and \code{\link{map.pick}}.
#' @param merge \code{numeric} - parameter passed to \code{\link{treeMap.merge}}.
#' @param positions_only \code{logical} - if \code{TRUE} returns a 3D \code{LAS} object, otherwise returns a 2D tree map as a \code{data.table}.
#' @return signed \code{LAS} or \code{data.table}.
#' @template example-tree-map
#' @export
treeMap = function(las, method = map.hough(), merge=0.2, positions_only=FALSE){

  isLAS(las)

  if(!hasAttribute(method, 'tls_map_mtd'))
    stop('invalid method: check ?treeMap')

  preCheck(las)

  if(hasField(las, 'Classification'))
    las %<>% filter_poi(Classification != 2)

  map = method(las) %>% setAttribute('tree_map')

  if(merge > 0){
    map = treeMap.merge(map, merge)
  }

  if(positions_only){
    map = treeMap.positions(map, FALSE)
  }

  return(map)
}


#' Convert a tree map to a 2D \code{data.table}
#' @description Extracts the tree XY positions from a \emph{tree_map} \code{LAS} object.
#' @param las \code{LAS} object - output from \code{\link{treeMap}}.
#' @param plot \code{logical} - plot the tree map?
#' @return \code{data.table} of tree IDs and XY coordinates with \emph{tree_map_dt} signature.
#' @template example-tree-map
#' @export
treeMap.positions = function(las, plot=TRUE){

  if(!hasAttribute(las, 'tree_map') && !hasAttribute(las, 'tree_map_dt'))
    stop('las is not a tree_map object: check ?treeMap')

  if(hasAttribute(las, 'tree_map_dt')){
    pos = las
  }else{
    if(hasField(las, 'TreePosition')){
      las %<>% filter_poi(TreePosition)
    }

    pos = las@data[,.(X=median(X), Y=median(Y)),by=TreeID]
    pos = pos[order(TreeID)]

    pos %<>% setAttribute('tree_map_dt')
  }

  if(plot){
    pos %$% plot(Y ~ X, cex=3, pch=20, main='tree map', xlab='X', ylab='Y')
  }

  return(pos)
}


#' Merge nearby trees in a \emph{tree_map}
#' @description Check all tree neighborhoods and merge TreeIDs which are too close in a \code{tree_map} object.
#' @param las \code{LAS} or \code{data.table} object generated by \code{\link{treeMap}}.
#' @param d \code{numeric} - distance criteria to merge
#' @importFrom nabor knn
#' @export
treeMap.merge = function(las, d=.2){

  if( !(hasAttribute(las, 'tree_map') || hasAttribute(las, 'tree_map_dt')) )
    stop('las is not a tree_map object: check ?treeMap')

  if(d < 0)
    stop('d must be a positive number')

  is_data_table = hasAttribute('tree_map_dt')
  nxy = if(is_data_table) las else treeMap.positions(las, plot = F)
  nn = knn(nxy[,-1], k=2)
  dst = nn$nn.dists[,2] %>% sort %>% unique
  step = dst[-1] - dst[-length(dst)]
  lg_step = which(step > d)
  lg_step = lg_step[ lg_step/length(step) < .5 ]

  if(length(lg_step) == 0) return(las)

  lg_step = 1:max(lg_step)

  id_step = which(nn$nn.dists[,2] %in% dst[lg_step])
  id_mat = nn$nn.idx[id_step,]

  if(nrow(id_mat) == 0) return(las)

  for(rid in 1:nrow(id_mat)){
    idx = id_mat[rid,]
    if(is_data_table){
      las[TreeID == nxy$TreeID[idx[2]],]$TreeID = nxy$TreeID[idx[1]]
    }else{
      las@data[TreeID == nxy$TreeID[idx[2]],]$TreeID = nxy$TreeID[idx[1]]
    }
  }

  return(las)
}


#' Classify areas of influence for individual trees
#' @description Identify the \code{TreeID} of every point on a \code{LAS} object based on coordinates from a \code{\link{treeMap}}. Tree region segmentation methods are prefixed by \code{trp}.
#' @template param-las
#' @param map output from \code{\link{treeMap}} or \code{\link{treeMap.positions}}.
#' @param method segmentation algorithm - currently available: \code{\link{trp.voronoi}} and \code{\link{trp.crop}}.
#' @template return-las
#' @examples
#' file = system.file("extdata", "pine_plot.laz", package="TreeLS")
#' tls = readTLS(file)
#' @export
treePoints = function(las, map, method = trp.voronoi()){

  isLAS(las)

  if(hasAttribute(map, 'tree_map_dt')){
    # map = map
  }else if(hasAttribute(map, 'tree_map')){
    map %<>% treeMap.positions(F)
  }else{
    stop('map is not a tree_map object: check ?treeMap')
  }

  if( map$TreeID %>% duplicated %>% any )
    stop('input map must have unique TreeIDs')

  if(!hasAttribute(method, 'tpt_mtd'))
    stop('invalid method: check ?treePoints')

  las %<>% method(map)
  return(las)

}


#' Stem points classification
#' @description Classify stem points of all trees in a \strong{normalized} point cloud. Stem denoising methods are prefixed by \code{stm}.
#' @template param-las
#' @param method stem denoising algorithm - currently available: \code{\link{stm.hough}}, \code{\link{stm.eigen.knn}} and \code{\link{stm.eigen.voxel}}.
#' @return signed \code{LAS} object.
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
stemPoints = function(las, method = stm.hough()){

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

  las = method(las)
  las@data[is.na(Stem)]$Stem = FALSE
  las@data$Stem %<>% as.logical

  return(las)

}


#' Stem segmentation
#' @description Measure stem segments from a classified point cloud. Stem segmentation methods are prefixed by \code{sgt}.
#' @template param-las
#' @param method stem segmentation algorithm - currently available: \code{\link{sgt.ransac.circle}}, \code{\link{sgt.ransac.cylinder}}, \code{\link{sgt.irls.circle}} and \code{\link{sgt.irls.cylinder}}.
#' @return \code{data.table} of stem segments with \emph{stem_dt} signature.
#' @template example-segmentation
#' @export
stemSegmentation = function(las, method=sgt.ransac.circle()){

  isLAS(las)

  if(!hasAttribute(method, 'stem_sgmt_mtd'))
    stop('invalid method: check ?stemSegmentation')

  return(method(las))

}


#' Filter points based on gpstime
#' @description This is a simple wrapper to \code{\link[lidR:filter_poi]{filter_poi}} that takes as inputs relative values instead of absolute time stamp values for filtering a point cloud object based on the gpstime field. This function is particularly useful to check narrow intervals of point cloud frames from mobile scanning data.
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
  las %<>% filter_poi(gpstime >= qts[1] & gpstime <= qts[2])

  return(las)

}


#' Rotate point cloud to fit a horizontal plane
#' @description Check for ground points and rotates the point cloud aligning its ground surface to a horizontal plane (XY).
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

  ground = las@data[,c('X','Y','Z')] %>%
    toLAS %>%
    classify_ground(csf(class_threshold = .2), F) %>%
    filter_poi(Classification == 2)

  center = ground@data[,.(median(X), median(Y))] %>% as.double
  ground = tlsCrop(ground, center[1], center[2], 5) %>% las2xyz

  az = planeAngle(ground, 'z')
  ax = planeAngle(ground, 'x')
  ay = planeAngle(ground, 'y')

  rz = ifelse(az > pi/2, pi-az, -az)
  rx = ifelse(ay < pi/2, -ax, ax)

  rot = rotationMatrix(0, rz, rx) %>% as.matrix
  xy_back = rotationMatrix(0,0,-rx) %>% as.matrix

  minXYZ = apply(las@data[,1:3], 2, min) %>% as.double

  las@data$X = las@data$X - minXYZ[1]
  las@data$Y = las@data$Y - minXYZ[2]
  las@data$Z = las@data$Z - minXYZ[3]

  las@data[,c('X','Y','Z')] = (las2xyz(las) %*% rot) %*% xy_back %>% as.data.table

  las@data$X = las@data$X + minXYZ[1]
  las@data$Y = las@data$Y + minXYZ[2]
  las@data$Z = las@data$Z + minXYZ[3]

  return(las)
}


#' Alter point cloud's coordinates
#' @description Apply transformations to the XYZ axes of a point cloud.
#' @template param-las
#' @param xyz \code{character} vector of length 3 - \code{las}' columns to be reassigned as XYZ, respectively.
#' Use minus signs to mirror an axis` coordinates - more details in the sections below.
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
tlsTransform = function(las, xyz = c('X', 'Y', 'Z'), bring_to_origin = FALSE, rotate = FALSE){

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

  if(rotate){
    las %<>% tlsRotate()
  }

  if(bring_to_origin){
    las = LAS(las@data)

    mincoords = las@data[,.(min(X), min(Y), min(Z))] %>% as.double

    las@data$X = las@data$X - mincoords[1]
    las@data$Y = las@data$Y - mincoords[2]
    las@data$Z = las@data$Z - mincoords[3]
  }

  return(las)

}


#' Point cloud circle fit
#' @description Fits a 2D horizontally-aligned circle on a set of 3D points.
#' @template param-las
#' @param method method for estimating the circle parameters. Currently available: \code{"qr"}, \code{"nm"}, \code{"irls"} and \code{"ransac"}.
#' @template param-n-ransac
#' @template param-inliers
#' @template param-conf
#' @param n_best \code{numeric} - Applicable when \code{method == "ransac"}. For \code{n_best == 0}, takes the best RANSAC iteration as estimates, otherwise, it estimates the parameters as the average of \code{n_best} RANSAC iterations.
circleFit = function(las, method = 'irls', n=5, inliers=.8, p=.99, n_best = 0){
  if(nrow(las@data) < 3) return(NULL)
  if(method == 'ransac' & nrow(las@data) <= n) method = 'qr'
  pars = cppCircleFit(las %>% las2xyz, method, n, p, inliers, n_best)
  pars[3] = pars[3]
  names(pars)[1:4] = c('X','Y','radius', 'err')
  if(length(pars) == 5) names(pars)[5] = 'err2'
  pars = pars %>% t %>% as.data.frame
  return(pars)
}


#' Point cloud cylinder fit
#' @description Fits a cylinder on a set of 3D points.
#' @template param-las
#' @param method method for estimating the cylinder parameters. Currently available: \code{"nm"}, \code{"irls"}, \code{"ransac"} and \code{"bf"}.
#' @template param-n-ransac
#' @template param-inliers
#' @template param-conf
#' @param max_angle \code{numeric} - used when \code{method == "bf"}. The maximum tolerated deviation, in degrees, from an absolute vertical line (Z = c(0,0,1)).
cylinderFit = function(las, method = 'ransac', n=5, inliers=.9, p=.95, max_angle=30, n_best=20){
  if(nrow(las@data) < 3) return(NULL)
  if(method == 'ransac' & nrow(las@data) <= n) method = 'nm'
  pars = cppCylinderFit(las %>% las2xyz, method, n, p, inliers, max_angle, n_best)
  if(method == 'bf'){
    pars[3] = pars[3]
    names(pars) = c('x','y','radius', 'err', 'ax', 'ay')
  }else{
    pars[5] = pars[5]
    pars %<>% c(apply(las@data[,.(X,Y,Z)], 2, function(x) sum(range(x))/2) %>% as.double)
    names(pars) = c('rho','theta','phi', 'alpha', 'radius', 'err', 'px', 'py', 'pz')
  }
  pars = pars %>% t %>% as.data.frame
  return(pars)
}


#' Point cloud cylinder/circle fit
#' @description Fits a 3D cylinder or 2D circle on a set of 3D points, retrieving the optimized parameters.
#' @param stem_segment \code{NULL} or a \code{\link[lidR:LAS]{LAS}} object with a single stem segment. When \code{NULL} returns a parameterized function to be used as input in other functions (e.g. \code{\link{tlsInventory}}).
#' @param shape \code{character} - either \code{circle} or \code{cylinder}.
#' @param algorithm optimization method for estimating the shape parameters. Currently available: \code{"ransac"}, \code{"irls"}, \code{"nm"}, \code{"qr"} (circle only) ,\code{"bf"} (cylinder only).
#' @template param-n-ransac
#' @template param-conf
#' @template param-inliers
#' @template param-n-best
#' @template param-z-dev
#' @export
shapeFit = function(stem_segment=NULL, shape='circle', algorithm='ransac', n=10, conf=0.95, inliers=0.9, n_best = 10, z_dev = 30){

  ls_shapes = c('circle', 'cylinder')
  ls_algo   = c('ransac', 'irls', 'nm', 'qr', 'bf')

  if(!(shape %in% ls_shapes)){
    stop(glue::glue("'{shape}' not available. Enter a valid shape name."))
  }

  if(!(algorithm %in% ls_algo)){
    stop(glue::glue("'{algorithm}' not available. Enter a valid algorithm name."))
  }

  if(algorithm == 'qr' && shape == 'cylinder'){
    stop('qr algorithm only available for circle shapes')
  }else if(algorithm == 'bf' && shape == 'circle'){
    stop('bf algorithm only available for cylinder shapes')
  }

  params = list(
    n = n,
    conf = conf,
    inliers = inliers,
    n_best = n_best,
    z_dev = z_dev
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

  if(shape == 'circle'){
    func = function(las){
      fit = circleFit(las, algorithm, n, inliers, conf, n_best)
      fit = fit[1:4]
      colnames(fit) = c('X', 'Y', 'Radius', 'Error')
      return(as.data.table(fit))
    }
  }else if(shape == 'cylinder'){
    func = function(las){
      fit = cylinderFit(las, algorithm, n, inliers, conf, z_dev, n_best)
      if(algorithm=='bf'){
        colnames(fit) = c('X', 'Y', 'Radius', 'Error', 'DX', 'DY')
      } else {
        colnames(fit)[5:9] = c('Radius', 'Error', 'PX', 'PY', 'PZ')
      }
      return(as.data.table(fit))
    }
  }

  func %<>% setAttribute('shape_fit_method')

  if(!is.null(stem_segment)){
    isLAS(stem_segment)
    return(func(stem_segment))
  }else{
    return(func)
  }
}


#' Extract forest inventory metrics from a point cloud
#' @description Estimation of diameter and height tree-wise.
#' @param las normalized \code{\link[lidR:LAS]{LAS}} object - single tree or forest plot with stem points marked by \code{\link{stemPoints}}.
#' @param dh \code{numeric} - height layer (above ground) from which to estimate stem diameters.
#' @param dw \code{numeric} - height layer width.
#' @param hp \code{numeric} - height percentile to extract per tree (0-1). Use 1 for top height.
#' @param d_method parameterized \code{\link(shapeFit)} function - method to use for diameter estimation.
#' @export
tlsInventory = function(las, dh = 1.3, dw = 0.5, hp = 1, d_method = shapeFit(shape = 'circle', algorithm='ransac', n=15, n_best = 20)){

  preCheck(las)

  if(!hasField(las, 'Stem')){
    stop("las must be a normalized point cloud with highlighted stem points ('Stem' field) - see ?stemPoints")
  }

  params = list(
    dh = dh,
    dw = dw,
    hp = hp
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

  if(!hasAttribute(d_method, 'shape_fit_method')){
    stop('d_method must be a parameterized shapeFit function')
  }

  hfunc = function(Z,p) as.double(quantile(Z, p))
  dfunc = function(X,Y,Z) d_method(suppressMessages(LAS(data.table(X,Y,Z))))

  dlas = filter_poi(las, Stem & Z > (dh - dw/2) & Z < (dh + dw/2))

  if(hasField(las, 'TreeID')){
    h = las@data[TreeID > 0, .(H = hfunc(Z, hp)), by='TreeID']
    d = dlas@data[,dfunc(X,Y,Z),by=TreeID]
    dh_tab = merge(d,h,by='TreeID')[order(TreeID)]
  } else {
    dh_tab = dlas@data %$% dfunc(X,Y,Z)
    dh_tab$H= hfunc(las$Z, hp)
  }

  dh_tab$h_radius = dh
  dh_tab %<>% setAttribute('tls_inventory_dt')
  return(dh_tab)

}


#' Plot \emph{TreeLS} outputs
#' @description Plot the \code{LAS} outputs of tls functions on the same scene using \code{rgl}. Check ?stemSegmentation
#' for usage examples.
#' @param las \code{LAS} object - ideally an output from \code{\link{stemPoints}}.
#' @param sgmt optional \code{data.table} - output from \code{\link{stemSegmentation}}.
#' @param map optional \code{LAS} object - output from \code{\link{treeMap}}.
#' @param tree_id optional \code{numeric} - single \emph{TreeID} to extract from \code{las}.
#' @param color optional - color of the plotted stem segment representations.
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
tlsPlot = function(..., fast=FALSE, tree_id = NULL, segment = NULL){

  tls = as.list(match.call(expand.dots = F))[['...']]

  for(i in 1:length(tls)){
    obj = tls[[i]] %>% as.character %>% get
    tp = class(obj)[1]
    if(!(tp %in% c('LAS', 'data.table', 'data.frame'))){
      stop('input data must be either LAS or data.table objects.')
    }
    tls[[i]] = obj
  }

  seg_ids = FALSE
  tree_ids = FALSE
  pt_cex = 1.5

  tid_plot = function(obj){
    if(tree_ids || !is.null(tree_id)) return(NULL)
    if(!hasField(obj, 'TreeID')) return(NULL)
    add_treeIDs(0, obj, color='yellow')
    tree_ids <<- TRUE
  }

  sid_plot = function(obj){
    if(tree_ids || seg_ids) return(NULL)
    add_segmentIDs(0, obj, color='yellow', pos=4)
    seg_ids <<- TRUE
  }

  h_plot = function(las){
    colors = lidR:::set.colors(las$Z, lidR::height.colors(50))
    rgl.points(las %>% las2xyz, color=colors, size=pt_cex)
  }

  open3d()
  bg3d('black')

  for(las in tls){
    if(class(las)[1] == 'LAS'){

      if(hasField(las, 'TreeID') && !is.null(tree_id)){
        las = filter_poi(las, TreeID == tree_id)
      }

      if(hasField(las, 'Segment') && !is.null(segment)){
        las = filter_poi(las, Segment == segment)
      }

      if(hasField(las, 'Stem')){
        add_stemPoints(0, las, color='white', size=pt_cex)
        if(hasField(las, 'TreeID')) add_treePoints(0, las, size=pt_cex) else h_plot(las)
        tid_plot(las)
        sid_plot(las)
      }else if(hasAttribute(las, 'tree_map')){
        add_treeMap(0, las, color='yellow')
        tid_plot(las)
      }else if(hasField(las, 'TreeID')){
        add_treePoints(0, las, size=pt_cex)
        tid_plot(las)
      }else{
        h_plot(las)
      }
      pt_cex = pt_cex + .5
    }else{

      if(hasField(las, 'TreeID') && !is.null(tree_id)){
        las = las[TreeID == tree_id]
      }

      if(hasField(las, 'Segment') && !is.null(segment)){
        las = las[Segment == segment]
      }

      if(hasAttribute(las, 'tls_inventory_dt')){
        add_tlsInventory(0, las, fast=fast)
        tid_plot(las)
      }else if(hasAttribute(las, 'single_stem_dt') || hasAttribute(las, 'multiple_stems_dt')){
        add_stemSegments(0, las, fast=fast)
        tid_plot(las)
        sid_plot(las)
      }else if(hasAttribute(las, 'tree_map_dt')){
        add_treeMap(0, las, color='yellow')
        tid_plot(las)
      }
    }
  }

  pan3d(2)

}


#########################################


#' EXPERIMENTAL: Point cloud multiple circle fit
#' @description Searches and fits circles multiple 2D circles on a point cloud layer from a single tree.
#' @param dlas \code{\link[lidR:LAS]{LAS}} object.
#' @template param-pixel-size
#' @template param-max-radius
#' @param votes_percentile \code{numeric} - votes filter criterion for isolating nearby stems.
#' @template param-min-density
#' @param plot \code{logical} - plot the results?
#' @export
shapeFit.forks = function(dlas, pixel_size = .02, max_radius = .2, votes_percentile = .7, min_density = .25, plot=FALSE){

  isLAS(dlas)

  params = list(
    pixel_size = pixel_size,
    max_radius = max_radius,
    votes_percentile = votes_percentile,
    min_density = min_density
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

  hg = getHoughCircle(dlas %>% las2xyz, pixel_size, rad_max = max_radius, min_den = min_density, min_votes = 2) %>% do.call(what=rbind) %>% as.data.table
  names(hg) = c('x','y','r','v')
  hg = hg[v > quantile(v, votes_percentile)]
  hg$clt = 1

  hough_clusters = hg
  centers = hg[v == max(v)] %>% apply(2,mean) %>% t %>% as.data.table

  k = 2
  repeat{
    km = kmeans(hg[,1:3], k)
    hg$clt = km$cluster
    mxs = hg[,.(x=mean(x), y=mean(y), r=mean(r), v=mean(v)), by=clt]
    dst = mxs[,c('x','y')] %>% dist
    combs = combn(nrow(mxs), 2)
    rst = apply(combs, 2, function(x){mxs$r[x[1]] + mxs$r[x[2]]}) - pixel_size
    is_forked = all(min(dst) > rst)

    if(!is_forked) break

    hough_clusters$clt = km$cluster
    centers = mxs[,.(x=mean(x),y=mean(y),r=mean(r),v=mean(v)),by=clt][order(clt)]
    centers$ssRatio = (km$withinss / km$size) / min(km$withinss / km$size)
    centers$nRatio = km$size / max(km$size)
    centers$absRatio = centers$ssRatio / centers$nRatio
    v_ratio = centers %$% ifelse(min(absRatio) > 5, absRatio, 5)
    centers = centers[absRatio <= v_ratio][order(-v)]
    k=k+1
  }

  if(plot){
    plot(dlas$Y ~ dlas$X, cex=.5, asp=1, pch=20, ylab='Y (m)', xlab='X (m)')#, ...)
    vcols = lidR:::set.colors(hough_clusters$v, height.colors(hough_clusters$v %>% unique %>% length))
    vcols = lidR:::set.colors(hough_clusters$clt, height.colors(hough_clusters$clt %>% unique %>% length))
    points(hough_clusters$x, hough_clusters$y, col=vcols, pch=20, cex=1)
  }

  estimates = data.frame()
  dlas@data$StemID = 0
  for(i in 1:nrow(centers)){
    temp = centers[i]

    dsts = dlas@data %$% sqrt( (X - temp$x)^2 + (Y - temp$y)^2 )
    pts = dsts < temp$r + 2*pixel_size & dlas@data$StemID == 0
    dlas@data$StemID[pts] = i

    cld = filter_poi(dlas, StemID == i)
    est = circleFit(cld, 'ransac', n=15, inliers = .7, n_best = 30)
    if(is.null(est)) next

    dst = cld@data %$% sqrt( (X - est$X)^2 + (Y - est$Y)^2 )
    rad = est$d / 200
    thickness = ifelse(rad/2 > pixel_size, pixel_size,  rad/2)

    inner = dst <= thickness
    area_inner = pi*thickness^2
    den_inner = (inner %>% which %>% length) / area_inner

    strip = dst > (rad - thickness/2) & dst < (rad + thickness/2)
    area_strip = pi*((rad + thickness/2)^2 - (rad - thickness/2)^2)
    den_strip = (strip %>% which %>% length) / area_strip

    outer = dst < (rad + thickness) & dst >= (rad + thickness/2)
    area_outer = pi*((rad + thickness)^2 - (rad + thickness/2)^2)
    den_outer = (outer %>% which %>% length) / area_outer

    est$ratioInner = den_inner / den_strip
    est$ratioOuter = den_outer / den_strip
    est$StemID = i
    est$TreeID = dlas$TreeID[1]
    est$score = 1

    if(est$ratioInner > 1){
      est$score = 3
    }else if(est$ratioInner > .8 | est$ratioOuter > .8){
      est$score = 2
    }

    estimates %<>% rbind(est)

    if(plot){
      cld@data %$% points(X,Y, col=ifelse(strip, 'darkblue', ifelse(inner | outer, 'darkred', 'black')), cex=.75, pch=20)

      angs = seq(0,pi*2,length.out = 36)
      xh = temp$x + cos(angs) * temp$r
      yh = temp$y + sin(angs) * temp$r
      lines(xh,yh,col='orange',lwd=2)

      xr = est$X + cos(angs) * est$radius
      yr = est$Y + sin(angs) * est$radius
      lines(xr,yr,col='green',lwd=2)
    }
  }

  check = estimates$score == 3
  if(!all(check)){
    estimates = estimates[ estimates$score < 3 ,]
  }

  return(estimates)
}

