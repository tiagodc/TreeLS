#' Tree mapping algorithm: Hough Transform
#' @description This function is meant to be used inside \code{\link{treeMap}}. It applies an adapted version of the Hough Transform for circle search. Mode details are given in the sections below.
#' @template param-hmin-hmax
#' @template param-hstep
#' @template param-pixel-size
#' @template param-max-radius
#' @template param-min-density
#' @template param-min-votes
#' @section \code{LAS@data} Special Fields:
#'
#' Each point in the \code{LAS} object output represents a pixel center that is
#' \emph{possibly} also a stem cross-section center.
#'
#' The variables describing each point in the output are:
#'
#' \itemize{
#' \item \code{Intensity}: number of votes received by that point
#' \item \code{PointSourceID}: unique stem segment ID (among all trees)
#' \item \code{Keypoint_flag}: if \code{TRUE}, the point is the most likely circle center
#' of its stem segment (\code{PointSourceID})
#' \item \code{Radii}: approximate radius estimated by that point - always a multiple of the \code{pixel_size}
#' \item \code{TreeID}: unique tree ID of the point
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
#' @export
map.hough = function(hmin = 1, hmax = 3, hstep = 0.5, pixel_size = 0.025, max_radius = 0.25, min_density = 0.1, min_votes = 3){

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

  func = function(las){

    rgz = las$Z %>% range

    if(hmax < rgz[1])
      stop('hmax is too low - below the point cloud')

    if(hmin > rgz[2])
      stop('hmin is too high - above the point cloud')

    map = stackMap(las %>% las2xyz, hmin, hmax, hstep, pixel_size, max_radius, min_density, min_votes) %>%
      do.call(what=cbind) %>% as.data.table

    map$Intensity %<>% as.integer
    map$Keypoint_flag %<>% as.logical
    map$PointSourceID %<>% as.integer
    map$TreePosition %<>% as.logical
    map %<>% LAS %>% setHeaderTLS

    return(map)
  }

  func %<>% setAttribute('tls_map_mtd')

  return(func)
}
