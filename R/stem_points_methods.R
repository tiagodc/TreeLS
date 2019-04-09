#' Stem denoising algorithm: Hough Transform
#' @description This function is meant to be used inside \code{\link{stemPoints}}. It applies an adapted version of the Hough Transform for circle search. Mode details are given in the sections below.
#' @template param-hstep
#' @template param-max-radius
#' @template param-hbase
#' @template param-pixel-size
#' @template param-min-density
#' @template param-min-votes
#' @section \code{LAS@data} Special Fields:
#'
#' Meaninful fields in the output:
#'
#' \itemize{
#' \item \code{TreeID}: unique tree ID of the point - available when a \emph{tree_map} is provided
#' \item \code{Stem}: \code{TRUE} for stem points
#' \item \code{Segment}: stem segment number (from bottom to top)
#' \item \code{Radius}: approximate radius of the point's stem segment estimated by the Hough Transform - always a multiple of the \code{pixel_size}
#' \item \code{Votes}: votes received by the stem segment's center through the Hough Transform
#' }
#'
#' @template section-hough-transform
#' @template reference-thesis
#' @examples
#' file = system.file("extdata", "spruce.laz", package="TreeLS")
#' tls = readTLS(file)
#'
#' ### identify stem points
#' tls = stemPoints(tls, method = stem.hough(max_radius=.2))
#' plot(tls, color='Stem')
#'
#' @export
stem.hough = function(hstep=0.5, max_radius=0.25, hbase = c(1,2.5), pixel_size=0.025, min_density=0.1, min_votes=3){

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

  func = function(las, map){

    if(min(las$Z) < 0)
      message("points with Z below 0 will be ignored")

    if(min(las$Z) > 5)
      message("point cloud doesn't look normalized (Z values too high) - check ?tlsNormalize")

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

      if(!hasAttribute(las, 'tree_points')){
        las@data$TreeID[!groundPts] = results$TreeID
        las@data$TreeID[las@data$TreeID %>% is.na] = 0
      }
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

    if(map %>% is.null){
      las %<>% setAttribute("single_stem_points")
    }else{
      las %<>% setAttribute("multiple_stem_points")
    }

    return(las)

  }

  func %<>% setAttribute('stem_pts_mtd')

  return(func)

}
