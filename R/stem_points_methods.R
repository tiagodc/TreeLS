#' @template param-hstep
#' @template param-max-radius
#' @template param-hbase
#' @template param-pixel-size
#' @template param-min-density
#' @template param-min-votes
#' @template section-hough-transform
#' @template reference-thesis
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
