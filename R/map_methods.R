#' @template param-hmin-hmax
#' @template param-hstep
#' @template param-pixel-size
#' @template param-max-radius
#' @template param-min-density
#' @template param-min-votes
#' @template reference-thesis
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
