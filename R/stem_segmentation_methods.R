#' @param tol \code{numeric} - tolerance offset between absolute radii estimates and hough transform estimates.
#' @param n \code{integer} - number of points selected on every RANSAC iteration.
#' @param conf \code{numeric} - confidence level.
#' @param inliers \code{numeric} - expected proportion of inliers among stem segments' point cloud chunks.
sgmt.ransac.circle = function(tol=0.025, n = 10, conf = 0.99, inliers = 0.8){

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

  func = function(las){

    if(las %>% hasAttribute('single_stem_points')){

      message('performing single stem segmentation')

      las %<>% lasfilter(Stem)

      estimates = ransacStem(las %>% las2xyz, las@data$Segment, las@data$Radius, n, conf, inliers, tol) %>% do.call(what = rbind) %>% as.data.table
      names(estimates) = c('X', 'Y', 'Radius', 'Error', 'Segment')

      z = las@data[, .(AvgHeight = mean(Z), N = .N), Segment]

      estimates %<>% merge(z, by='Segment')

      # Segment = NULL
      estimates = estimates[order(Segment)]

      estimates %<>% setAttribute("single_stem_dt")

    }else if(las %>% hasAttribute('multiple_stem_points')){

      message('performing multiple stems segmentation')

      las %<>% lasfilter(Stem)

      estimates = ransacPlot(las %>% las2xyz, las$TreeID, las$Segment, las$Radius, n, conf, inliers, tol) %>% sapply(do.call, what=rbind) %>% do.call(what = rbind) %>% as.data.table
      names(estimates) = c('X', 'Y', 'Radius', 'Error', 'Segment', 'TreeID')

      z = las@data[, .(AvgHeight = mean(Z), N = .N), .(TreeID, Segment)]

      estimates %<>% merge(z, by=c('TreeID', 'Segment'))

      # TreeID = Segment = NULL
      estimates = estimates[order(TreeID, Segment)]

      estimates %<>% setAttribute("multiple_stems_dt")

    }else{
      stop('stem points identifier missing - check ?stemPoints')
    }

    return(estimates)
  }

  func %>% setAttribute('stem_sgmt_mtd') %>% return()

}
