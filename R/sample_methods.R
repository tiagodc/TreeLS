#' @useDynLib TreeLS, .registration = TRUE
voxelize = function(spacing = 0.05){

  if(spacing <= 0)
    stop('spacing must be a positive number')

  func = function(las){
    las %>% las2xyz %>% thinCloud(spacing) %>% return()
  }

  func %<>% setAttribute('tls_sample_mtd')

  return(func)
}

#' @importFrom stats rbinom
randomize = function(p = 0.5){

  if(p <= 0)
    stop('p must be a positive number')

  if(p >= 1)
    stop('p must be a number between 0 and 1 for random sampling')

  func = function(las){
    n = nrow(las@data)
    keep = rbinom(n, 1, p) == 1
    return(keep)
  }

  func %<>% setAttribute('tls_sample_mtd')

  return(func)
}
