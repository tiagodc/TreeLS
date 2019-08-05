#' Classify point cloud according to voronoi polygons.
#' @description Performs point cloud segmentation according to regions of influence based on voronoi polygons, calculated from a tree map object.
#' @examples
#' file = system.file("extdata", "pine.laz", package="TreeLS")
#' tls = readTLS(file)
#' @export
trp.voronoi = function(){

  func = function(las, xymap){
    xt = extent(las) + c(-1,1,-1,1)
    vPoly = dismo::voronoi(xymap[,2:3], xt)
    vPoly$id = xymap$TreeID
    names(vPoly) = 'TreeID'
    crs(vPoly) = crs(las)

    xysp = las@data[,.(X,Y)] %>% SpatialPoints %>% over(y = vPoly)
    xysp = xysp[,1] %>% as.double
    xysp[is.na(xysp)] = 0

    las@data$TreeID = xysp
    las@data$TreeID[las@data$Classification == 2] = 0

    las %<>% setAttribute('tree_points')

    return(las)
  }

  func %<>% setAttribute('tpt_mtd')
  return(func)
}

trp.clip = function(r = 1, circle=T){
  func = function(las, xymap){

    las@data$TreeID = 0
    for(i in 1:nrow(xymap)){
      x = xymap$X[i]
      y = xymap$Y[i]
      id = xymap$TreeID[i]
      bool = RCropCloud(las %>% las2xyz, x, y, r, circle, F)
      las@data$TreeID[bool] = id
    }
    las %<>% setAttribute('tree_points')
    return(las)
  }

  func %<>% setAttribute('tpt_mtd')
  return(func)
}
