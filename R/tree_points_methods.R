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

trp.crop = function(l = 1, circle=T){
  func = function(las, xymap){
    las@data$TreeID = treeIdsFromMap(las@data[,.(X,Y)] %>% as.matrix, xymap[,.(X,Y)] %>% as.matrix, xymap$TreeID %>% as.integer, l, circle)
    las %<>% setAttribute('tree_points')
    return(las)
  }
  func %<>% setAttribute('tpt_mtd')
  return(func)
}
