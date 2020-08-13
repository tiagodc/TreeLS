#' Tree points algorithm: voronoi polygons.
#' @description This function is meant to be used inside \code{\link{treePoints}}. Assign **all** points to a \code{TreeID} based on their closest \code{\link{treeMap}} coordinate.
#' @importFrom dismo voronoi
#' @importFrom sp over SpatialPoints
#' @export
trp.voronoi = function(){

  func = function(las, xymap){
    xt = extent(las) + c(-1,1,-1,1)
    v_poly = voronoi(xymap[,2:3], xt)
    v_poly$id = xymap$TreeID
    names(v_poly) = 'TreeID'
    crs(v_poly) = crs(las)

    xysp = las@data[,.(X,Y)] %>% SpatialPoints
    crs(xysp) = crs(las)
    xysp %<>% over(y = v_poly)

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


#' Tree points algorithm: fixed size patches.
#' @description This function is meant to be used inside \code{\link{treePoints}}. Assign points to a \code{TreeID} inside circles/squares of fixed area around \code{\link{treeMap}} coordinates.
#' @param l \code{numeric} - circle radius or square side length.
#' @param circle \code{logical} - assign \code{TreeID}s to circular (\code{TRUE}) or squared (\code{FALSE}) patches.
#' @export
trp.crop = function(l = 1, circle=TRUE){
  func = function(las, xymap){
    las@data$TreeID = treeIdsFromMap(las@data[,.(X,Y)] %>% as.matrix, xymap[,.(X,Y)] %>% as.matrix, xymap$TreeID %>% as.integer, l, circle)
    las %<>% setAttribute('tree_points')
    return(las)
  }
  func %<>% setAttribute('tpt_mtd')
  return(func)
}
