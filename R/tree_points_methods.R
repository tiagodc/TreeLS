#' Tree points algorithm: assign tree points by voronoi polygons.
#' @description This function is meant to be used inside \code{\link{treePoints}}. Assign all points to a \emph{TreeID}, based on the closest \code{\link{treeMap}} coordinate.
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


#' Tree points algorithm: assign tree points by fixed area.
#' @description This function is meant to be used inside \code{\link{treePoints}}. Assign points to a \emph{TreeID} inside circles/squares of fixed area around \code{\link{treeMap}} coordinates.
#' @param l \code{numeric} - circle radius (for circle == TRUE), or square side length (for circle == FALSE).
#' @param circle \code{logical} - assign \emph{TreeID}s inside circles (TRUE) or squares (FALSE)?
#' @export
trp.crop = function(l = 1, circle=T){
  func = function(las, xymap){
    las@data$TreeID = treeIdsFromMap(las@data[,.(X,Y)] %>% as.matrix, xymap[,.(X,Y)] %>% as.matrix, xymap$TreeID %>% as.integer, l, circle)
    las %<>% setAttribute('tree_points')
    return(las)
  }
  func %<>% setAttribute('tpt_mtd')
  return(func)
}
