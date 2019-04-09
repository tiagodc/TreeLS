trees.voronoi = function(){

  func = function(las, xymap){
    xt = extent(las) + c(-1,1,-1,1)
    vPoly = dismo::voronoi(xymap[,2:3], xt)
    vPoly$id = xymap$TreeID
    names(vPoly) = 'TreeID'

    las %<>% lasmergespatial(vPoly, 'TreeID')
    las@data$TreeID[las@data$Classification == 2] = 0

    las %<>% setAttribute('tree_points')

    return(las)
  }

  return(func)
}

