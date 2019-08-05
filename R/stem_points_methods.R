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
#' Meaninful new fields in the output:
#'
#' \itemize{
#' \item \code{Stem}: \code{TRUE} for stem points
#' \item \code{Segment}: stem segment number (from bottom to top and nested with TreeID)
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
#' tls = stemPoints(tls, method = stm.hough(max_radius=.2))
#' plot(tls, color='Stem')
#'
#' @export
stm.hough = function(hstep=0.5, max_radius=0.25, hbase = c(1,2.5), pixel_size=0.025, min_density=0.1, min_votes=3){

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

  func = function(las){

    if(min(las$Z) < 0)
      message("points with Z below 0 will be ignored")

    if(min(las$Z) > 5)
      message("point cloud doesn't look normalized (Z values too high) - check ?tlsNormalize")

    groundPts = if(las %>% hasField('Classification')){
      las$Classification == 2
    }else{
      rep(F, las@data %>% nrow)
    }

    if(!hasField(las, 'TreeID')){
      message('no TreeID field found with tree_points signature: performing single stem point classification')
      results = houghStemPoints(las2xyz(las)[!groundPts,], hbase[1], hbase[2], hstep, max_radius, pixel_size, min_density, min_votes)
    }else{
      message('performing point classification on multiple stems')
      results = houghStemPlot(las2xyz(las)[!groundPts,], las@data$TreeID[!groundPts], hbase[1], hbase[2], hstep, max_radius, pixel_size, min_density, min_votes)
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

    if(!hasField(las, 'TreeID')){
      las %<>% setAttribute("single_stem_points")
    }else{
      las %<>% setAttribute("multiple_stem_points")
    }

    return(las)

  }

  func %<>% setAttribute('stem_pts_mtd')

  return(func)

}


stm.eigen.knn = function(hstep = .5, pln = .2, vrt = 20, vxl = .025, max_d = .5, dvt = .2, v3d = F){

  func = function(las){

    las@data$PointID = 1:nrow(las@data)
    mtrlst = c('N', 'Planarity', 'Verticality', 'EigenVectors')

    # if(hasField(las, 'TreeID')){
    #   las@data = las@data %>% split(las@data$TreeID) %>% lapply(LAS) %>%
    #     lapply(pointMetrics, method = ptm.knn(), metrics_list=mtrlst) %>%
    #     lapply(function(x) x@data) %>% do.call(what=rbind) %>% as.data.table
    # }

    las = pointMetrics(las, ptm.knn(), mtrlst)

    las@data$Stem = with(las@data, Classification != 2 & N > 3 & Planarity < pln & abs(Verticality - 90) < vrt)

    stemSeg = seq(0, max(las$Z)+hstep, hstep)
    las@data$Segment = cut(las$Z, stemSeg, include_lowest=T, right=F, ordered_result=T) %>% as.integer
    las@data$Segment[las@data$Z < 0] = 0

    if(hasField(las, 'TreeID')){
      points = las@data[Stem & order(TreeID, Segment, PointID), .(TreeID, Segment, PointID, X, Y, Z, EigenVector13, EigenVector23, EigenVector33)]

      tds = points$TreeID
      sgs = points$Segment
      ids = points$PointID

      votes = points[,-c(1:3)] %>% as.matrix %>% plotEigenHough(ids, tds, sgs, vxl, max_d/2, !v3d, F) %>%
        lapply(function(x) x %>% do.call(what=cbind)) %>% do.call(what=rbind) %>% as.data.table

      colnames(votes) = c('Votes','Radius','PointID', 'Segment', 'TreeID')

    }else{
      points = las@data[Stem & order(Segment, PointID), .(Segment, PointID, X, Y, Z, EigenVector13, EigenVector23, EigenVector33)]

      sgs = points$Segment
      ids = points$PointID

      votes = points[,-c(1:2)] %>% as.matrix %>% treeEigenHough(ids, sgs, vxl, max_d/2, !v3d, F) %>%
        lapply(function(x) x %>% do.call(what=cbind)) %>% do.call(what=rbind) %>% as.data.table

      colnames(votes) = c('Votes','Radius','PointID', 'Segment')
    }

    keepcols = !(colnames(las@data) %in% c('Votes', 'Radius', 'MaxVotes'))
    keepcols = colnames(las@data)[keepcols]
    las@data = las@data[,..keepcols]
    las@data = merge(las@data, votes[,.(PointID, Votes, Radius)], by='PointID', sort=F, all.x=T)
    las@data[!las@data$Stem, c('Votes', 'Radius')] = 0

    if(hasField(las, 'TreeID')){
      maxVotes = las@data[,.(MaxVotes = max(Votes)), by=TreeID]
      las@data = merge(las@data, maxVotes, by='TreeID', sort=F, all.x=T)
      las@data$VotesWeight = las@data$Votes / las@data$MaxVotes
      las@data$MaxVotes = NULL
    }else{
      las@data$VotesWeight = las@data$Votes / max(las@data$Votes)
    }

    las@data$Stem = las@data$VotesWeight > dvt

    return(las)
  }

  func %<>% setAttribute('stem_pts_mtd')
  return(func)
}


stm.eigen.voxel = function(hstep = .5, pln = .2, vrt = 20, vxl = .1, max_d = .5, dvt = .2, v3d = F){

  func = function(las){

    mtrlst = c('N', 'Planarity', 'Verticality', 'EigenVectors')

    las = pointMetrics(las, ptm.voxels(vxl), mtrlst)

    las@data$Stem = with(las@data, Classification != 2 & N > 3 & Planarity < pln & abs(Verticality - 90) < vrt)

    stemSeg = seq(0, max(las$Z)+hstep, hstep)
    las@data$Segment = cut(las$Z, stemSeg, include_lowest=T, right=F, ordered_result=T) %>% as.integer
    las@data$Segment[las@data$Z < 0] = 0

    if(hasField(las, 'TreeID')){
      voxels = las@data[Stem & order(TreeID, Segment), .(TreeID = mean(TreeID), Segment = mean(Segment), X = mean(X), Y = mean(Y), Z = mean(Z), EigenVector13 = mean(EigenVector13), EigenVector23 = mean(EigenVector23), EigenVector33 = mean(EigenVector33)), by = VoxelID]

      tds = voxels$TreeID
      sgs = voxels$Segment
      ids = voxels$VoxelID

      votes = voxels[,-c(1:3)] %>% as.matrix %>% plotEigenHough(ids, tds, sgs, vxl, max_d/2, !v3d, F) %>%
        lapply(function(x) x %>% do.call(what=cbind)) %>% do.call(what=rbind) %>% as.data.table

      colnames(votes) = c('Votes','Radius','VoxelID', 'Segment', 'TreeID')

    }else{
      voxels = las@data[Stem & order(Segment), .(Segment = mean(Segment), X = mean(X), Y = mean(Y), Z = mean(Z), EigenVector13 = mean(EigenVector13), EigenVector23 = mean(EigenVector23), EigenVector33 = mean(EigenVector33)), by = VoxelID]

      sgs = voxels$Segment
      ids = voxels$VoxelID

      votes = voxels[,-c(1:2)] %>% as.matrix %>% treeEigenHough(ids, sgs, vxl, max_d/2, !v3d, F) %>%
        lapply(function(x) x %>% do.call(what=cbind)) %>% do.call(what=rbind) %>% as.data.table

      colnames(votes) = c('Votes','Radius','VoxelID', 'Segment')
    }

    keepcols = !(colnames(las@data) %in% c('Votes', 'Radius', 'MaxVotes'))
    keepcols = colnames(las@data)[keepcols]
    las@data = las@data[,..keepcols]
    las@data = merge(las@data, votes[,.(VoxelID, Votes, Radius)], by='VoxelID', sort=F, all.x=T)
    las@data[!las@data$Stem, c('Votes', 'Radius')] = 0

    if(hasField(las, 'TreeID')){
      maxVotes = las@data[,.(MaxVotes = max(Votes)), by=TreeID]
      las@data = merge(las@data, maxVotes, by='TreeID', sort=F, all.x=T)
      las@data$VotesWeight = las@data$Votes / las@data$MaxVotes
      las@data$MaxVotes = NULL
    }else{
      las@data$VotesWeight = las@data$Votes / max(las@data$Votes)
    }

    las@data$Stem = las@data$VotesWeight > dvt

    return(las)
  }

  func %<>% setAttribute('stem_pts_mtd')
  return(func)
}
