# ===============================================================================
#
# Developers:
#
# Tiago de Conto - tdc.florestal@gmail.com -  https://github.com/tiagodc/
#
# COPYRIGHT: Tiago de Conto, 2019
#
# This piece of software is open and free to use, redistribution and modifications
# should be done in accordance to the GNU General Public License >= 3
# Use this software as you wish, but no warranty is provided whatsoever. For any
# comments or questions on TreeLS, please contact the developer (prefereably through my github account)
#
# If publishing any work/study/research that used the tools in TreeLS,
# please don't forget to cite the proper sources!
#
# Enjoy!
#
# ===============================================================================

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
stm.hough = function(h_step=0.5, max_radius=0.25, h_base = c(1,2.5), pixel_size=0.025, min_density=0.1, min_votes=3){

  if(length(h_base) != 2)
    stop('hbase must be a numeric vector of length 2')

  if(diff(hbase) <= 0)
    stop('hbase[2] must be larger than hbase[1]')

  params = list(
    h_step = h_step,
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

    surveyPts = if(las %>% hasField('Classification')){
      las$Classification != 2
    }else{
      rep(T, las@data %>% nrow)
    }

    if(!hasField(las, 'TreeID')){
      message('no TreeID field found with tree_points signature: performing single stem point classification')
      results = houghStemPoints(las2xyz(las)[surveyPts,], hbase[1], hbase[2], h_step, max_radius, pixel_size, min_density, min_votes)
    }else{
      message('performing point classification on multiple stems')
      surveyPts = surveyPts & las$TreeID > 0
      results = houghStemPlot(las2xyz(las)[surveyPts,], las@data$TreeID[surveyPts], hbase[1], hbase[2], h_step, max_radius, pixel_size, min_density, min_votes)
    }

    las@data$Stem = F
    las@data$Stem[surveyPts] = results$Stem

    las@data$Segment = 0
    las@data$Segment[surveyPts] = results$Segment

    las@data$Radius = 0
    las@data$Radius[surveyPts] = results$Radius

    las@data$Votes = 0
    las@data$Votes[surveyPts] = results$Votes

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


#' Stem denoising algorithm: KNN geometry + voxel voting
#' @description This function is meant to be used inside \code{\link{stemPoints}}. It filters points based on their nearest neighborhoods geometries (check \code{\link{pointMetrics}}) and assign them to stem patches if reaching a voxel with enough votes.
#' @template param-h_step
#' @template param-max-planarity
#' @template param-max-verticality
#' @template param-voxel-spacing
#' @template param-max-d
#' @template param-votes-weight
#' @template param-v3d
#' @export
stm.eigen.knn = function(h_step = .5, max_planarity = .2, max_verticality = 20, voxel_spacing = .025, max_d = .5, votes_weight = .2, v3d = F){

  params = list(
    h_step = h_step,
    max_planarity = max_planarity,
    max_verticality = max_verticality,
    voxel_spacing = voxel_spacing,
    max_d = max_d,
    votes_weight = votes_weight
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

  if(max_planarity > 1) stop('max_planarity must be a number between 0 and 1')
  if(max_verticality > 180) stop('max_verticality must be a number between 0 and 180')
  if(votes_weight > 1) stop('votes_weight must be a number between 0 and 1')

  func = function(las){

    las@data$PointID = 1:nrow(las@data)
    mtrlst = c('N', 'Planarity', 'Verticality', 'EigenVectors')
    mtrnames = c(mtrlst[-4], 'EigenVector13', 'EigenVector23', 'EigenVector33')

    # if(hasField(las, 'TreeID')){
    #   las@data = las@data %>% split(las@data$TreeID) %>% lapply(LAS) %>%
    #     lapply(pointMetrics, method = ptm.knn(), metrics_list=mtrlst) %>%
    #     lapply(function(x) x@data) %>% do.call(what=rbind) %>% as.data.table
    # }

    checkPtMetrics = mtrnames %>% sapply(function(x) hasField(las, x)) %>% as.logical %>% all
    if(!checkPtMetrics){
      message('Calculating knn pointMetrics')
      las = pointMetrics(las, ptm.knn(), mtrlst)
    }

    las@data$Stem = with(las@data, N > 3 & Planarity < max_planarity & abs(Verticality - 90) < max_verticality)
    if(hasField(las, 'Classification')){
      las@data$Stem = las@data$Stem & las@data$Classification != 2
    }

    stemSeg = seq(0, max(las$Z)+h_step, h_step)
    las@data$Segment = cut(las$Z, stemSeg, include_lowest=T, right=F, ordered_result=T) %>% as.integer
    las@data$Segment[las@data$Z < 0] = 0

    if(hasField(las, 'TreeID')){
      las@data$Stem = las@data$Stem & las@data$TreeID > 0
      points = las@data[Stem & order(TreeID, Segment, PointID), .(TreeID, Segment, PointID, X, Y, Z, EigenVector13, EigenVector23, EigenVector33)]

      tds = points$TreeID
      sgs = points$Segment
      ids = points$PointID

      votes = points[,-c(1:3)] %>% as.matrix %>% plotEigenHough(ids, tds, sgs, voxel_spacing, max_d/2, !v3d, F) %>%
        lapply(function(x) x %>% do.call(what=cbind)) %>% do.call(what=rbind) %>% as.data.table

      colnames(votes) = c('Votes','Radius','PointID', 'Segment', 'TreeID')

    }else{
      points = las@data[Stem & order(Segment, PointID), .(Segment, PointID, X, Y, Z, EigenVector13, EigenVector23, EigenVector33)]

      sgs = points$Segment
      ids = points$PointID

      votes = points[,-c(1:2)] %>% as.matrix %>% treeEigenHough(ids, sgs, voxel_spacing, max_d/2, !v3d, F) %>%
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

    las@data$Stem = las@data$VotesWeight > votes_weight

    return(las)
  }

  func %<>% setAttribute('stem_pts_mtd')
  return(func)
}


#' Stem denoising algorithm: Voxel geometry + voting
#' @description This function is meant to be used inside \code{\link{stemPoints}}. It filters points based on their voxel geometries (check \code{\link{pointMetrics}}) and assign them to stem patches if reaching a voxel with enough votes.
#' @template param-h_step
#' @template param-max-planarity
#' @template param-max-verticality
#' @template param-voxel-spacing
#' @template param-max-d
#' @template param-votes-weight
#' @template param-v3d
#' @export
stm.eigen.voxel = function(h_step = .5, max_planarity = .2, max_verticality = 20, voxel_spacing = .1, max_d = .5, votes_weight = .2, v3d = F){

  params = list(
    h_step = h_step,
    max_planarity = max_planarity,
    max_verticality = max_verticality,
    voxel_spacing = voxel_spacing,
    max_d = max_d,
    votes_weight = votes_weight
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

  if(max_planarity > 1) stop('max_planarity must be a number between 0 and 1')
  if(max_verticality > 180) stop('max_verticality must be a number between 0 and 180')
  if(votes_weight > 1) stop('votes_weight must be a number between 0 and 1')

  func = function(las){

    mtrlst = c('N', 'Planarity', 'Verticality', 'EigenVectors')
    mtrnames = c(mtrlst[-4], 'EigenVector13', 'EigenVector23', 'EigenVector33', 'VoxelID')

    checkPtMetrics = mtrnames %>% sapply(function(x) hasField(las, x)) %>% as.logical %>% all
    if(!checkPtMetrics){
      message('Calculating voxel pointMetrics')
      las = pointMetrics(las, ptm.voxel(voxel_spacing), mtrlst)
    }

    las@data$Stem = with(las@data, Classification != 2 & N > 3 & Planarity < max_planarity & abs(Verticality - 90) < max_verticality)

    stemSeg = seq(0, max(las$Z)+h_step, h_step)
    las@data$Segment = cut(las$Z, stemSeg, include_lowest=T, right=F, ordered_result=T) %>% as.integer
    las@data$Segment[las@data$Z < 0] = 0

    if(hasField(las, 'TreeID')){
      voxels = las@data[Stem & order(TreeID, Segment), .(TreeID = mean(TreeID), Segment = mean(Segment), X = mean(X), Y = mean(Y), Z = mean(Z), EigenVector13 = mean(EigenVector13), EigenVector23 = mean(EigenVector23), EigenVector33 = mean(EigenVector33)), by = VoxelID]

      tds = voxels$TreeID
      sgs = voxels$Segment
      ids = voxels$VoxelID

      votes = voxels[,-c(1:3)] %>% as.matrix %>% plotEigenHough(ids, tds, sgs, voxel_spacing, max_d/2, !v3d, F) %>%
        lapply(function(x) x %>% do.call(what=cbind)) %>% do.call(what=rbind) %>% as.data.table

      colnames(votes) = c('Votes','Radius','VoxelID', 'Segment', 'TreeID')

    }else{
      voxels = las@data[Stem & order(Segment), .(Segment = mean(Segment), X = mean(X), Y = mean(Y), Z = mean(Z), EigenVector13 = mean(EigenVector13), EigenVector23 = mean(EigenVector23), EigenVector33 = mean(EigenVector33)), by = VoxelID]

      sgs = voxels$Segment
      ids = voxels$VoxelID

      votes = voxels[,-c(1:2)] %>% as.matrix %>% treeEigenHough(ids, sgs, voxel_spacing, max_d/2, !v3d, F) %>%
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

    las@data$Stem = las@data$VotesWeight > votes_weight

    return(las)
  }

  func %<>% setAttribute('stem_pts_mtd')
  return(func)
}
