# ===============================================================================
#
# Developers:
#
# Tiago de Conto - tdc.florestal@gmail.com -  https://github.com/tiagodc/
#
# COPYRIGHT: Tiago de Conto, 2020
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
#' @template param-h_step
#' @template param-max-d
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
#' @template reference-olofsson
#' @template reference-thesis
#' @export
stm.hough = function(h_step=0.5, max_d=0.5, h_base = c(1,2.5), pixel_size=0.025, min_density=0.1, min_votes=3){

  if(length(h_base) != 2)
    stop('h_base must be a numeric vector of length 2')

  if(diff(h_base) <= 0)
    stop('h_base[2] must be larger than h_base[1]')

  params = list(
    h_step = h_step,
    max_d = max_d,
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

    survey_points = if(las %>% hasField('Classification')){
      las$Classification != 2
    }else{
      rep(T, las@data %>% nrow)
    }

    if(!hasField(las, 'TreeID')){
      message('no TreeID field found with tree_points signature: performing single stem point classification')
      results = houghStemPoints(las2xyz(las)[survey_points,], h_base[1], h_base[2], h_step, max_d/2, pixel_size, min_density, min_votes)
    }else{
      message('performing point classification on multiple stems')
      survey_points = survey_points & las$TreeID > 0
      results = houghStemPlot(las2xyz(las)[survey_points,], las@data$TreeID[survey_points], h_base[1], h_base[2], h_step, max_d/2, pixel_size, min_density, min_votes)
    }

    las@data$Stem = F
    las@data$Stem[survey_points] = results$Stem

    las@data$Segment = 0
    las@data$Segment[survey_points] = results$Segment

    las@data$Radius = 0
    las@data$Radius[survey_points] = results$Radius

    las@data$Votes = 0
    las@data$Votes[survey_points] = results$Votes

    las = cleanFields(las, c('Radius', 'Votes'))

    return(las)

  }

  func %<>% setAttribute('stem_pts_mtd')

  return(func)

}


#' Stem denoising algorithm: KNN eigen decomposition + point normals intersections voting
#' @description This function is meant to be used inside \code{\link{stemPoints}}. It filters points based on their nearest neighborhood geometries (check \code{\link{fastPointMetrics}}) and assign them to stem patches if reaching a voxel with enough votes.
#' @template param-h_step
#' @template param-max-curvature
#' @template param-max-verticality
#' @param voxel_spacing \code{numeric} - voxel (or pixel) spacing for counting point normals intersections.
#' @template param-max-d
#' @template param-votes-weight
#' @template param-v3d
#' @template section-eigen-decomposition
#' @template section-normals-voting
#' @template reference-liang
#' @template reference-olofsson-2016
#' @export
stm.eigen.knn = function(h_step = .5, max_curvature = .1, max_verticality = 10, voxel_spacing = .025, max_d = .5, votes_weight = .2, v3d = FALSE){

  params = list(
    h_step = h_step,
    max_curvature = max_curvature,
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

  if(max_curvature > 1) stop('max_curvature must be a number between 0 and 1')
  if(max_verticality > 180) stop('max_verticality must be a number between 0 and 180')
  if(votes_weight > 1) stop('votes_weight must be a number between 0 and 1')

  func = function(las){

    las@data$PointID = 1:nrow(las@data)
    mtrlst = c('N', 'Curvature', 'Verticality', 'EigenVector13', 'EigenVector23', 'EigenVector33')

    # if(hasField(las, 'TreeID')){
    #   las@data = las@data %>% split(las@data$TreeID) %>% lapply(LAS) %>%
    #     lapply(fastPointMetrics, method = ptm.knn(), metrics_list=mtrlst) %>%
    #     lapply(function(x) x@data) %>% do.call(what=rbind) %>% as.data.table
    # }

    check_point_metrics = mtrlst %>% sapply(function(x) hasField(las, x)) %>% as.logical %>% all
    if(!check_point_metrics){
      message('Calculating knn fastPointMetrics')
      las = fastPointMetrics(las, ptm.knn(), mtrlst)
    }

    las@data$Stem = with(las@data, N > 3 & Curvature < max_curvature & abs(Verticality - 90) < max_verticality)
    if(hasField(las, 'Classification')){
      las@data$Stem = las@data$Stem & las@data$Classification != 2
    }

    stem_seg = seq(0, max(las$Z)+h_step, h_step)
    las@data$Segment = cut(las$Z, stem_seg, include_lowest=T, right=F, ordered_result=T) %>% as.integer
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
    las@data = las@data[,keepcols,with=F]
    las@data = merge(las@data, votes[,.(PointID, Votes, Radius)], by='PointID', sort=F, all.x=T)
    las@data[!las@data$Stem, c('Votes', 'Radius')] = 0

    by_cols = 'Segment'
    if(hasField(las, 'TreeID')){
      max_votes = las@data[,.(MaxVotes = max(Votes)), by=TreeID]
      las@data = merge(las@data, max_votes, by='TreeID', sort=FALSE, all.x=TRUE)
      las@data$VotesWeight = las@data$Votes / las@data$MaxVotes
      las@data$MaxVotes = NULL
      by_cols %<>% c('TreeID')
    }else{
      las@data$VotesWeight = las@data$Votes / max(las@data$Votes)
    }

    las@data$Stem = las@data$VotesWeight > votes_weight
    radii = las@data[Radius > 0 & Stem, .(med_rad = median(Radius)), by=by_cols]
    las@data = merge(las@data, radii, by=by_cols, sort=FALSE, all.x=TRUE)
    las@data$Radius = las@data$med_rad
    las@data$med_rad = NULL
    las@data[Stem == F]$Radius = 0
    las = cleanFields(las, c('Votes', 'VotesWeight', 'Radius'))

    return(las)
  }

  func %<>% setAttribute('stem_pts_mtd')
  return(func)
}


#' Stem denoising algorithm: Voxel eigen decomposition + point normals intersections voting
#' @description This function is meant to be used inside \code{\link{stemPoints}}. It filters points based on their voxel geometries (check \code{\link{fastPointMetrics}}) and assign them to stem patches if reaching a voxel with enough votes.
#' @template param-h_step
#' @template param-max-curvature
#' @template param-max-verticality
#' @template param-voxel-spacing
#' @template param-max-d
#' @template param-votes-weight
#' @template param-v3d
#' @template section-eigen-decomposition
#' @template section-normals-voting
#' @template reference-liang
#' @template reference-olofsson-2016
#' @export
stm.eigen.voxel = function(h_step = .5, max_curvature = .1, max_verticality = 10, voxel_spacing = .025, max_d = .5, votes_weight = .2, v3d = FALSE){

  params = list(
    h_step = h_step,
    max_curvature = max_curvature,
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

  if(max_curvature > 1) stop('max_curvature must be a number between 0 and 1')
  if(max_verticality > 180) stop('max_verticality must be a number between 0 and 180')
  if(votes_weight > 1) stop('votes_weight must be a number between 0 and 1')

  func = function(las){

    mtrlst = c('N', 'Curvature', 'Verticality', 'EigenVector13', 'EigenVector23', 'EigenVector33', 'VoxelID')

    check_point_metrics = mtrlst %>% sapply(function(x) hasField(las, x)) %>% as.logical %>% all
    if(!check_point_metrics){
      message('Calculating voxel fastPointMetrics')
      las = fastPointMetrics(las, ptm.voxel(voxel_spacing), mtrlst)
    }

    las@data$Stem = with(las@data, Classification != 2 & N > 3 & Curvature < max_curvature & abs(Verticality - 90) < max_verticality)

    stem_seg = seq(0, max(las$Z)+h_step, h_step)
    las@data$Segment = cut(las$Z, stem_seg, include_lowest=T, right=F, ordered_result=T) %>% as.integer
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
    las@data = las@data[,keepcols,with=F]
    las@data = merge(las@data, votes[,.(VoxelID, Votes, Radius)], by='VoxelID', sort=F, all.x=T)
    las@data[!las@data$Stem, c('Votes', 'Radius')] = 0

    by_cols = 'Segment'
    if(hasField(las, 'TreeID')){
      max_votes = las@data[,.(MaxVotes = max(Votes)), by=TreeID]
      las@data = merge(las@data, max_votes, by='TreeID', sort=F, all.x=T)
      las@data$VotesWeight = las@data$Votes / las@data$MaxVotes
      las@data$MaxVotes = NULL
      by_cols %<>% c('TreeID')
    }else{
      las@data$VotesWeight = las@data$Votes / max(las@data$Votes)
    }

    las@data$Stem = las@data$VotesWeight > votes_weight
    radii = las@data[Radius > 0 & Stem, .(med_rad = median(Radius)), by=by_cols]
    las@data = merge(las@data, radii, by=by_cols, sort=FALSE, all.x=TRUE)
    las@data$Radius = las@data$med_rad
    las@data$med_rad = NULL
    las@data[Stem == F]$Radius = 0
    las = cleanFields(las, c('Votes', 'VotesWeight', 'Radius'))

    return(las)
  }

  func %<>% setAttribute('stem_pts_mtd')
  return(func)
}
