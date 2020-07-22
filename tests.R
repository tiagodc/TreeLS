Rcpp::sourceCpp('src/r_interface.cpp', rebuild = T)
source('R/methods.R')
source('R/sample_methods.R')
source('R/map_methods.R')
source('R/stem_points_methods.R')
source('R/stem_segmentation_methods.R')
source('R/tree_points_methods.R')
source('R/point_metrics_methods.R')
source('R/plot_extras.R')

require(magrittr)
require(lidR)
require(rgl)
require(data.table)
require(dismo)
require(deldir)
require(nabor)

rm(list = c('.', 'X', 'Y', 'Z', 'Classification', 'TreePosition', 'TreeID', 'Stem', 'Segment', 'gpstime', 'AvgHeight', 'Radius'))

###################

las = readTLS('inst/extdata/pine.laz')
las = tlsNormalize(las, keep_ground = F)
# las = pointMetrics(las, ptm.knn(10))
las = stemPoints(las, stm.eigen.voxel(voxel_spacing = .05))
plot(las,color='Stem')
rad = las@data[Stem == T, mean(Radius), by='Segment']
segs = stemSegmentation(las, sgt.ransac.circle(.1))
segs

seg = filter_poi(las, Segment == 3)
plot(seg, color='Radius', clear_artifacts=F, size=2)
temp = seg %>% filter_poi(Votes == max(Votes))
seg@data$Radius %>% median


cols = seg@data$Votes %>% max %>% height.colors
# plot(temp, clear_artifacts=F)
for(i in 1:nrow(temp@data)){
  row = temp@data[i,]
  spheres3d(row$X, row$Y, row$Z, .005, col='white')
  lines3d(
    c(row$X - row$Radius * row$EigenVector13, row$X + row$Radius * row$EigenVector13),
    c(row$Y - row$Radius * row$EigenVector23, row$Y + row$Radius * row$EigenVector23),
    c(row$Z - row$Radius * row$EigenVector33, row$Z + row$Radius * row$EigenVector33),
    color=cols[row$Votes], lwd=1
  )
}

rgl.points(seg@data[,.(X,Y,Z)], size=2, color='red')

points = las@data[Stem & order(Segment, PointID), .(Segment, PointID, X, Y, Z, EigenVector13, EigenVector23, EigenVector33)]
votes = points[,-c(1:2)] %>% as.matrix %>% treeEigenHough(points$PointID, points$Segment, 0.025, .25, F, T) %>%
  lapply(function(x) x %>% do.call(what=rbind))# %>% do.call(what=rbind) %>% as.data.table

vt = toLAS(votes[[5]], c('Votes', 'X', 'Y', 'Z'))
seg = tail(vt@data, 1)[[1]]
seg = filter_poi(las, Segment == seg)
vt@data = vt@data[-nrow(vt@data),][X > -10]
vt@data %>% apply(2,range)
plot(vt, color='Votes', clear_artifacts=F)

rgl.points(seg@data[,.(X,Y,Z)], color='white', size=1)

