Rcpp::sourceCpp('src/r_interface.cpp', rebuild = T)
source('R/sample_methods.R')
source('R/map_methods.R')
source('R/stem_points_methods.R')
source('R/stem_segmentation_methods.R')
source('R/tree_points_methods.R')
source('R/point_metrics_methods.R')
source('R/plot_extras.R')
source('R/methods.R')

require(magrittr)
require(lidR)
require(rgl)
require(data.table)
require(dismo)
require(deldir)
require(nabor)
require(benchmarkme)
require(glue)

rm(list = c('.', 'X', 'Y', 'Z', 'Classification', 'TreePosition', 'TreeID', 'Stem', 'Segment', 'gpstime', 'AvgHeight', 'Radius'))

###################
las = readLAS('test_data/Parcela.las', filter='-keep_random_fraction 0.25', select='xyzi')
las = readLAS('inst/extdata/pine.laz', select='xyzi')
nrow(las@data)

las = fastPointMetrics(las, ptm.voxel())
las@data$KnnDensity %>% hist

plot(las, color='VoxelID', colorPalette=pastel.colors(10000))
las = tlsNormalize(las, keep_ground = F)
# thin = tlsSample(las, smp.voxelize(.025))
# map = treeMap(thin, map.hough())
# las = treePoints(las, map, trp.crop(1))
las = stemPoints(las, stm.eigen.knn())
plot(las,color='ZRange')
segs = stemSegmentation(las, sgt.bf.cylinder(n = 10, inliers = .95))
inv = tlsInventory(las)

tlsPlot(las, inv, fast=F)

vm = voxelMetrics(las2xyz(las), las@data$VoxelID, POINT_METRICS_NAMES %in% 'N') %>% do.call(what = rbind) %>% unlist %>% as.double

cols = seg@data$Votes %>% max %>% height.colors
# plot(temp, clear_artifacts=F)
for(i in 1:nrow(seg@data)){
  row = seg@data[i,]
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

cloud = 1:12 %>% matrix(ncol=3)
dists = c()

Rcpp::sourceCpp('test.cpp', rebuild = T)
pointDistances()

for(i in 1:3){
  for(j in (i+1):4){
    sumsq = 0;
    for(k in 1:3) sumsq = sumsq + (cloud[j,k] - cloud[i,k])^2;
    dists = c(dists, sqrt(sumsq))
  }
}
