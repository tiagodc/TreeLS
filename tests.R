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
require(benchmarkme)
require(glue)

rm(list = c('.', 'X', 'Y', 'Z', 'Classification', 'TreePosition', 'TreeID', 'Stem', 'Segment', 'gpstime', 'AvgHeight', 'Radius'))

###################
las = readTLS('test_data/Parcela.las', filter='-keep_random_fraction 0.25', select='xyzi')
las = readTLS('inst/extdata/pine.laz', select='xyzi')
plot(las)
las = tlsNormalize(las, keep_ground = F)
thin = tlsSample(las, smp.voxelize(.025))
map = treeMap(thin, map.hough())
map %>%
  # treeMap.merge(.1) %>%
  treeMap.positions()
las = treePoints(las, map, trp.crop(1))
plot(las, color='TreeID')
las = stemPoints(las, stm.hough())
plot(las,color='Stem')
segs = stemSegmentation(las, sgt.ransac.cylinder(n = 10, inliers = .98))
inv = tlsInventory(las)

seg = filter_poi(las, Segment == 3 & TreeID == 3 & Stem)
plot(seg, color='Votes', clear_artifacts=F, size=2)
temp = seg %>% filter_poi(Votes == max(Votes))
seg@data$Radius %>% median


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

