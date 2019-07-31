Rcpp::sourceCpp('src/r_interface.cpp')
source('R/methods.R')
source('R/sample_methods.R')
source('R/map_methods.R')
source('R/stem_points_methods.R')
source('R/stem_segmentation_methods.R')
source('R/tree_points_methods.R')
source('R/point_metrics_methods.R')

require(magrittr)
require(lidR)
require(rgl)
require(data.table)
require(dismo)
require(deldir)
require(RANN)

rm(list = c('.', 'X', 'Y', 'Z', 'Classification', 'TreePosition', 'TreeID', 'Stem', 'Segment', 'gpstime', 'AvgHeight', 'Radius'))

###################

las = readTLS('test_data/ento_u_clip.laz')#, filter='-keep_random_fraction 0.025')
las = readTLS('inst/extdata/pine.laz')

las %<>% tlsTransform(c('x','z','-y'), T, T)
las %<>% tlsNormalize(keep_ground = F)

map = treeMap(las, map.eigen.knn())
# plot(map)

# map$TreeID %>% unique
# tree = map %>% lasfilter(TreeID == 5)

tree=map

voxel_size = .05
max_radius = .25

tree %<>% pointMetrics(ptm.knn())
toco = lasfilter(tree, Z > 5 & Z < 7)
plot(toco)

# a = toco@data[N > 3,.(X = mean(X), Y=mean(Y), Z=mean(Z), e1=mean(EigenVector13), e2=mean(EigenVector23), e3=mean(EigenVector33)), by=VoxelID]
# a = a[,-1] %>% as.matrix
a = toco@data[,.(X,Y,Z,EigenVector13,EigenVector23,EigenVector33)] %>% as.matrix

b = hough3d(a, voxel_size) %>% do.call(what = rbind)
px = cbind(b[,-1]) %>% toLAS()
px@data$hough = b[,1]
px = lasfilter(px, hough > 100 & X < 1e4 & Y < 1e4)

px@data %>% apply(2,range)

plot(px, color='hough', clear_artifacts=F, size=3)
rgl.points(toco@data)

