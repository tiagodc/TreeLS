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
# las = readTLS('inst/extdata/pine.laz')

las %<>% tlsAlter(c('-x','z','y'), T, T)
xy = las@data[,.(mean(X),mean(Y))] %>% as.double
vx = tlsCrop(las, xy[1], xy[2], 5)
vx %<>% tlsNormalize(keepGround = F)
vx %<>% ptm.voxels(d = 0.1)
vx %<>% ptm.knn(30)
# vx %<>% ptm.radius(metrics_list = 'MeanDistance')

pt = vx@data$VoxelID %>% unique %>% length %>% pastel.colors

plot(vx, color='Intensity')#, colorPalette=pt)


vx@data$MedianDistance %>% hist
trunk = with(vx@data, N > 3 & MeanDistance < .05 & Planarity < .2 & Verticality > 80 & Verticality < 100)
vx %<>% lasadddata(trunk,'trunk')

temp = lasfilter(vx, N > 3 & MeanDistance < .05 & Planarity < .2 & Verticality > 80 & Verticality < 100)

plot(vx, size=.5, color='trunk')
plot(temp)


cb = lasfilter(vx, VoxelID == 19)
cb %>% las2xyz %>% temp_func
cb %>% las2xyz %>% cov %>% eigen()
cb@data[1,]

n = 20
d = 0.2
search = c('knn', 'sphere', 'voxel')
search = search[1]
stype = ifelse(search == 'knn', 'standard', ifelse(search == 'sphere', 'radius', search))

knn = RANN::nn2(tree %>% las2xyz(), k = n+1, treetype = 'kd', searchtype = stype, radius = d)

kid = knn$nn.idx[,-1]
kds = knn$nn.dists[,-1]

ptm = pointMetrics(tree %>% las2xyz, kid) %>% do.call(what = rbind) %>% as.data.table
colnames(ptm) = c('planarity', 'verticality', 'linearSaliency', 'planarSaliency', 'scattering', 'anisotropy', 'zrange', 'zsd')

ptm %>% apply(2,range)

bole = lasfilter(tree, ptm$planarSaliency > .5 & ptm$verticality > 80 & ptm$verticality < 100)
plot(bole)

for(i in 1:ncol(ptm)){
  tree %<>% lasadddata(ptm[,i], 'temp')
  plot(tree, color='temp')
}

# temp = lasfilter(tree, Stem)
# ps = seq(min(tree$Z)-1, max(tree$Z)+1, by = .5)
# bool = cut(tree$Z, ps) %>% as.integer
#
# tlist = lapply(bool %>% unique, function(x) lasfilter(tree, bool == x))
#
# for(i in tlist){
#   plot(i, clear_artifacts=F, color='Stem', size=2)
#   axes3d(col='white')
# }
