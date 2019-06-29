Rcpp::sourceCpp('src/r_interface.cpp')
source('R/methods.R')
source('R/sample_methods.R')
source('R/map_methods.R')
source('R/stem_points_methods.R')
source('R/stem_segmentation_methods.R')
source('R/tree_points_methods.R')

require(magrittr)
require(lidR)
require(rgl)
require(data.table)
require(dismo)
require(deldir)
require(RANN)

rm(list = c('.', 'X', 'Y', 'Z', 'Classification', 'TreePosition', 'TreeID', 'Stem', 'Segment', 'gpstime', 'AvgHeight', 'Radius'))

###################

# nn = RANN::nn2(da ta = tls %>% las2xyz, k=20, treetype = 'kd', searchtype = 'standard')

# tls = readTLS('test_data/zeb.laz', filter='-keep_circle -66 249 10')
#
# tls %<>% tlsNormalize(keepGround = F)
#
# thin = tlsSample(tls, voxelize(.025))
# map = treeMap(thin, map.hough(min_votes=5, min_density = .05))
#
# tls %<>% treePoints(map)
#
# tls %<>% stemPoints()
#
# df = stemSegmentation(tls)
#
# tlsPlot(tls, df, treeID = 47)
# tree = lasfilter(tls, TreeID == 47)
#
# plot(tree)

tree = readTLS('inst/extdata/pine.laz')

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
