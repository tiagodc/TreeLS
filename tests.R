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

n = 30
d = 0.2
search = c('knn', 'sphere', 'voxel')
search = search[1]
stype = ifelse(search == 'knn', 'standard', ifelse(search == 'sphere', 'radius', search))

knn = RANN::nn2(tree %>% las2xyz(), k = n+1, treetype = 'kd', searchtype = stype, radius = d)

kid = knn$nn.idx[,-1]
kds = knn$nn.dists[,-1]

a = temp(tree %>% las2xyz, kid)

b = las2xyz(tree)[kid[2,], ]
eigen(b %>% cov)
a[[2]]

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
