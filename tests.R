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

# las = readTLS('test_data/zeb.laz', filter='-keep_random_fraction 0.025')
las = readTLS('inst/extdata/pine.laz')

las %<>% las2xyz %>% toLAS

vx = ptm.voxels(las, d = 0.15)
idx = split(1:length(vx), vx)

temp = voxelMetrics(las2xyz(las), idx) %>% do.call(what = rbind) %>% as.data.table
colnames(temp) = c('planarity', 'verticality', 'linearSaliency', 'planarSaliency', 'scattering', 'anisotropy', 'zrange', 'zsd')
temp$VoxelID = names(idx) %>% as.double

las@data$VoxelID = vx #%>% as.factor %>% as.double
las@data = merge(las@data, temp, by='VoxelID')

cpa = lidR::pastel.colors(las@data$VoxelID %>% unique %>% length)
plot(las, color='VoxelID', colorPalette = cpa, size=2)



las@data$trunk = with(las@data, {planarity < .05})
plot(las, color='trunk')

vxl = lasfilter(las, voxel == sample(unique(vx), 1) )
apply(vxl %>% las2xyz, 2, function(x) diff(range(x)))
plot(vxl) ; axes3d(col='white')


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
