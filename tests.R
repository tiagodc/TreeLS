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

las %<>% tlsTransform(c('-x','z','y'), T, T)
las %<>% tlsNormalize(keep_ground = F)

map = treeMap(las, map.eigen.voxel(vxl = .05))

las = treePoints(las, map)

hstep = .5
pln = .2
vrt = 20
vxl = .05
max_d = .5
min_h = 2
min_n = 100

lasfull = las

las = lasfilter(las, Classification != 2)
las = pointMetrics(las, ptm.voxels(vxl, F), c('N', 'Planarity', 'Verticality', 'EigenVectors'))

las = lasfilter(las, N > 3 & Planarity < pln & abs(Verticality - 90) < vrt) %>% nnFilter(vxl*2, 10)
# las = lasfilter(las, TreeID == 36)

stemSeg = seq(0, max(las$Z)+hstep, hstep)
las@data$Segment = cut(las$Z, stemSeg, include_lowest=T, right=F, ordered_result=T) %>% as.integer
las@data$Segment[las@data$Z < 0] = 0

voxels = las@data[order(Segment, VoxelID), .(X=mean(X), Y=mean(Y), Z=mean(Z), e1=mean(EigenVector13), e2=mean(EigenVector23), e3=mean(EigenVector33)), by=.(Segment, VoxelID)]
dups = duplicated(voxels$VoxelID)
voxels = voxels[!dups]

ids = voxels$Segment
a = voxels[,-c(1:2)] %>% as.matrix

b = treeEigenHough2d(a, ids, vxl/2, max_d/2)

sids = ids %>% unique %>% sort
segs = 1:length(b) %>% lapply(function(x){
  data.table(VoxelID = voxels[Segment == sids[[x]]]$VoxelID, Votes = b[[x]][[1]], Radius = b[[x]][[2]])
}) %>% do.call(what=rbind)

voxels = merge(voxels, segs, by='VoxelID', sort=F)
las@data = merge(las@data, voxels[,.(VoxelID, Votes, Radius)], by='VoxelID', sort=F)


plot(las, size=2, color='Votes')
voxels$Votes %>% hist

temp = lasfilter(las, Votes > 20)
plot(temp, size=3, color='Radius')

b[[1]][[1]] %>% length

rgl.points(voxels[,4:6], size=.5)

las@data$TreeID = 0
maxdst = max_d*2
i = 1
while(any(las@data$TreeID == 0)){
  xy = las@data[TreeID == 0,.(X,Y)][1,] %>% as.double
  dst = las@data[,sqrt( (X-xy[1])^2 + (Y-xy[2])^2 )] %>% as.double
  las@data[TreeID == 0 & dst < maxdst]$TreeID = i
  i=i+1
}

hn = las@data[,.(H=max(Z) - min(Z), .N), by=TreeID]
las = lasfilter(las, !(TreeID %in% hn$TreeID[hn$H < min_h] | hn$N[TreeID] < min_n) )

