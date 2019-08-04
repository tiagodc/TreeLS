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
require(nabor)

rm(list = c('.', 'X', 'Y', 'Z', 'Classification', 'TreePosition', 'TreeID', 'Stem', 'Segment', 'gpstime', 'AvgHeight', 'Radius'))

###################

las = readTLS('test_data/ento_u_clip.laz')#, filter='-keep_random_fraction 0.025')
# las = readTLS('inst/extdata/pine.laz')

las %<>% tlsTransform(c('-x','z','y'), T, T)
las %<>% tlsNormalize(keep_ground = F)

map = treeMap(las, map.eigen.knn())

las = treePoints(las, map)

t1 = Sys.time()
las = stemPoints(las, stm.eigen.knn(dvt = .1))
t2 = Sys.time()
print(t2-t1)

th = las@data[las@data$Stem, .(max(Z), .N), by='TreeID']
temp = lasfilter(las, Stem)
plot(temp, color="VotesWeight")

hstep = .5
pln = .2
vrt = 20
vxl = .05
max_d = .5
min_h = 2
min_n = 100

lasfull = las
# las = lasfilter(lasfull, TreeID == 10)

las@data$PointID = 1:nrow(las@data)
las = lasfilter(las, Classification != 2)
# las = pointMetrics(las, ptm.voxels(vxl, F), c('N', 'Planarity', 'Verticality', 'EigenVectors'))

t1 = Sys.time()
laslist = las@data %>% split(las@data$TreeID) %>% lapply(LAS) %>%
  lapply(pointMetrics, method = ptm.knn(), metrics_list=c('N', 'Planarity', 'Verticality', 'EigenVectors')) %>%
  lapply(function(x) x@data) %>% do.call(what=rbind) %>% LAS
t2 = Sys.time()
print(t2-t1)

gc()

t1 = Sys.time()
las = pointMetrics(las, ptm.knn(), c('N', 'Planarity', 'Verticality', 'EigenVectors'))
t2 = Sys.time()
print(t2-t1)

las = lasfilter(las, N > 3 & Planarity < pln & abs(Verticality - 90) < vrt) %>% nnFilter(vxl*2, 10)
# las = lasfilter(las, TreeID == 36)

stemSeg = seq(0, max(las$Z)+hstep, hstep)
las@data$Segment = cut(las$Z, stemSeg, include_lowest=T, right=F, ordered_result=T) %>% as.integer
las@data$Segment[las@data$Z < 0] = 0

# voxels = las@data[order(Segment, VoxelID), .(X=mean(X), Y=mean(Y), Z=mean(Z), e1=mean(EigenVector13), e2=mean(EigenVector23), e3=mean(EigenVector33)), by=.(Segment, VoxelID)]
# dups = duplicated(voxels$VoxelID)
# voxels = voxels[!dups]

voxels = las@data[order(TreeID, Segment, PointID), .(TreeID, Segment, PointID, X, Y, Z, EigenVector13, EigenVector23, EigenVector33)]

sgs = voxels$Segment
ids = voxels$PointID
tds = voxels$TreeID
a = voxels[,-c(1:3)] %>% as.matrix

# a = split(voxels, trids) %>% lapply(function(x){
#   as.matrix(x[,-c(1:3)]) %>% treeEigenHough(x$Segment, vxl, max_d, F)
# })

b = plotEigenHough(a, ids, tds, sgs, vxl, max_d/2, T, F)
# b = treeEigenHough(a, ids, sgs, vxl, max_d/2, F, F)
# b = b[[1]]

sids = ids %>% unique %>% sort
segs = b %>% lapply(function(x) x %>% do.call(what=cbind)) %>% do.call(what=rbind) %>% as.data.table
colnames(segs) = c('Votes','Radius','PointID', 'Segment', 'TreeID')

# segs = lapply(b, function(x) rev(x)[-1] %>% do.call(what=rbind) %>% cbind(rev(x)[1])) %>% do.call(what=rbind) %>% as.data.table()

voxels = merge(voxels, segs, by='PointID', sort=F)
las@data = merge(las@data, voxels[,.(PointID, Votes, Radius)], by='PointID', sort=F)

plot(las, size=2, color='Votes', clear_artifacts=F)
rgl.points(tlsSample(las2, randomize(.2))@data, size=.5)
pan3d(2)
voxels$Votes %>% hist

temp = lasfilter(las, Votes > 100)
plot(temp, size=2, color='Votes', clear_artifacts=F)
rgl.points(tlsSample(las2, randomize(.2))@data, size=.5)
pan3d(2)

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

