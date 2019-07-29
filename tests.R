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

las %<>% lasfilterduplicates()

las %<>% tlsTransform(c('-x','z','y'), T, T)
las %<>% tlsNormalize(keep_ground = F)


map = treeMap(las, map.eigen.voxel(.2,20,.07))
# plot(map)

tree = lasfilter(map, TreeID == 18)
# plot(tree)

voxel_size = .07
max_radius = .25

minxyz = tree@data[,1:3] %>% apply(2,min) %>% as.double

vxid = tree@data[,.N,by=VoxelID]
vxid = vxid$VoxelID[vxid$N %>% which.max]

vx = lasfilter(tree, VoxelID == vxid)
normal = eigen(vx@data[,1:3] %>% cov)$vectors[,3] %>% as.double
center = vx@data[,1:3] %>% colMeans %>% as.double

dists = seq(-max_radius, max_radius+voxel_size, by=voxel_size)
xs = center[1] + dists*normal[1]
ys = center[2] + dists*normal[2]
zs = center[3] + dists*normal[3]

xv = floor( (xs - minxyz[1]) / voxel_size )
yv = floor( (ys - minxyz[2]) / voxel_size )
zv = floor( (zs - minxyz[3]) / voxel_size )

xyzv = cbind(xv,yv,zv)
vx = floor( (center - minxyz) / voxel_size )

plot(vx, clear_artifacts=F, size=3)
lines3d(rbind(minpt,maxpt),color='white',lwd=2)

tree %<>% pointMetrics(ptm.voxels(.07))
a = tree@data[N > 3,.(X = mean(X), Y=mean(Y), Z=mean(Z), e1=mean(EigenVector13), e2=mean(EigenVector23), e3=mean(EigenVector33)), by=VoxelID]

