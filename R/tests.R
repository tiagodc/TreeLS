# require(magrittr)
# require(lidR)
# require(rgl)
# require(TreeLS)
# las = 'd:/Projects/TLStools/test_clouds/gerdau.laz'
# las %<>% readTLS
# las %<>% tlsSample(val = 0.025, by = 'voxel')
# las %<>% tlsNormalize()
# stc = treeMap(las, pixel = 0.025)
# map = treePositions(stc)
#
# i = 23
# tree = map %$% lasclipCircle(las, X[i], Y[i], 1) %>% lasfilter(Classification != 2)
#
# tree %<>% stemPoints()
# plot(tree, color='Votes', size=2)
#
# tree@data$Radius
#
# a = lasfilter(tree, Stem)
# plot(a)
#
# plot(las), color='Classification')
# axes3d(col='white')
#
# plot(stc)
#
# stc@data %>% head
#
# rgl.points(stc@data, color='red')
#
# a = stc@data %>% as.matrix

