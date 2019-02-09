# require(magrittr)
# require(lidR)
# require(rgl)
# require(TreeLS)
#
# las = 'd:/Projects/TLStools/test_clouds/gerdau.laz'
#
# las %<>% readTLS
# #
# # dim(las@data)
# las %<>% tlsSample(val = 0.1, by = 'voxel')
# # dim(las@data)
# #
# # las@data %>% nrow
# # las %<>% lasfilter(a)
# # las@data %>% nrow
# #
# las %<>% tlsNormalize()
# #
# # plot(las), color='Classification')
# # axes3d(col='white')
# #
# stc = treeMap(las)
# # plot(stc)
# #
# # stc@data %>% head
# #
# # rgl.points(stc@data, color='red')
# #
# # a = stc@data %>% as.matrix
#
# las = system.file("extdata", "pine.laz", package="TreeLS")
