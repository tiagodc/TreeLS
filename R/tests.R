# require(TreeLS)
# require(magrittr)
# require(lidR)
# require(rgl)
# las = 'inst/extdata/model_boles.laz'
# las %<>% readTLS
#
# las %<>% tlsSample(val = 0.025, by = 'voxel')
# las %<>% tlsNormalize()
# stc = treeMap(las, min_den = 0.05)
# las %<>% stemPoints_plot(stc)
#
# plot(stc)
# plot(las, color='Stem')
#
#
# for(i in 1:nrow(map)){
#   print(i)
#   tree = map %$% lasclipCircle(las, X[i], Y[i], 1) %>% lasfilter(Classification != 2)
#   tree %<>% stemPoints()
#   stem = lasfilter(tree, Stem)
#   rgl.points(stem@data[,1:3], color='red', size=2)
#   # rgl.points(tree@data[,1:3], color='green', size=.5)
# }
#
# rgl.points(las@data[,1:3], size=.5, color='gray')
# # plot(tree, color='Stem', size=2)
#
# a = las %>% lasfilter(Classification != 2)
# a = las@data[,1:3] %>% as.matrix
#
# b = houghStemPlot(a, map %>% as.matrix, h2 = 2.5) %>% as.data.frame
# las = las@data %>% cbind(b) %>% LAS
#
# plot(las, color='Stem', clear_artifacts=F, size=.5, colorPalette = c('gray','red'))
# # rgl.points(stc@data[,1:3], color='green')
