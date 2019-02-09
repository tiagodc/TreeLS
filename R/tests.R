# require(magrittr)
# require(lidR)
# require(rgl)
require(TreeLS)
las = 'd:/Projects/TLStools/test_clouds/gerdau.laz'
las %<>% readTLS
las %<>% tlsSample(val = 0.025, by = 'voxel')
las %<>% tlsNormalize()
stc = treeMap(las, pixel = 0.025)
map = treePositions(stc)

i = 12
tree = map %$% lasclipCircle(las, X[i], Y[i], 1) %>% lasfilter(Classification != 2)

tree %<>% stemPoints()
plot(tree, color='Stem', size=2)

a = lasfilter(tree, Stem)
plot(a)

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
#
require(TreeLS)
las = system.file("extdata", "spruce.laz", package="TreeLS") %>% readTLS

las %<>% stemPoints()

plot(las, color='Stem')



