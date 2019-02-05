require(magrittr)
require(lidR)
require(rgl)

las = 'd:/Projects/TLStools/sample_data/square.las'

las %<>% readLAS(select='XYZ', filter='-keep_z 2 3 -thin_with_voxel 0.01')

plot(las)

m = singleStack(las@data %>% as.matrix, 0.025) %>% do.call(what=cbind) %>% as.data.frame %>% LAS

clear3d()
rgl.points(las@data, size=.5)
rgl.points(m@data, col='red')

dim(m@data)
m %<>% lasfilterduplicates()
dim(m@data)
