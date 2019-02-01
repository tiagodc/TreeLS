require(magrittr)
require(lidR)

las = 'd:/Projects/TLStools/test_clouds/gerdau.laz'

las %<>% readLAS(select='XYZ', filter='-keep_z 1 2')

plot(las)

m = singleStack(las@data %>% as.matrix, 0.025) %>% do.call(what=cbind) %>% as.data.frame %>% LAS

plot(m)
dim(m@data)
m %<>% lasfilterduplicates()
dim(m@data)

m %<>% lasfilterduplicates
