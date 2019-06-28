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

# nn = RANN::nn2(da ta = tls %>% las2xyz, k=20, treetype = 'kd', searchtype = 'standard')

las = readTLS('test_data/zeb.laz', filter='-keep_circle -66 249 10') #%>% lasfilter(Z > 2 & Z < 4)

# xy = las@data[,1:2] %>% apply(2,mean) %>% as.double
# las %<>% tlsCrop(xy[1], xy[2], 10)


