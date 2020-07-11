Rcpp::sourceCpp('src/r_interface.cpp', rebuild = T)
source('R/methods.R')
source('R/sample_methods.R')
source('R/map_methods.R')
source('R/stem_points_methods.R')
source('R/stem_segmentation_methods.R')
source('R/tree_points_methods.R')
source('R/point_metrics_methods.R')
source('R/plot_extras.R')

require(magrittr)
require(lidR)
require(rgl)
require(data.table)
require(dismo)
require(deldir)
require(nabor)

rm(list = c('.', 'X', 'Y', 'Z', 'Classification', 'TreePosition', 'TreeID', 'Stem', 'Segment', 'gpstime', 'AvgHeight', 'Radius'))

###################

las = readTLS('./test_data/Parcela.las', filter='-keep_random_fraction 0.5')

las = tlsNormalize(las, keep_ground = F)

map = treeMap(las, map.eigen.knn(.05, 5, .05), merge = .1)

las = treePoints(las, map, trp.crop(1.5, circle = F))
plot(las, color='TreeID', colorPalette=pastel.colors(100))

las@data[,23:39] = NULL
las@data %>% names
las = stemPoints(las, stm.eigen.knn())
las %>% tlsSample(smp.randomize(.2)) %>% plot(color='VotesWeight')

segs = stemSegmentation(las, sgt.ransac.circle())
