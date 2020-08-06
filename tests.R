Rcpp::sourceCpp('src/r_interface.cpp', rebuild = T)
source('R/sample_methods.R')
source('R/map_methods.R')
source('R/stem_points_methods.R')
source('R/stem_segmentation_methods.R')
source('R/tree_points_methods.R')
source('R/point_metrics_methods.R')
source('R/plot_extras.R')
source('R/methods.R')

require(magrittr)
require(lidR)
require(rgl)
require(data.table)
require(dismo)
require(deldir)
require(nabor)
require(benchmarkme)
require(glue)

rm(list = c('.', 'X', 'Y', 'Z', 'Classification', 'TreePosition', 'TreeID', 'Stem', 'Segment', 'gpstime', 'AvgHeight', 'Radius'))

###################
require(TreeLS)

tls = readTLS('inst/extdata/pine.laz') %>%
  tlsNormalize() %>%
  stemPoints(stm.eigen.knn())

segs = stemSegmentation(tls, sgt.bf.cylinder())

chunk = filter_poi(tls, Segment == 5 & Stem)
shapeFit(chunk, 'circle', 'irls')
