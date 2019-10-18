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

files = dir('~/Desktop/bracell/', pattern = '\\.laz$', full.names = T)

for(f in files){
  print(f)

  las = readLAS(f)
  las = tlsTransform(las, c('z','x','y'), rotate = T)

  if(grepl('_v\\.laz$',f)) las = tlsRotate(las)

  nm = sub('\\.laz$', '_fix.laz', f)
  writeLAS(las, nm)
}

