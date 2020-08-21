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

rm(list = c('X','Y','Z','Classification','TreePosition','TreeID','Stem','Segment','gpstime','AvgHeight','Radius','N','Curvature','Verticality','MeanDist','PX','PY','PZ','h_radius','PointID','VoxelID','StemID','EigenVector13','EigenVector23','EigenVector33','Votes','absRatio','clt','r','v','x','y'))

###################
# require(TreeLS)
