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

cyl = tlsCylinder(1000, .5, .12, .04)
cyl@data$X = cyl@data$X + 37
cyl@data$Y = cyl@data$Y - 57
cyl@data$Z = cyl@data$Z + 12

cyl@data[,1:3] = as.data.table(as.matrix(cyl@data[,1:3]) %*% rotationMatrix(20*pi/180, 12*pi/180, 0))

est = cylinderFit(cyl, 'bf', n=15) %T>% print
tlsPlot.dh(cyl, est, clear = T)

