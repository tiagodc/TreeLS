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
# require(glue)

library(TreeLS)

### overview of some new methods on v2.0
file = system.file("extdata", "pine.laz", package="TreeLS")
tls = readTLS(file) %>% tlsNormalize()

# calculate some point metrics
tls = fastPointMetrics(tls, ptm.knn())
x = plot(tls, color='Verticality')

# get its stem points
tls = stemPoints(tls, stm.eigen.knn(voxel_spacing = .02))
add_stemPoints(x, tls, size=3, color='red')

# get dbh and height
dbh_algo = shapeFit(shape='cylinder', algorithm = 'bf', n=15, inliers=.95, z_dev=10)
inv = tlsInventory(tls, hp = .95, d_method = dbh_algo)
add_tlsInventory(x, inv)

# segment the stem usind 3D cylinders and getting their directions
seg = stemSegmentation(tls, sgt.irls.cylinder(n=300))
add_stemSegments(x, seg, color='blue')

# check out a specific tree segment
tlsPlot(seg, tls, segment = 3)
