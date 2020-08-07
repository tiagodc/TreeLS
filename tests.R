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

file = 'test_data/Nelder Faz_Jugui 13-09-19.las'
hd = readLASheader(file)

xs = seq(hd@PHB$`Min X`, hd@PHB$`Max X`, length.out = 4)
ys = seq(hd@PHB$`Min Y`, hd@PHB$`Max Y`, length.out = 4)

j = i = 1
for(i in 1:3){
  for(j in 1:3){
    filt = glue('-keep_x {xs[i]} {xs[i+1]} -keep_y {ys[j]} {ys[j+1]}')
  }
}

tls = readTLS(file, filter=filt, select='xyz')
tls = tlsNormalize(tls, min_res = .5, keep_ground = F)
thin = tlsSample(tls, smp.voxelize(0.025))
map = treeMap(thin, map.hough(min_density = .2, max_d = .4))
tls = treePoints(tls, map, trp.crop(.75))
tls = stemPoints(tls, stm.hough(max_d=.4))
inv = tlsInventory(tls)

tiny = tlsSample(tls, smp.randomize(.1))
tlsPlot(tiny, inv, tree_id = 31)

inv$H %>% hist
hist(inv$Radius*200)

# .rs.restartR()
