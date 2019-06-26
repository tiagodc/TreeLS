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

las = readTLS('test_data/cema01_52_16_u45_norm_plot.laz') #%>% lasfilter(Z > 2 & Z < 4)

max_size = 1E6
d = 0.1
n = 3

npts = nrow(las@data)
zclass = 0

if(npts > max_size){
  nparts = ceiling(npts/max_size)
  probs = seq(0, 1, by = 1/nparts)
  probs = probs[probs > 0 & probs < 1]
  zqts = quantile(las$Z, probs) %>% as.double
  hs = c(min(las$Z)-1, zqts, max(las$Z)+1)
  zclass = cut(las$Z, hs) %>% as.integer
}

keep = rep(T, nrow(las@data))

for(i in unique(zclass)){
  bool = zclass == i
  xyz = las@data[bool,.(X,Y,Z)]
  rnn = RANN::nn2(data = xyz, radius = d, treetype = 'kd', searchtype = 'radius')$nn.idx %>% as.data.frame

  knn = rep(0, nrow(rnn))
  for(j in 2:ncol(rnn)){
    temp = ifelse(rnn[,j] > 0, 1, 0) %>% as.double
    knn = knn + temp
  }

  keep[bool] = knn > n
}

las = lasfilter(las, keep)
