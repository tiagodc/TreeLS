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

files = dir('/media/tiago/DATA/Bracell', full.names = T)

f = 16
# for(f in 5:length(files)){
  gc() ; dev.off()
  las = readTLS(files[f])

  # plot(las)

  las = tlsNormalize(las, .2, F)
  thin = tlsSample(las, smp.voxelize(.02))

  map = treeMap(thin, map.hough(max_d=.3, min_h = 2.5, max_h = 4.5, min_density = .15))

  treePositions(map)
  map = treeMapAggregate(map)
  tp = treePositions(map)

  plot(map, clear_artifacts=F)
  rgl.points(lasfilter(thin, Z < 6)@data[,.(X,Y,Z)], size=.5)

  las = treePoints(las, map, trp.crop(1,F))
  las = stemPoints(las, stm.hough(.5, h_base = c(1,3), max_radius = .15, pixel_size = .025, min_density = .2))

  h = las@data[TreeID > 0, .(H = max(Z)), by='TreeID']
  las_dbh = lasfilter(las, Z > 1.2 & Z < 1.4 & Stem)

  par(mfrow=c(2,2))
  dbh = data.table()
  for(d in unique(las_dbh$TreeID)){
    temp = lasfilter(las_dbh, TreeID == d)
    est = circleFit(temp, 'ransac', 15, n_best = 30)
    tlsPlot.dh(temp, est, F)
    est$TreeID = d
    # est = robustDiameter(temp)
    dbh = rbind(dbh,est)
  }

  rgl.points(lasfilter(thin, Z < 5 & Z > .5)@data[,.(X,Y,Z)], size=1)
  text3d(tp$X, tp$Y, 0, tp$TreeID, size=1.5, col='yellow')

  inv = merge(dbh, h, by='TreeID')
  inv = inv[!(TreeID %in% c(202))]

  nm = sub('_clip\\.la[sz]', '_inv.csv', files[f])
  nm = sub('.+/(.+\\.csv)', '../bracell_results/\\1', nm)

  hist(inv$d, main=nm, xlab=mean(inv$d) %>% round(2))
  write.table(inv, nm, quote = F, row.names = F, sep=',')
# }


v = c()
for( i in 58:73){
  tls = read.csv(paste0('../bracell_results/madura_11_',i,'_inv.csv')) %>% as.data.table
  validation = read.csv('../bracell_results/011_madura.csv', dec=',') %>% as.data.table
  val = validation[nm_parcela==i & nm_dap != 0]
  v = c(v, mean(val$nm_dap) - mean(tls$d))
}

mean(v)
