Rcpp::sourceCpp('src/r_interface.cpp')
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

files = dir('../../plantar/results/', full.names = T)
projs = sub('(.+?)_.+', '\\1', files) %>% unique

for(pr in projs[7]){
  paste('.. project', pr) %>% print

  las = pr %>% paste0('_h_stem.rds') %>% readRDS
  vlas = pr %>% paste0('_v.rds') %>% readRDS
  gc()

  estimates = data.frame()
  for(id in unique(las@data[TreeID > 0]$TreeID)){
    paste('.. .. dbh', id) %>% print

    # id = 7
    tree = lasfilter(las, TreeID == id)
    # clear3d() ; rgl.points(tree@data[,.(X,Y,Z)], size=.5)

    dbh = lasfilter(tree, Z > 1.2 & Z < 1.4)
    # plot(dbh)

    # est = estimates[estimates$TreeID == id,]
    est = robustDiameter(dlas=dbh, plot=F, max_d = .25, pixel_size = .02)
    # par(mfrow=c(1,nrow(est)))
    # est[,1:4] %>% apply(1, function(x) tlsPlot.dh(dbh, x, clear = T))
    estimates %<>% rbind(est)
  }

  heights = vlas@data[TreeID > 0, .(h=max(Z)), by=TreeID]
  estimates %<>% merge(heights, by='TreeID')

  out = paste0(pr, '_res.csv')
  write.csv(estimates, out, quote = F, row.names = F)

  par(mfrow=c(1,2), oma=c(0,0,1,0))
  estimates$d %>% hist(main='dbh')
  estimates$h %>% hist(main='h')
  title(main=pr, outer = T)
}



estimates %<>% as.data.table
estimates$d %>% hist
estimates$d %>% mean
estimates$d %>% sd
estimates$d %>% median


estimates[score > 1][order(score)]

estimates[order(TreeID),.N,by=TreeID][N > 1]
plot(las %>% lasfilter(Z > .5), clear_artifacts=F, size=.5) ; pan3d(2)
pos = las@data[,.(mean(X),mean(Y)),by=TreeID]
text3d(pos[,-1] %>% cbind(0), texts = pos$TreeID, color='white')
