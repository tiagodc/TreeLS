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

# for(pr in projs[7]){
  pr = projs[6]
  paste('.. project', pr) %>% print

  las = pr %>% paste0('_h_stem.rds') %>% readRDS
  # vlas = pr %>% paste0('_v.rds') %>% readRDS
  res = pr %>% paste0('_res.csv') %>% read.csv %>% as.data.table

  dev.off()
  par(mfrow=c(1,2), oma=c(0,0,1,0))
  res$d %>% hist(main='d')
  res$h %>% hist(main='h')
  title(main=pr, outer = T)

  res$d %>% mean
  res$d %>% median
  res$d %>% sd

  las %<>% tlsSample(randomize(.33)) %>% lasfilter(Z > .5 & Z < 3)
  plot(las, clear_artifacts=F, size=.5, colorPalette='white') ; pan3d(2)
  txt = las@data[TreeID > 0, .(X=mean(X), Y=mean(Y), Z=min(Z)-.5),by=TreeID]
  text3d(txt[,-1], col='yellow', cex=1.5, texts = txt$TreeID)

  par(mfrow=c(1,1))
  for(i in 1:nrow(res)){
    pars = res[i,]
    dbh = lasfilter(las, TreeID == pars$TreeID & Z > 1.2 & Z < 1.4)
    tlsPlot.dh(dbh, pars[,2:5], F)
  }

  redone = data.frame()
  # check = c(5,34,37,38,171,93,479,107,102,109,69,241,87)
  # res[score > 1]
  # res[TreeID %in% check]

  # res = res[StemID == 1]
  # pr %>% paste0('_res_ok.csv') %>% write.csv(x = res)

  # estimates = data.frame()
  # for(id in unique(las@data[TreeID > 0]$TreeID)){
  # paste('.. .. dbh', id) %>% print

    id = 113
    tree = lasfilter(las, TreeID == id)
  #   # clear3d() ; rgl.points(tree@data[,.(X,Y,Z)], size=.5)
  #
    dbh = lasfilter(tree, Z > 1.2 & Z < 1.4)
    # dbh %<>% gpsTimeFilter(to=.4)
    # plot(dbh)
  #
    est = res[TreeID == id,2:5]
    par(mfrow=c(1,nrow(est))) ; clear3d()
    est[,1:4] %>% apply(1, function(x) tlsPlot.dh(dbh, x, clear = F))

    est = robustDiameter(dlas=dbh, plot=T, max_d = .12, pixel_size = .01, min_den = .33) ; est %>% nrow
    # est = est[1,]
    par(mfrow=c(1,nrow(est))) ; clear3d()
    plot(tree, clear_artifacts=F, size=.5, colorPalette='white')
    est[,1:4] %>% apply(1, function(x) tlsPlot.dh(dbh, x, clear = F))
    # estimates %<>% rbind(est)
  # }

    redone %<>% rbind(est)

    redone %<>% merge(res[,.(TreeID,h)], by='TreeID')

    res = res[!(TreeID %in% redone$TreeID)]
    res %<>% rbind(redone)
    pr %>% paste0('_res_ok.csv') %>% write.csv(x = res)

    res$d %>% hist

  # heights = vlas@data[TreeID > 0, .(h=max(Z)), by=TreeID]
  # estimates %<>% merge(heights, by='TreeID')
  #
  # out = paste0(pr, '_res.csv')
  # write.csv(estimates, out, quote = F, row.names = F)
  #
  # par(mfrow=c(1,2), oma=c(0,0,1,0))
  # estimates$d %>% hist(main='dbh')
  # estimates$h %>% hist(main='h')
  # title(main=pr, outer = T)
# }

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
