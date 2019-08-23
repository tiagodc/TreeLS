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

files = dir('/media/tiago/DATA/plantar/results', full.names = T)
projs = sub('(.+?)_.+', '\\1', files) %>% unique

# for(pr in projs[7]){
  pr = projs[10]
  paste('.. project', pr) %>% print

  las = pr %>% paste0('_h_stem.rds') %>% readRDS
  # vlas = pr %>% paste0('_v.rds') %>% readRDS
  res = pr %>% paste0('_res.csv') %>% read.csv %>% as.data.table
  # res = res[,-1]

  dev.off()
  par(mfrow=c(1,2), oma=c(0,0,1,0))
  res$d %>% hist(main='d')
  res$h %>% hist(main='h')
  title(main=pr, outer = T)

  res$d %>% mean
  res$d %>% median
  res$d %>% sd

  las %<>% lasfilter(Z > .5 & Z < 3)
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
  #   paste('.. .. dbh', id) %>% print

    id = 16
    tree = lasfilter(las, TreeID == id)
  #   # clear3d() ; rgl.points(tree@data[,.(X,Y,Z)], size=.5)
  #
    dbh = lasfilter(tree, Z > 1.2 & Z < 1.4)
    # dbh %<>% gpsTimeFilter(from=.6)
    # dbh = tlsCrop(dbh, 2.2,10.4, .5)
    # plot(dbh)
  #
    # est = res[TreeID == id,2:5]
    # par(mfrow=c(1,nrow(est))) ; clear3d()
    # est[,1:4] %>% apply(1, function(x) tlsPlot.dh(dbh, x, clear = F))

    est = robustDiameter(dlas=dbh, plot=T, max_d = .12, pixel_size = .02, min_den = .1) ; est %>% nrow
    # est = est[1,]
    par(mfrow=c(1,nrow(est))) ; clear3d()
    rgl.points(tree@data[,.(X,Y,Z)], size=.5, color='white')
    est[,1:4] %>% apply(1, function(x) tlsPlot.dh(dbh, x, clear = F))
  #   estimates %<>% rbind(est)
  # }

    redone %<>% rbind(est)
    # res %<>% as.data.table
    redone %<>% merge(res[,.(TreeID,h)], by='TreeID')

    res = res[!(TreeID %in% redone$TreeID)]
    res %<>% rbind(redone)
    res %<>% unique
    pr %>% paste0('_res_ok.csv') %>% write.csv(x = res)


    for(i in projs){
      # map = i %>% paste0('_h_map.rds') %>% readRDS
      # file = i %>% paste0('_v.rds') %>% readRDS %>% treePoints(map, trp.crop(2.5))
      # hs = file@data[TreeID > 0,.(h=max(Z)),by=TreeID][order(TreeID)]
      #
      # res = i %>% paste0('_res_ok.csv') %>% read.csv %>% as.data.table
      # res = res[,-1][order(TreeID)]
      # res$h = NULL
      #
      # res = merge(res, hs, by='TreeID')
      # res$h[res$h < 13] = mean(res$h)
      #
      # i %>% paste0('_final.csv') %>% write.csv(x = res, quote = F, row.names = F)
      #
      res = i %>% paste0('_final.csv') %>% read.csv()
      i %>% paste0('_hists.png') %>% png(20,30,'cm',res = 300)
      par(oma=c(0,0,1,0))

      layout(c(1,2,3,3) %>% matrix(ncol=2,byrow = T), widths = c(1,1), heights = c(1,2))
      md = res$d %>% mean %>% round(2)
      hist(res$d, main=md %>% paste('cm'), xlab='DAP (cm)', ylab='Frequência')
      abline(v=md, lwd=3, lty=2)

      mh = res$h %>% mean %>% round(2)
      hist(res$h, main=mh %>% paste('m'), xlab='H (m)', ylab='Frequência')
      abline(v=mh, lwd=3, lty=2)

      idRes = as.data.table(res)[,.(X=mean(X),Y=mean(Y)), by=TreeID]
      plot(Y~X,res,pch=20,cex=.8,main='Mapa da Parcela', xlab='X (m)', ylab='Y (m)', asp=1)
      idRes %$% text(X,Y,TreeID, pos=1, offset=.3)
      legend('topright', legend = 'ID da árvore', pch=20)

      # plot(res$h ~ res$d, ylim=c(0,25), xlim=c(0,25), pch=20, cex=1.5, xlab='DAP (cm)', ylab='H (m)')
      # dh = lm(h~d,data=res)
      # dest = seq(min(res$d), max(res$d), length.out = 20)
      # prd = predict.lm(dh, list(d=dest))
      # lines(dest, prd, col='red', lwd=2)

      title(main= rev(strsplit(i,'/')[[1]])[1] ,outer=T)
      dev.off()
    }


tab = projs %>% paste0('_final.csv') %>% lapply(read.csv)
for(i in 1:length(tab)){
  pt = rev( strsplit(projs[i], '/')[[1]] )[1]
  tab[[i]]$Parcela = pt
  tab[[i]] = tab[[i]][,c('Parcela', 'TreeID', 'h', 'StemID', 'd', 'score', 'X', 'Y')]
  names(tab[[i]]) = c('Parcela', 'Árvore', 'Altura', 'Fuste', 'DAP', 'Confiança', 'X', 'Y')
}

tab = do.call(rbind, tab) %>% unique
write.csv(tab, 'platar.csv', quote = F, row.names = F)
