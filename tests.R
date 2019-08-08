Rcpp::sourceCpp('src/r_interface.cpp')
source('R/methods.R')
source('R/sample_methods.R')
source('R/map_methods.R')
source('R/stem_points_methods.R')
source('R/stem_segmentation_methods.R')
source('R/tree_points_methods.R')
source('R/point_metrics_methods.R')

require(magrittr)
require(lidR)
require(rgl)
require(data.table)
require(dismo)
require(deldir)
require(nabor)

rm(list = c('.', 'X', 'Y', 'Z', 'Classification', 'TreePosition', 'TreeID', 'Stem', 'Segment', 'gpstime', 'AvgHeight', 'Radius'))

###################

las = readTLS('test_data/ento_u_clip.laz') %>% tlsTransform(c('-x','z','y'), T, T) %>% tlsNormalize(keep_ground = F)

# las = pointMetrics(las, ptm.knn())
# thin = tlsSample(las, voxelize(.025))

# map1 = treeMap(las, map.eigen.knn(mds = .1)) #eigen.knn(mds = .1))
# map2 = treeMap(thin, map.hough()) #eigen.knn(mds = .1))


# las = stemPoints(las, stm.hough(max_radius = .3, hstep = .3))
# las = stemPoints(las, stm.eigen.knn(hstep = .3, max_d = .6))
# las@data$Stem = with(las@data, Votes > 3 & VotesWeight > .2)

# df = stemSegmentation(las, sgmt.ransac.circle())

# tlsPlot(las, df) ; pan3d(2)

# pal = las@data$TreeID %>% unique %>% length %>% pastel.colors
# plot(las, color='TreeID', colorPalette=pal)

# plot(map2, clear_artifacts=F)
# rgl.points(map1@data[,.(X,Y,Z)])
# rgl.points(thin@data[,.(X,Y,Z)], size=.5, color='grey')
# pan3d(2)


#####################################

map1 = readRDS('../case_tls/eigen_map_ento.rds')
las = treePoints(las, map1, trp.clip(2, F))
las = pointMetrics(las, ptm.knn())

inv = read.csv('../case_tls/inv_entomologia.csv', sep=';', dec=',')
names(inv)[1] = 'FID'

inv = inv[inv$CAP > 0,]

gps = c(.5, .33)
stmpts = c('hough', 'eigen knn')
dbhmode = c('section extraction', 'stem segmentation')
dbhshape = c('cylinder', 'circle')
dbhseg = c('irls', 'ransac', 'nm', 'qr')
dbhint = list(c(1,1.6), c(1.1, 1.5), c(1.2,1.4))
nransac = c(5,10,15,20)
iterransac = c(1, 10, 30, 50, 100)
dbhdiff = c(999,4,3,2)
dbhstat = c('mean', 'median')

# h = gps[1]
# i = stmpts[1]
# j = dbhmode[1]
# k = dbhshape[2]
# l = dbhseg[1]
# m = dbhint[[1]]
# n = nransac[3]
# o = iterransac[1]
# p = dbhdiff[1]
# q = dbhstat[1]

results = list()
count = 1

for(h in gps){
  paste('gps time filter:', h, '\n') %>% cat

  gpslas = gpsTimeFilter(las, to=h)

  for(i in stmpts){
    paste('.. stem classification:', i, '\n') %>% cat

    if(i == 'hough'){
      gpslas = stemPoints(gpslas, stm.hough(hstep = .3, max_radius = .3, hbase = c(1.5,3)))
    }else{
      gpslas = stemPoints(gpslas, stm.eigen.knn(.3, max_d = .6))
    }

    for(j in dbhmode){
      paste('.. .. dbh origin:', j, '\n') %>% cat

      for(k in dbhshape){
        if(h == .5 & k == "cylinder") next
        paste('.. .. .. shape fit:', k, '\n') %>% cat

        for(l in dbhseg){
          if(j == 'stem segmentation' & l %in% c('nm', 'qr')) next
          if(k == 'cylinder' & l == 'qr') next
          paste('.. .. .. .. algorithm:', l, '\n') %>% cat

          for(m in dbhint){
            if(j == 'stem segmentation' & m[[1]] > 1) next
            paste('.. .. .. .. .. Z interval:', paste(m, collapse = ' - '), '\n') %>% cat

            for(n in nransac){
              if(l != 'ransac' & n > 5) next
              paste('.. .. .. .. .. .. ransac n:', n, '\n') %>% cat

              for(o in iterransac){
                if(j == 'stem segmentation' & o > 1) next
                if(l != 'ransac' & o > 1) next
                if(k == 'cylinder' & o > 30) next
                paste('.. .. .. .. .. .. .. ransac iterations:', o, '\n') %>% cat

                if(j == 'stem segmentation'){
                  if(k == 'circle'){
                    mtd = if(l == 'irls') sgmt.irls.circle() else sgmt.ransac.circle(n = n)
                  }else{
                    mtd = if(l == 'irls') sgmt.irls.cylinder() else sgmt.ransac.cylinder(n = n)
                  }

                  segdata = stemSegmentation(gpslas, mtd)
                  dbhdata = segdata[AvgHeight > m[1] & AvgHeight < m[2], .(d = mean(Radius)*200), by= TreeID]

                }else{
                  ids = unique(gpslas@data$TreeID)
                  ids = ids[ids %in% inv$Rmap]
                  mtd = if(k == 'circle') circleFit else cylinderFit
                  segdata = lapply(ids, function(x) lasfilter(gpslas, Z > m[[1]] & Z < m[[2]] & Stem & TreeID == x))

                  ids = ids[sapply(segdata, function(x) nrow(x@data)) >= 3]

                  if(o == 1){
                    dbhdata = lapply(segdata, function(x) mtd(x, method=l, n=n)) %>% do.call(what = rbind) %>% as.data.table
                    dbhdata$TreeID = ids
                  }else{
                    dbhlist = lapply(segdata, function(x){ lapply(1:o, function(y) mtd(x, method=l, n=n)) %>% do.call(what = rbind) })
                  }
                }

                for(q in dbhstat){
                  if(j == 'stem segmentation' & q == 'median') next
                  if(l != 'ransac' & q == 'median') next
                  if(o == 1 & q == 'median') next
                  paste('.. .. .. .. .. .. .. .. dbh center:', q, '\n') %>% cat

                  if(o > 1 & j == 'section extraction'){
                    dbhdata = lapply(dbhlist, function(x) apply(x,2,get(q))) %>% do.call(what = rbind) %>% as.data.table
                    dbhdata$TreeID = ids
                  }

                  for(p in dbhdiff){
                    paste('.. .. .. .. .. .. .. .. .. outlier cut:', p, '\n') %>% cat

                    tlsInv = merge(inv, dbhdata, by.x='Rmap', by.y='TreeID')
                    tlsInv$diff = with(tlsInv, d - DAP)
                    tlsInv = tlsInv[abs(tlsInv$diff) < p,]

                    nf = nrow(inv)
                    nr = nrow(tlsInv)
                    np = nr / nf
                    mae  = mean(tlsInv$diff)
                    rmse = sqrt( sum(tlsInv$diff^2) / nr )

                    fieldMean   = mean(tlsInv$DAP)
                    fieldMedian = median(tlsInv$DAP)
                    lidarMean   = mean(tlsInv$d)
                    lidarMedian = median(tlsInv$d)

                    pmae  = 100 * mae / fieldMean
                    prmse = 100 * rmse / fieldMean

                    tempResult = data.frame(
                      gpstime_filter = h,
                      stem_classification = i,
                      dbh_pick = j,
                      shape = k,
                      algorithm= l,
                      height_interval = m[[2]] - m[[1]],
                      ransac_n = n,
                      ransac_iter = o,
                      outlier_cut = p,
                      ransac_summary = q,
                      field_mean = fieldMean,
                      field_median = fieldMedian,
                      lidar_mean = lidarMean,
                      lidar_median = lidarMedian,
                      n_trees = nf,
                      n_trees_lidar = nr,
                      n_trees_ratio = np,
                      rmse = rmse,
                      prmse = prmse,
                      mae = mae,
                      pmae = pmae
                    )

                    results[[count]] = list(summary = tempResult, raw = tlsInv)
                    count = count+1

                    paste('\n\n\n',
                          'time stamp percentile:    ', h, '\n',
                          'Stem classification:      ', i, '\n',
                          'DBH measure process:      ', j, '\n',
                          'Shape pattern:            ', k, '\n',
                          'Shape fitting algorithm:  ', l, '\n',
                          'DBH search height:        ', paste(m, collapse = '-'), 'm \n',
                          'RANSAC samples:           ', n, '\n',
                          'RANSAC trials:            ', o, '\n',
                          'outlier cut:              ', p, '\n',
                          'RANSAC summary statistic: ', q, '\n',
                          '\nmean: \n',
                          '... field:                ', fieldMean %>% round(2), '\n',
                          '... LiDAR:                ', lidarMean %>% round(2), '\n',
                          '\nmedian: \n',
                          '... field:                ', fieldMedian %>% round(2), '\n',
                          '... LiDAR:                ', lidarMedian %>% round(2), '\n',
                          '\nn trees: \n',
                          '... field:                ', nf, '\n',
                          '... LiDAR:                ', nr, '\n',
                          '... ratio:                ', np %>% round(2), '\n',
                          '\nRMSE: \n',
                          '... absolute:             ', rmse %>% round(2), '\n',
                          '... relative (%):         ', prmse %>% round(2), '\n',
                          '\nMAE: \n',
                          '... absolute:             ', mae %>% round(2), '\n',
                          '... relative (%):         ', pmae %>% round(2), '\n\n',
                          Sys.time(),
                    '\n\n\n') %>% cat

                    par(mfrow=c(1,3), oma=c(0,0,1,0))

                    plot(tlsInv$d ~ tlsInv$DAP, pch=20, main='', xlab='Field DBH (cm)', ylab='LiDAR DBH (cm)', cex=.25)
                    text(tlsInv$DAP, tlsInv$d, tlsInv$Rmap)
                    abline(0, 1, col='red', lwd=2)

                    hist(tlsInv$diff, main='', xlab='DBH residuals (cm)')

                    maxd = tlsInv$DAP %>% c(tlsInv$d) %>% max
                    brk = seq(0, maxd + 2, by=2)

                    hist(tlsInv$DAP, brk, col=rgb(1,0,0,.25), xlim=c(5,maxd), ylim=c(0,.2), freq=F, xlab='DBH (cm)', main='')
                    hist(tlsInv$d, brk, col=rgb(0,0,1,.25), add=T, freq=F)

                    fieldQuant = qnorm(1:999 / 1000, mean(tlsInv$DAP), sd(tlsInv$DAP))
                    fieldDen   = dnorm(fieldQuant, mean(tlsInv$DAP), sd(tlsInv$DAP))
                    lines(fieldQuant, fieldDen, col='red', lwd=2)

                    lidarQuant = qnorm(1:999 / 1000, mean(tlsInv$d), sd(tlsInv$d))
                    lidarDen   = dnorm(lidarQuant, mean(tlsInv$d), sd(tlsInv$d))
                    lines(lidarQuant, lidarDen, col='blue', lwd=2)

                    legend('topright', fill = c('red','blue'), legend = c('Field', 'LiDAR'), cex=1.2)
                    title(main=paste(l, k), outer=T)

                    gc()

                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
