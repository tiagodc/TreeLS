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

las = pointMetrics(las, ptm.knn())
thin = tlsSample(las, voxelize(.025))

map1 = treeMap(las, map.eigen.knn(mds = .1)) #eigen.knn(mds = .1))
# map2 = treeMap(thin, map.hough()) #eigen.knn(mds = .1))
eigen_map_ento.rds


las = treePoints(las, map1, trp.clip(2, F))

las = stemPoints(las, stm.hough(max_radius = .3, hstep = .3))
# las@data$Stem = with(las@data, Votes > 3 & VotesWeight > .2)

df = stemSegmentation(las, sgmt.ransac.circle())

tlsPlot(las, df) ; pan3d(2)

# pal = las@data$TreeID %>% unique %>% length %>% pastel.colors
# plot(las, color='TreeID', colorPalette=pal)

# plot(map2, clear_artifacts=F)
# rgl.points(map1@data[,.(X,Y,Z)])
# rgl.points(thin@data[,.(X,Y,Z)], size=.5, color='grey')
# pan3d(2)


#####################################

inv = read.csv('../case_tls/inv_entomologia.csv', sep=';', dec=',')
names(inv)[1] = 'FID'

inv = inv[inv$CAP > 0,]

df = stemSegmentation(las, sgmt.irls.cylinder())
tlsInv = df[AvgHeight > 1 & AvgHeight < 1.6, .(dbh = 200*mean(Radius)), by=TreeID]
tlsInv = tlsInv[TreeID %in% inv$Rmap,]

allInv = merge(inv, tlsInv, by.x='Rmap', by.y='TreeID', all=T)
allInv$DAP %>% mean
allInv$dbh %>% mean(na.rm=T)
allInv$diff = allInv$dbh - allInv$DAP

allInv = allInv[ !(allInv$Rmap %in% c(-1, 7, 27, 48, 26)) ,]

mae = (sum(allInv$diff, na.rm=T) / nrow(allInv)) %T>% print
rmse = sqrt(sum( allInv$diff ^2, na.rm=T ) / nrow(allInv)) %T>% print

plot(dbh ~ DAP, data=allInv, cex=0)
text(allInv$DAP, allInv$dbh, allInv$Rmap)
abline(0,1,col='red',lwd=2)

100 * rmse / mean(allInv$DAP)
100 * mae / mean(allInv$DAP)
