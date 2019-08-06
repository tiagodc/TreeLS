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
map1 = readRDS('../case_tls/eigen_map_ento.rds')

# las = pointMetrics(las, ptm.knn())
# thin = tlsSample(las, voxelize(.025))

# map1 = treeMap(las, map.eigen.knn(mds = .1)) #eigen.knn(mds = .1))
# map2 = treeMap(thin, map.hough()) #eigen.knn(mds = .1))

las = treePoints(las, map1, trp.clip(2, F))

# las = stemPoints(las, stm.hough(max_radius = .3, hstep = .3))
las = stemPoints(las, stm.eigen.knn(hstep = .3, max_d = .6))
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

df = stemSegmentation(las, sgmt.irls.circle())
tlsInv = df[AvgHeight > 1 & AvgHeight < 1.6, .(dbh = 200*mean(Radius)), by=TreeID]
tlsInv = tlsInv[TreeID %in% inv$Rmap,]

allInv = merge(inv, tlsInv, by.x='Rmap', by.y='TreeID', all=T)
allInv$DAP %>% mean
allInv$dbh %>% mean(na.rm=T)
allInv$diff = allInv$dbh - allInv$DAP

# allInv = allInv[ !(allInv$Rmap %in% c(-1, 7, 27, 48, 26)) ,]

mae = (sum(allInv$diff, na.rm=T) / nrow(allInv)) %T>% print
rmse = sqrt(sum( allInv$diff ^2, na.rm=T ) / nrow(allInv)) %T>% print

plot(dbh ~ DAP, data=allInv, cex=0)
text(allInv$DAP, allInv$dbh, allInv$Rmap)
abline(0,1,col='red',lwd=2)

100 * rmse / mean(allInv$DAP)
100 * mae / mean(allInv$DAP)

dbhClouds = list()
tids = inv$Rmap

for(i in tids){
  print(i)
  if(i == -1) next
  dbhClouds[[i]] = lasfilter(las, TreeID == i & Z < 1.4 & Z > 1.2)
}

dbhEsts = list()
var = 'V3'
for(i in 1:length(dbhClouds)){
  temp = dbhClouds[[i]]
  if(is.null(temp)) next
  print(i)
  temp %<>% lasfilter(Stem)

  # pars = lapply(1:100, function(x) cppCircleFit(temp %>% las2xyz, 'ransac', n = 15)) %>% do.call(what=rbind) %>% as.data.frame
  # pars$V5 = i
  # pars[,var] = pars[,var] * 200

  pars = cppCylinderFit(temp %>% las2xyz, 'ransac', n = 15)
  pars[5] = pars[5] * 200
  pars = c(pars, i)
  dbhEsts[[i]] = pars #%>% apply(2,median)
}

dbhEsts %<>% do.call(what=rbind) %>% as.data.frame

inv$Rmap
n = 32
inv[inv$Rmap == n,]
dbhEsts[[n]]$V3 %>% hist
dbhEsts[[n]]$V3 %>% mean
dbhEsts[[n]]$V3 %>% median

plot(Y~X, lasfilter(dbhClouds[[n]], Stem)@data, pch=20, cex=.5, asp=1)
plot(dbhClouds[[n]], color="Stem")

var='V5'
allInv = merge(inv, dbhEsts, by.x='Rmap', by.y='V7', all=T)
allInv = allInv[ !(allInv$Rmap %in% c(26, 7, 27, 23)) ,]
# allInv[,var] = allInv[,var]*200
allInv$DAP %>% mean
allInv[,var] %>% mean(na.rm=T)
allInv[,var] %>% median(na.rm=T)
allInv$diff = allInv[,var] - allInv$DAP

plot(allInv[,var] ~ allInv$DAP, cex=0)
text(allInv$DAP, allInv[,var], allInv$Rmap)
abline(0,1,col='red',lwd=2)

hist(allInv$diff)

mae = (sum(allInv$diff, na.rm=T) / nrow(allInv)) %T>% print
rmse = sqrt(sum( allInv$diff ^2, na.rm=T ) / nrow(allInv)) %T>% print
