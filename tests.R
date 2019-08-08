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
las = treePoints(las, map1, trp.clip(2, F))
las = pointMetrics(las, ptm.knn())

las = stemPoints(las, stm.eigen.knn(.3, max_d = .6))
plot(las, color="Stem")

s1 = readRDS('../case_tls/simulations.rds')
s2 = readRDS('../case_tls/simulations_2.rds')

simList = s1
simList[(length(s1)+1):(length(s1) + length(s2))] = s2
simNames = sapply(simList, function(x) x[[1]][,1:10] %>% lapply(as.character) %>% paste(collapse='_') %>% gsub(pattern = '\\s', replacement = ""))
simList = simList[!duplicated(simNames)]
simNames = simNames[!duplicated(simNames)]

names(simList) = simNames

simulations = c(lapply(s1, function(x) x[[1]]), lapply(s2, function(x) x[[1]])) %>% do.call(what=rbind) %>% as.data.table %>% unique

simulations$n_trees_ratio = simulations$n_trees_lidar / 36

temp = simulations[outlier_cut > 2 & gpstime_filter < 2]
temp = temp[order(rmse / n_trees_ratio)] %>% head(20) %>% print
temp$n_trees_ratio %>% mean
temp$rmse %>% mean
temp$mae %>% mean

# gpstime_filter
stem_classification = 'eigen knn'
dbh_pick = 'section extraction'
# shape
algorithm = 'ransac'
height_interval < .6
ransac_n >= 10
ransac_iter >= 10
# outlier_cut
# ransac_summary

# i = 1
# hash = temp[i,1:10] %>% sapply(as.character) %>% paste(collapse = '_') %>% gsub(pattern = '\\s', replacement =  '')
# hashData = simList[[hash]]

inv = read.csv('../case_tls/inv_entomologia.csv', sep=';', dec=',')
names(inv)[1] = 'FID'
inv = inv[inv$CAP > 0,]

ids = unique(las@data$TreeID)
ids = ids[ids %in% inv$Rmap]

dbhClouds = lapply(ids, function(x) lasfilter(las, Stem & TreeID == x & Z > 1.2 & Z < 1.4 & VotesWeight > .33))
names(dbhClouds) = ids

dbhEsts = lapply(dbhClouds, function(x) lapply(1:10, function(y) cylinderFit(x, method='ransac', n=15, inliers = .9)) %>% do.call(what=rbind) %>% apply(2,mean)) %>% do.call(what=rbind) %>% as.data.table
# dbhEsts = lapply(dbhClouds, cylinderFit, method='ransac', n=15) %>% do.call(what=rbind)
dbhEsts$TreeID = names(dbhClouds) %>% as.double

tlsInv = merge(inv, dbhEsts, by.x='Rmap', by.y='TreeID') %>% as.data.table
tlsInv$diff = tlsInv$d - tlsInv$DAP

tlsInv[order(-abs(diff))]

plot(d ~ DAP, data=tlsInv, cex=.3, pch=20)
abline(0,1,col='red',lwd=2)
text(tlsInv$DAP, tlsInv$d, tlsInv$Rmap)

for(i in 1:nrow(tlsInv)){
  row = tlsInv[i,]
  cld = dbhClouds[[ row[[1]] %>% as.character ]]
  plot(Y~X, data=cld@data,pch=20,cex=1, asp=1, main=paste('tree',row$Rmap))
  points(row$X, row$Y, col='red', cex=2, pch=3)
  xp = cos(seq(0, 2*pi, length.out = 36)) * (row$d/200) + row$X
  yp = sin(seq(0, 2*pi, length.out = 36)) * (row$d/200) + row$Y
  lines(xp, yp, col='red', lwd=2)
}

mean(tlsInv$DAP)
(mean(tlsInv$d) + median(tlsInv$d))/2
mean(tlsInv$diff)
sqrt(sum(tlsInv$diff^2)/nrow(tlsInv))

100 * length(tlsInv[abs(diff) < 4]$diff) / 36

outliers = c(7,23,26,27,35,50)
tlsInv[Rmap %in% outliers]

lastemp = lasfilter(las, Z > 1.2 & Z < 1.4 & TreeID == 50 & Stem & VotesWeight > .55)

# lastemp@data %<>% merge(tlsInv, by.x='TreeID', by.y='Rmap', sort=F)
# lastemp@data$absdiff = lastemp@data$diff %>% abs
# names(lastemp@data)[3:4] = c('X','Y')
# plot(lastemp, color='Votes', size=2, clear_artifacts=F)
# pan3d(2)
plot(lastemp@data[,.(X,Y)], pch=20, cex=.5, asp=1)

tlsInv[Rmap == 50]
pars = circleFit(lastemp, 'ransac', n=15, inliers = .9) %T>% print



