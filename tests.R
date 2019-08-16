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

files = dir('../../plantar/', full.names = T)

projs = sub('.*/plantar/(\\d+P\\d+).+las$', '\\1', files) %>% unique
projs = projs[ !grepl('^\\.', projs) ]

j = 3
i = projs[j]

pfiles = files[grep(i, files)]

k = 6
las = readTLS(pfiles[k] %T>% print)
las = tlsTransform(las, bring_to_origin = T, rotate = T)
las = tlsNormalize(las, keep_ground = F)

las = pointMetrics(las, ptm.knn())

thin = tlsSample(las, voxelize(.03))
map = treeMap(thin, map.hough(hmin=2, hmax=4, max_d=.35))
tp = treePositions(map)

# map = treeMap(las, map.eigen.knn(max_d=.35, pln = .15, vrt = 15, min_n = 100, min_h = 1, mds = .1))
# tp = treePositions(map)

plot(Y ~ X, data=tp, asp=1, cex=1, pch=20)
text(tp$X, tp$Y, tp$TreeID)

rmTrees = c(122,263,32,120, 185,26, 13)
map@data = map@data[!(TreeID %in% rmTrees)]
tp = treePositions(map)

las = treePoints(las, map, trp.crop(r = 2, circle = F))

col = las$TreeID %>% unique %>% length %>% pastel.colors
plot(las, color='TreeID', colorPalette=col)

las = stemPoints(las, stm.hough(.3, max_radius = .15, hbase = c(2.5,5)))
# las = stemPoints(las, stm.eigen.knn(.3, max_d=.35))
plot(las, color='Stem')

dbhSegs = lasfilter(las, Stem & Z > 1.2 & Z < 1.4)
dbhSegs = split(dbhSegs@data, dbhSegs@data$TreeID) %>% lapply(LAS)

estimates = data.frame()
for(tid in names(dbhSegs)){
  temp = dbhSegs[[tid]]
  est = lapply(1:30, function(x) circleFit(temp, 'ransac', n = 15, inliers = .7, n_best = 20)) %>% do.call(what=rbind) %>% apply(2, mean)
  tlsPlot.dh(temp, est)
  estimates %<>% rbind(c(as.double(tid), est))
}
names(estimates) = c('TreeID', names(est))

estimates$d %>% hist
estimates$d %>% mean
estimates$d %>% median

v = 2
vlas = readTLS(pfiles[v] %T>% print)
plot(vlas)

pfiles[k] %>% print
saveRDS(las, '../../plantar/results/217P03_h_stem.rds')
saveRDS(map, '../../plantar/results/217P03_h_map.rds')
write.csv(estimates, '../../plantar/results/217P03_dbh.csv', quote = F, row.names = F)
