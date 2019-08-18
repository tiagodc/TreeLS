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

j = 7
i = projs[j]

pfiles = files[grep(i, files)]

k = 5
las = readTLS(pfiles[k] %T>% print)
las = tlsTransform(las, bring_to_origin = T, rotate = T)
las = tlsNormalize(las, keep_ground = F)

las = pointMetrics(las, ptm.knn())

thin = tlsSample(las, voxelize(.025))
map = treeMap(thin, map.hough(hmin=1, hmax=4, max_d=.35, min_density = .075, min_votes = 4, hstep = .5))
tp = treePositions(map)

map = treeMap(las, map.eigen.knn(max_d=.35, pln = .17, vrt = 16, min_n = 100, min_h = 1, mds = .1))
tp = treePositions(map)

plot(Y ~ X, data=tp, asp=1, cex=1, pch=20)
text(tp$X, tp$Y, tp$TreeID)

plot(las, clear_artifacts=F, size=.5) ; pan3d(2) ; axes3d(col='white')

rmTrees = c(115, 128, 101, 29, 42, 46, 82, 53, 88, 94, 66, 19, 15, 14, 12, 20, 58, 49, 77, 107, 32)
map@data = map@data[!(TreeID %in% rmTrees)]
tp = treePositions(map)

las = treePoints(las, map, trp.crop(r = 2, circle = F))

col = las$TreeID %>% unique %>% length %>% pastel.colors
plot(las %>% lasfilter(TreeID > 0), color='TreeID', colorPalette=col, size=.5) ; pan3d(2)

# las = stemPoints(las, stm.hough(.3, max_radius = .15, hbase = c(2.5,5)))
las = stemPoints(las, stm.eigen.knn(.3, max_d=.35, dvt = .33))
plot(las %>% lasfilter(Stem), color='VotesWeight')

dbhSegs = lasfilter(las, Stem & Z > 1.2 & Z < 1.4)
plot(dbhSegs)
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
saveRDS(las, '../../plantar/results/401P04_h_stem.rds')
saveRDS(map, '../../plantar/results/401P04_h_map.rds')
write.csv(estimates, '../../plantar/results/401P04_dbh.csv', quote = F, row.names = F)
