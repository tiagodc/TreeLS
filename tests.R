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

las = readTLS('../entomologia_zeb_cubagem.laz') %>% tlsNormalize(.1, F)
plot(las)

dev.off()
par(mfrow=c(2,2))
for(i in seq(min(las$Z), max(las$Z), .5)){
  cyl = lasfilter(las, Z > i & Z < i+.5)
  # a = circleFit(cyl, method = "ransac", 15, .8, .99, 5)
  # a = cylinderFit(cyl %>% toLAS, 'bf', 10, .8, .95, 30)
  a = robustDiameter(cyl, max_d=.5)
  if(is.null(a)) next
  tlsPlot.dh(cyl, a, clear = F)
}

vals = bruteForceRansacCylinder(cyl %>% as.matrix, 10, .99, .7, 5, 25)
vals = do.call(rbind, vals) %>% as.data.table
vals[order(V4),]

angs = seq(-45,45,3)

res = data.table()
for(i in seq(-45,45,3)){
  for(j in seq(-45,45,3)){
    tpc = as.matrix(cyl) %*% rotationMatrix(i * pi/180, j * pi/180, 0)
    pars = circleFit(tpc %>% toLAS, 'ransac', 15, .8, .99, n_best = 5)
    pars$i = i
    pars$j = j
    res %<>% rbind(pars)
  }
}

res[order(abs(err)),]


clear3d() ; axes3d()
plot3d(cyl, aspect = F)
lines3d(c(0,0), c(0,0), c(0,1), col='red', lwd=5)

for(i in seq(-15,15,3)){
  for(j in seq(-15,15,3)){
    ax = (c(0,0,1) %*% rotationMatrix(i * pi/180, j * pi/180, 30))
    lines3d(c(0,ax[1]), c(0,ax[2]), c(0,ax[3]), col='blue', lwd=2)
  }
}

