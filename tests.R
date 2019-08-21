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

las = readRDS('../plantar/results/399P02_h_stem.rds')
las = pointMetrics(las, ptm.knn(10))
las = stemPoints(las, stm.eigen.knn(.3, max_d=.3))

dbhSeg = lasfilter(las, Z > 1.2 & Z < 1.4 & TreeID > 0 & Stem)
plot(dbhSeg)
dbhSeg = lapply(dbhSeg$TreeID %>% unique, function(x) lasfilter(dbhSeg, TreeID == x))

estimates = data.frame()
for(i in dbhSeg){
  est = circleFit(i, 'ransac', n_best=30)
  tlsPlot.dh(i, est)
  estimates %<>% rbind(est)
}

for(id in unique(las@data[TreeID > 0]$TreeID)){
# 6,10,12,14,26
# id = 26
tree = lasfilter(las, TreeID == id)
# plot(tree, size=.5) ; pan3d(2)

dbh = lasfilter(tree, Z > 1.2 & Z < 1.4)
# plot(dbh, color='gpstime')

px = .02; dmax = .2; p = .67; d = .25

  hg = getHoughCircle(dbh %>% las2xyz, px, rad_max = dmax/2, min_den = d) %>% do.call(what=rbind) %>% as.data.table
  names(hg) = c('x','y','r','v')
  hg = hg[v > quantile(v, p)]
  hg$clt = 1

  houghClusters = hg
  centers = hg[v == max(v)] %>% apply(2,mean) %>% t %>% as.data.table

  k = 2
  repeat{
    km = kmeans(hg[,1:3], k)
    hg$clt = km$cluster
    mxs = hg[,.(x=mean(x), y=mean(y), r=mean(r), v=mean(v)), by=clt]
    dst = mxs[,c('x','y')] %>% dist
    combs = combn(nrow(mxs), 2)
    rst = apply(combs, 2, function(x){mxs$r[x[1]] + mxs$r[x[2]]}) - px
    isForked = all(min(dst) > rst)
    # isForked = all(abs(min(dst) - rst) < 2*px)

    if(!isForked) break

    houghClusters$clt = km$cluster
    centers = mxs[,.(x=mean(x),y=mean(y),r=mean(r),v=mean(v)),by=clt][order(clt)]
    centers$ssRatio = (km$withinss / km$size) / min(km$withinss / km$size)
    centers$nRatio = km$size / max(km$size)
    centers$absRatio = centers$ssRatio / centers$nRatio
    centers = centers[absRatio < 5][order(-v)]
    k=k+1
  }

  plot(dbh$Y ~ dbh$X, cex=.5, asp=1, pch=20, main=id, ylab='Y (m)', xlab='X (m)')
  vcols = lidR:::set.colors(houghClusters$v, height.colors(houghClusters$v %>% unique %>% length))
  vcols = lidR:::set.colors(houghClusters$clt, height.colors(houghClusters$clt %>% unique %>% length))
  points(houghClusters$x, houghClusters$y, col=vcols, pch=20, cex=1)
  #print(centers)

  dbh@data$StemID = 0
  # tree@data$StemID = 0
  for(i in 1:nrow(centers)){
    temp = centers[i]
    angs = seq(0,pi*2,length.out = 36)
    x = temp$x + cos(angs) * temp$r
    y = temp$y + sin(angs) * temp$r
    lines(x,y,col='orange',lwd=2)

    dsts = dbh@data %$% sqrt( (X - temp$x)^2 + (Y - temp$y)^2 )
    pts = dsts < temp$r + 2*px & dbh@data$StemID == 0
    dbh@data$StemID[pts] = i

    # tds = tree@data %$% sqrt( (X - temp$x)^2 + (Y - temp$y)^2 )
    # tps = tds < temp$r + 2*px & tree@data$StemID == 0
    # tree@data$StemID[tps] = i

    cld = lasfilter(dbh, StemID == i)
    est = circleFit(cld, 'ransac', n=15, inliers = .7, n_best = 50)
    x = est$X + cos(angs) * est$d/200
    y = est$Y + sin(angs) * est$d/200
    lines(x,y,col='green',lwd=2)
  }

}



