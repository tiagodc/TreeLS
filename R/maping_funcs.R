#' Hough transformation plot-wise
#' @description Applies the hough transformation on TLS plot, returning ALL circles found
#' @param raster output from \code{\link{makeRaster}}
#' @param rad interval of radii to fit circles
#' @param nRad number of radii to test in between \code{rad}
#' @param nAng number of points in a circle
#' @param min.val minimum pixel density value to be included in the circle fitting procedure
#' @param votes minimum number of votes for a circle to be included in the output
#' @param Plot if TRUE, saves a .png file showing all sample points generated for the input raster
#' @param img.prefix file name for Plot
#' @param ... arguments passed to \code{\link{plot}}
#' @return matrix of all circle centers found
#' @export
hough_plot = function(raster, rad = c(.025,.5), nRad = 50, nAng = 90, min.val = .1, votes=3 ,Plot = F, img.prefix = '', ...){

  rads = seq(rad[1],rad[2],length.out = nRad)
  angs = seq(0, 2*pi, length.out = nAng)

  combs = expand.grid(rads, angs)

  x = sin(combs[,2]) * combs[,1]
  y = cos(combs[,2]) * combs[,1]

  coords = cbind(x,y)
  index = rep(1:nRad, times = nAng)

  circles = split(as.data.frame(coords), index)

  survey = raster[[3]]
  survey[survey < min.val] = 0

  means.x = ( raster[[1]][-1] + raster[[1]][-length(raster[[1]])] ) / 2
  means.y = ( raster[[2]][-1] + raster[[2]][-length(raster[[2]])] ) / 2

  where = which(survey > 0, arr.ind = T)
  triang = cbind(means.x[where[,1]], means.y[where[,2]])

  circs=matrix(ncol=2,nrow=0)
  main.base = matrix(ncol = ncol(survey),nrow=0)
  base.index = c()
  rad.index = c()
  centers = matrix(ncol=3,nrow=0)
  for(i in 1:nRad){
    base = survey
    base[base>=0] = 0

    for(j in 1:nrow(triang)){
      temp = t(apply(circles[[i]], 1, function(u) u + triang[j,] ))

      cx = cut(temp[,1], breaks = raster[[1]])
      cy = cut(temp[,2], breaks = raster[[2]])

      tab.temp = table(cx,cy)
      tab.temp[tab.temp > 0] = 1
      base = base + tab.temp
      circs = rbind(circs, temp)
      rad.index = c(rad.index, rep(i, nrow(temp)))
    }
    if(Plot){
      png(filename = paste(img.prefix, i, '.png', sep = ''), width = 10, height = 10, units = 'cm',res = 200)
      par(mar=rep(2,4))
      plot(circs[rad.index == i,], main= round(rads[i], 2), pch=20, cex=.3, ...)
      dev.off()
    }
    main.base = rbind(main.base, base)
    base.index = c(base.index, rep(i, nrow(base)))

    if(sum(base) == 0) next

    bullseye = which(base >= votes, arr.ind = T)

    if(nrow(bullseye) == 0) next

    if(nrow(bullseye) == 1){
      x.cen = means.x[bullseye[1]]
      y.cen = means.y[bullseye[2]]
      xy.cen = c(x.cen, y.cen)
      radii = rads[i]
      append = c(xy.cen, radii)
    }else{
      x.cen = means.x[bullseye[,1]]
      y.cen = means.y[bullseye[,2]]
      xy.cen = cbind(x.cen, y.cen)
      radii = rep(rads[i], nrow(xy.cen))
      append = cbind(xy.cen, radii)
    }

    centers = rbind(centers, append)

  }
  return(centers)
}

#' Map stem coordinates in a TLS point cloud using the Hough transformation
#' @description Outputs the approximate stem positions in a TLS point cloud, estimated by circular patterns found through the hough transformation
#' @param xyz.cloud forest plot point cloud - \emph{xyz} matrix
#' @param z.lim height interval from which a slice of the point cloud will be taken
#' @param min.votes minimum number of votes for a circle to be used in the method
#' @param pixel.size pixel side length - passed to \code{\link{makeRaster}}
#' @param min.den minimum pixel density value to be surveyed in \code{\link{hough_plot}}
#' @param map.rad radius to aggregate circle centers to calculate a stem's coordinate
#' @param Plot if TRUE plots the cloud's slice with the estimated stem cooridinates
#' @param ... arguments passed to \code{\link{plot}}
#' @return xy matrix with all estimated stem coordinates
#' @export
map_HT = function(xyz.cloud, z.lim = 1:2, min.votes=2, pixel.size=.05, min.den=.1, map.rad = 1.5, Plot = T, ...){

  slice = xyz.cloud[ xyz.cloud[,3] >= z.lim[1] & xyz.cloud[,3] <= z.lim[2] ,]

  ras = makeRaster(slice, cell.size = pixel.size, image = F)
  ana = hough_plot(ras, nRad=50, nAng=90, votes = min.votes, min.val = min.den ,Plot = F)

  remaining = ana
  pos = matrix(nrow=0, ncol=2)
  while(nrow(remaining) > 0){

    pt = remaining[1,-3]
    dst = apply(remaining[,-3,drop=F], 1, function(u){ sqrt(sum((u-pt)^2)) })

    ind = which(dst <= map.rad)

    xy = apply(remaining[ind,-3,drop=F],2,mean)
    pos = rbind(pos, xy)

    remaining = remaining[-ind,,drop=F]

  }

  if(Plot){
    plot(slice[,-3], pch=20, cex=.5, ...)
    points(pos, col='red', cex=3, pch=3)
  }

  return(pos)

}

#' Map stem coordinates in a TLS point cloud through spectral decomposition filter
#' @description Outputs the approximate stem positions in a TLS point cloud, estimated from surface flatness and verticality obtained through spectral decomposition
#' @param xyz.cloud forest plot point cloud - \emph{xyz} matrix
#' @param z.lim height interval from which a slice of the point cloud will be taken
#' @param ppm points per meter to retain across the largest dimension (x or y)
#' @param k number of points to analyze in a neighborhood
#' @param flat.min minimum tolerated flatness value to use a point neighborehood in the analysis
#' @param ang.tol maximum angle deviation with \emph{z} tolerated to use a point neighborehood in the analysis
#' @param map.rad radius to aggregate point neighborhoods to calculate a stem's coordinate
#' @param Plot if TRUE plots the cloud's slice with the estimated stem cooridinates
#' @param ... arguments passed to \code{\link{plot}}
#' @return xy matrix with all estimated stem coordinates
#' @export
map_SD = function(xyz.cloud, z.lim = 1:2, ppm = 2000, k=30, flat.min=.9, ang.tol=20, map.rad=1.5, Plot=T, ...){

  slice = xyz.cloud[xyz.cloud[,3] >= z.lim[1] & xyz.cloud[,3] <= z.lim[2],]
  lens = apply(slice,2,function(u) diff(range(u)))
  mx = which(lens[1:2] == max(lens[1:2]))

  temp = slice
  temp[,mx] = slice[,3]
  temp[,3] = slice[,mx]
  temp = Vfilter(temp, 1, ppm)

  tl = Vsections(temp, l.int = 1, overlap = 1/3, Plot = F)

  tl = lapply(tl, function(u){

    if(!is.data.frame(u)) out = NULL else {

      slice = u
      slice[,3] = u[,mx]
      slice[,mx] = u[,3]

      out = SD.prefilt(slice, k, flat.min, ang.tol)
    }
    return(out)
  })

  map = do.call(rbind, tl)
  map = unique(map)
  row.names(map) = NULL

  coor = matrix(nrow = 0, ncol = 3)
  while(nrow(map)>0){
    first = map[1,-3]

    dst = apply(map[,-3], 1, function(u) sqrt(sum((u-first)^2)) )
    i = which(dst <= map.rad)

    if(length(i) < 3*k){ map=map[-i,] ; next }

    tmp = map[i,]
    coor = rbind(coor, apply(tmp, 2, mean))

    map=map[-i,]

  }

  coor2 = apply(coor, 1, function(u){
    a = apply(coor[,-3], 1, function(v) sqrt(sum((u[-3]-v)^2)))
    i = which(a <= map.rad)
    b = colMeans(coor[i,,drop=F])
    return(b)
  })
  coor = unique(t(coor2))

  if(Plot){
    plot(slice, pch=20, cex=.5, ...)
    points(coor, col='red', cex=3, pch=3)
  }

  return(coor)

}

