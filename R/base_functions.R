###Data sets

#' Scots Pine point cloud
#' @docType data
#' @name pine
#' @usage pine
NULL

#' Norway Spruce point cloud
#' @docType data
#' @name spruce
#' @usage spruce
NULL

#' Field measurements from sample trees
#' @docType data
#' @name proof
#' @usage proof
NULL

###Data simulation and analysis

#' Angle of a plane with \emph{z}
#' @description Calculates the angle between a plane's normal vector and the \emph{z} axis \code{c(0,0,1)}. The plane here is a point cloud sample.
#' @param XYZplane sample from a point cloud. Matrix n x 3, with columns containing x, y and z coordinates, respectively.
#' @return angle between the plane's normal vector and the \emph{z} axis \code{c(0,0,1)} (in degrees).
#' @export
ang.deg = function(XYZplane){
  e = eigen(cov(XYZplane))
  ang = ( e$vectors[,3] %*% c(0,0,1) ) / ( sqrt(sum(e$vectors[,3]^2)) * sqrt(sum(c(0,0,1)^2)) )
  ang = ang[,,drop=T]
  degs = acos(ang)*180/pi
  return(degs)
}

#' Angle between two vectors
#' @description Calculates the angle between two vectors.
#' @param a,b numerical vectors of same length
#' @return angle between \emph{a} and \emph{b} (in degrees).
#' @export
angle = function(a,b){
  prod = a %*% b
  lprod = sqrt(sum(a^2)) * sqrt(sum(b^2))
  ang = prod/lprod
  cang = acos(ang) * 180/pi
  return(cang[,,drop=T])
}

#' Rotates and dislocates a 3D dataset
#' @description Rotates and shifts the coordinates of a 3D point cloud
#' @param xyz.points point cloud, \emph{xyz} matrix object
#' @param rot.mat rotation matrix
#' @param shift dislocation in x, y and z directions, respectively - vector of length 3
#' @return point cloud with new coordinates
#' @seealso \code{\link{xyz.rotation.matrix}}
#' @export
change.coords = function(xyz.points, rot.mat=diag(3), shift=c(0,0,0)){
  xyz.points = as.matrix(xyz.points)
  rot = xyz.points %*% rot.mat
  mrot = t(apply(rot, 1, function(u) u + shift))
  return(mrot)
}

#' Random cylinder generation
#' @description Creates a random cylinder with specified paramenters
#' @param n number of points on the cylinder's surface
#' @param len length of the cylinder's axis
#' @param d cylinder diameter
#' @param dev deviation of points from surface, i.e. cylinder wall thickness
#' @param top.bottom if \code{TRUE}, creates only two circles with \code{n} points at bottom and top of the cylinder.
#' @return random cylinder point cloud - \emph{xyz} matrix
#' @export
cyl = function(n=10000, len=100, d=30, dev=NULL, top.bottom = F){

  if(is.null(dev)) rad = d/2 else rad = runif(n, d-dev, d+dev)/2

  if(top.bottom){
    angs = seq(0,2*pi, length.out = n)
    x = sin(angs)*rad
    y = cos(angs)*rad
    z = rep(c(0,len), each=n)

  } else {

    z=runif(n = n, min = 0, max = len)

    angs = runif(n, 0, 2*pi)
    x = sin(angs)*rad
    y = cos(angs)*rad

  }
  return(cbind(x,y,z))
}

#' Surface flatness
#' @description calculates the flatness of a point cloud
#' @param XYZplane point cloud - \emph{xyz} matrix
#' @return surface flatness - ranging from 0 (non-flat) to 1 (perfectly flat)
#' @export
FL = function(XYZplane){
  e = eigen(cov(XYZplane))
  flat = 1 - ( e$values[3] / sum(e$values) )
  return(flat)
}

#' Tukey's biweight function
#' @describe calculates weights for a dataset based on Tukey's biweight function
#' @param errors residuals from a fitted model
#' @param b efficiency constant (5 == 95\%)
#' @return list of length 2: \describe{\item{$Y}{residuals / MAD} \item{$weights}{weights}}
#' @export
tukey.estimator = function(errors, b = 5){
  s = mad(errors)
  Y = errors / s

  tue = Y
  larger = abs(tue) > b
  tue[larger] = 0#b^2 / 6
  tue[!larger] = (1-(tue[!larger]/b)^2)^2
  #tue[!larger] = (b^2 / 6) * (1-(1-(tue[!larger]/b)^2 )^3 )

  return(list(Y = Y, weights = tue))
}

#' Cross product between two vectors
#' @description calculates the cross product between two vectors - not to be confused with \code{crossprod}, which gives the dot product
#' @param a,b numerical vectors of same length
#' @return vector cross product
#' @export
Xprod = function(a,b){
  x = c(a[2]*b[3] - a[3]*b[2] ,
        a[3]*b[1] - a[1]*b[3] ,
        a[1]*b[2] - a[2]*b[1])
  return(x)
}

#' 3D rotation matrix
#' @description calculates the 3D rotation matrix for a \emph{xyz} dataset
#' @param ax angle of 1st rotation -  around the \emph{x} axis (in degrees)
#' @param az angle of 2nd rotation -  around the \emph{z} axis (in degrees)
#' @param ax2 angle of 3rd rotation - around the new \emph{x} axis (in degrees)
#' @return 3 x 3 \emph{xyz} rotation matrix
#' @export
xyz.rotation.matrix = function(ax, az, ax2){

  ax = ax * pi/180
  Rx = matrix(c(1,        0,        0,
                0,  cos(ax),  sin(ax),
                0, -sin(ax), cos(ax)),
              ncol=3, byrow = T)

  az = az * pi/180
  Rz = matrix(c(cos(az), 0, -sin(az),
                0, 1,        0,
                sin(az), 0,  cos(az)),
              ncol=3, byrow=T)

  ax2 = ax2 * pi/180
  Rx2 = matrix(c( cos(ax2), sin(ax2), 0,
                  -sin(ax2), cos(ax2), 0,
                  0,       0, 1),
               ncol=3, byrow=T)

  ro.mat = Rx2 %*% Rz %*% Rx
  return(ro.mat)
}


###General point cloud manipulation and visualization

#' Point cloud coloring by height
#' @description colors a point cloud by height intervals
#' @param xyz.cloud point cloud - \emph{xyz} matrix
#' @param pal color pallete - defaults to \code{rainbow}
#' @param n number of height intervals to color differently
#' @param quantile divide point cloud by height quantiles? - If FALSE, divides the point cloud in equal height intervals
#' @import rgl
#' @return plots \emph{xyz.cloud} using \code{rgl} and with the specified color pallete
#' @export
cloud.col = function(xyz.cloud, pal = rainbow, n=10, quantile = F){

  if(quantile){
    nq = quantile(xyz.cloud[,3], probs = seq(0,1,length.out = n+1))
    cl = cut(xyz.cloud[,3], breaks = nq, include.lowest = T)
  }else{
    rg = range(xyz.cloud[,3])
    cl = cut(xyz.cloud[,3], breaks = seq(rg[1], rg[2], by = (rg[2]-rg[1])/n), include.lowest = T)
  }
  cols = pal(n)[cl]
  return(cols)

}

#' Centralize of \emph{xy} to zero
#' @description reassigns the x and y coordinates of a point cloud considering zero as its center
#' @param xyz.cloud point cloud - \emph{xyz} matrix
#' @return point cloud with new coordinates
#' @export
center.zero = function(xyz.cloud){
  mn = colMeans(xyz.cloud)
  names(mn) = c('x','y','z')
  xyz.cloud[,1] = xyz.cloud[,1] - mn[1]
  xyz.cloud[,2] = xyz.cloud[,2] - mn[2]
  return(list(cloud = xyz.cloud, center = mn[1:2]))
}

#' Circular point cloud clip
#' @description clips a region of a point cloud in circular shape
#' @param cloud point cloud - \emph{xyz} matrix
#' @param rad circle radius
#' @param center x and y center coordinates for the circle
#' @return point cloud of all points inside the specified circle
#' @export
clip.XY = function(cloud, rad = 1.5, center = c(0,0)) {
  if (center[1] == 0 & center[2] == 0) {
    dists = sqrt(cloud[,1]^2+cloud[,2]^2)
    out = cloud[dists <= rad, ]
  } else {
    dists = sqrt(((cloud[,1]-center[1])^2 + (cloud[,2]-center[2])^2))
    out=cloud[dists<=rad, ]
  }
  return (out)
}

#' Plot 3D axes from the origin
#' @description plots the 3 main axes (\emph{xyz}) starting at \code{c(0,0,0)}
#' @param xyz length of each axis, x, y and z, respectively
#' @param cols color of x, y and z, respectively
#' @param ... further arguments passed to the \code{rgl.lines} function
#' @import rgl
#' @return plots 3D axes over the current rgl environment
#' @export
rglAXES = function(xyz = c(1,1,1), cols = c('red','green','blue'), ...){
  rgl.lines(c(0,xyz[1]), c(0,0), c(0,0), col=cols[1], ...)
  rgl.lines(c(0,0), c(0,xyz[2]), c(0,0), col=cols[2], ...)
  rgl.lines(c(0,0), c(0,0), c(0,xyz[3]), col=cols[3], ...)
}

## ## Plotting a stem model
## ## @description plots a stem model of stacked cylinders or circles, depending on the \emph{fitting} routine used to calculate the stem segments
## ## @param stem.out output from a stem fitting function - \code{\link{fit_RANSAC_circle}}, \code{\link{fit_RANSAC_cylinder}} or \code{\link{fit_IRTLS}}
## ## @param cyl.len optional - cylinder length for all stem segments
## ## @param col color pallete function or color string name to use to color the cylinders
## ## @param bg background color of the rgl environemnt
## ## @param alpha alpha value passed on to the \code{\link{ashape3d}} function
## ## @examples
## ##\dontrun{
## ## trunk <- pref_HT(spruce)
## ## stem <- fit_RANSAC_circle(trunk)
## ## obj3d <- stem.model(stem)
## ## rgl.points(spruce, size=1)
## ##}
## ## @return 3D stem model of stacked cylinders/circles
## ## @import alphashape3d
## ## @export
## stem.model = function(stem.out, cyl.len=NA, col=rainbow, bg='black', alpha=.5){
##
##   #require(alphashape3d)
##
##   st = stem.out[[1]]
##   ft = stem.out[[2]]
##
##   if(ncol(ft) == 8){
##
##     cln = if(is.na(cyl.len)) ft[,2]-ft[,1] else rep(cyl.len, nrow(ft))
##     cols = if(class(col)=='function') col(nrow(ft)) else rep(col, nrow(ft))
##
##     abs = list()
##     pts = list()
##     for(i in 1:nrow(ft)){
##
##       tp = st[st[,3] <= ft[i,'z2'] & st[,3] >= ft[i,'z1'],]
##
##       vcs = cyl.vectors(ft[i,3:7])
##       d = ft[i,'r']*2
##
##       a = vcs$a
##       #h = if(a[3] < 0) ft[i,'z1'] else ft[i,'z1']
##
##       zang = angle(a, c(0,0,1))
##       xang = angle(c(vcs$n[-3],0), c(1,0,0))
##
##       rot = xyz.rotation.matrix(0,zang,xang)
##
##       cl = cyl(n=1000, len=cln[i], d=d)
##
##       go = change.coords(cl, rot, shift = vcs$Q)
##       ed = sqrt(sum((colMeans(tp) - colMeans(go))^2))
##       if(a[3]<0) ed = -ed
##       go = t(apply(go, 1, function(x) x+a*ed))
##
##       #go[,3] = go[,3] + abs(h)-min(go[,3])
##
##       acl = ashape3d(go, alpha = alpha)
##
##       abs[[i]] = acl
##       pts[[i]] = tp
##
##     }
##
## } else {
##
##   if(ncol(ft) == 6){
##
##   if(is.na(cyl.len)) cln = .02 else cln = cyl.len
##   cols = if(class(col)=='function') col(nrow(ft)) else rep(col, nrow(ft))
##
##   abs = list()
##   pts = list()
##   for(i in 1:nrow(ft)){
##
##     tp = st[st[,3] <= ft[i,'z2'] & st[,3] >= ft[i,'z1'],]
##
##     xy = ft[i,3:4]
##     d = ft[i,'r']*2
##
##     if(d == 0 | d > 2) next
##
##     cl = cyl(n=1000, len=cln, d=d)
##
##     h = mean(ft[i,1:2])
##
##     go = t(apply(cl, 1, function(u) u + c(xy,h)))
##
##     acl = ashape3d(go, alpha = alpha)
##
##     abs[[i]] = acl
##     pts[[i]] = tp
##
##   }
##
##   }}
##
##   nulls = sapply(abs, is.null)
##   nl2 = sapply(pts, is.null)
##
##     bg3d(bg)
##     lapply(1:length(abs), function(u) if(!nulls[u]) plot.ashape3d(abs[[u]], clear=F, edges=F, vertices=F, col=cols[[u]]))
##     lapply(1:length(pts), function(u) if(!nl2[u]) rgl.points(pts[[u]], col=cols[[u]]) )
##
##     return(abs[!nulls])
##
## }

#' Height-based point cloud filter
#' @description reduces a point cloud's density processing different height intervals individually
#' @param XYZtree point cloud (not necessarily for a whole a tree) - \emph{xyz} matrix
#' @param l.int length of height intervals to split the point cloud into
#' @param thr maximum density threshold - i.e. maximum amount of points tolerated per height interval
#' @return point cloud with reduced density
#' @seealso \code{\link{Vsections}}
#' @export
Vfilter = function(XYZtree, l.int = .3, thr = 10000){

  #Description: reduces the point density of a point cloud

  #XYZtree == matrix or data frame with 3 columns containing x, y and z coordinates, respectively
  #l.int == passed to function Vsections. Length of height intervals
  #thr == threshold, maximum number of points to retain in each chunk (radomly selected when n of points > thr)

  #output == point cloud with reduced point density according to specified parameters

  a = Vsections(XYZtree, l.int = l.int, Plot = F)
  a = lapply(a , function(u){ if(nrow(u)>thr) u = u[sample(1:nrow(u), size = thr),] ; return(u) })
  a = do.call(rbind, a)
  return(a)
}

#' Split point cloud into height intervals
#' @description divides a point cloud in many smaller ones, according to height intervals
#' @param XYZtree point cloud (not necessarily a tree) - \emph{xyz} matrix
#' @param n.int number of height intervals (divided by point quantiles)
#' @param l.int optional - length of height intervals. If not NULL, splits a point cloud by fixed height intervals, instead of using \emph{n.int}
#' @param overlap optional - by how much should the height segments overlap?
#' @param Plot create a plot for every point cloud segment? TRUE or FALSE
#' @param units if \code{Plot == TRUE}, provide the units of measurement for labelling the plots
#' @param ... further arguments passed to \code{plot}
#' @return list with every compartment containing a section of the point cloud in fixed \emph{z} intervals
#' @export
Vsections = function(XYZtree, n.int = 100, l.int = NULL, overlap = NULL, Plot =T, units = 'm', ...){

  #Description: subdivides the point cloud in height intervals

  #XYZtree == matrix or data frame with 3 columns containing x, y and z coordinates, respectively
  #n.int == number of height intervals to subdivide the point cloud
  #l.int == if not NULL, n.int is disconsidered and the point cloud is subdivided in equal height invervals of length l.int
  #overlap == proportion of point cloud chunks to overlap in z, applied when overlap != NULL
  #Plot == if TRUE, plots the x and y coordinates of every chunk
  #units == spatial unit, applicable when Plot == TRUE
  #... == further 'plot' arguments

  #output == object of class 'list', with each compartment containing a slice of the input point cloud

  if(class(XYZtree) != 'data.frame') XYZtree = as.data.frame(XYZtree)
  rg = range(XYZtree[,3])

  if(is.null(l.int)){
    ints = seq(rg[1], rg[2], length.out = n.int+1)
  } else {
    ints = seq(rg[1]-(l.int/2), rg[2]+l.int, by = l.int)
  }

  classes = cut(XYZtree[,3], breaks = ints)
  section.list = split(XYZtree, f = classes)
  section.list = section.list[sapply(section.list,nrow)>0]

  if(!is.null(overlap) && length(section.list)>1){
    add = ints - l.int*overlap
    for(i in 2:length(section.list)){
      extra = section.list[[i-1]]
      extra = extra[extra[,3]>=add[i],]
      section.list[[i]] = rbind(section.list[[i]], extra)
      names(section.list)[i] = paste('(',add[i],',',ints[i+1],']', sep = '')
    }
  }

  if(Plot){
    lapply(section.list, function(u) plot(u[,2] ~ u[,1], pch=20, cex=.5,
                                          main = paste(round(range(u[,3]),digits = 2), units ,collapse = ' - ', sep=' '),
                                          xlab='x', ylab='y', ...))
  }

  return(section.list)
}


###Hough transformation

#' Tree base Hough transformation filter
#' @description identification of reference cylinder at the tree's base for filtering outliers using the Hough transformation
#' @param XYZmat single tree point cloud - \emph{xyz} matrix
#' @param z.int optional - height interval to take the reference cylinder. If not specified, the height interval adopted is from 5\% to 10\% of the tree's total height
#' @param rad.inf inflation factor to multuply the radius. All points (in the entire point cloud) outside a range of \emph{rad.inf * radius} from the reference cylinder's center will be deleted
#' @param cell.size pixel size for the Hough transformation
#' @param min.val passed on to \code{\link{hough}}
#' @param Plot plot the reference tree segment? TRUE or FALSE
#' @return  vector of length 5, containing the upper and lower height limits of the reference cylinder, xy coordinates of circle center and its radius
#' @seealso \code{\link{hough}}
#' @export
HT_base_filter = function(XYZmat, z.int = NULL, rad.inf = 2, cell.size = .025, min.val=.3, Plot=F){

  #Description:

  #XYZmat == matrix or data frame with 3 columns containing x, y and z coordinates, respectively
  #z.int == lower and upper height limits for extracting the reference cylinder
  #rad.inf == inflation factor to apply over the cylinder radius
  #Plot == if TRUE plots the xy coordinates of the extracted section of the input point cloud

  #output == vector of length 5, containing the upper and lower height limits of the reference cylinder,
  #xy coordinates of circle center and its radius

  if(is.null(z.int)) z.int = min(XYZmat[,3]) + c(.05,.1)*(max(XYZmat[,3]) - min(XYZmat[,3])) else z.int = min(XYZmat[,3]) + z.int

  chunk = XYZmat[XYZmat[,3] >= z.int[1] & XYZmat[,3] <= z.int[2],]

  ras = makeRaster(chunk, cell.size = cell.size, image = F)
  rad = hough(ras, pixel_size = cell.size, Plot = F, min.val = min.val)
  top = which(rad$centers[,4] == max(rad$centers[,4]))
  if(nrow(rad$centers[top,,drop=F]) > 1) goal = apply(rad$centers[top,], 2, mean) else goal = rad$centers[top,]

  if(Plot){
    angs = seq(0, 2*pi, length.out = 360)
    plot(chunk[,2] ~ chunk[,1], xlab='x', ylab='y', pch=20, cex=.5, main = paste(round(z.int,2), collapse = ' - '))
    points(x=goal[1], y=goal[2], col='blue', pch=3)
    lines(x = goal[1] + cos(angs)*rad.inf*goal[3], y = goal[2] + sin(angs)*rad.inf*goal[3], lwd=2, col='blue')
  }
  return(c(z = z.int, xy = goal[1:2], radius = rad.inf*goal[3]))
}

#' Circle hough transformation for a single circle
#' @description estimates the circle parameters for a point cloud using the Hough transformation
#' @param raster output from \code{\link{makeRaster}}
#' @param rad lower and upper limits of radii to survey
#' @param pixel_size pixel side length, in meters
#' @param min.val minimum pixel density or frequency that applies for testing
#' @param Plot if TRUE, saves a .png file showing all sample points generated for the input raster
#' @param img.prefix file name for \emph{Plot}
#' @param ... further arguments passed on to \code{\link[graphics]{plot}}
#' @seealso \code{\link{makeRaster}}
#' @return object of class 'list' with 3 compartments: \describe{\item{$centers}{
#' matrix with 4 columns, each row contains the x and y center coordinates,circle radius and number of 'votes' (overlapping peripheral circles)}
#' \item{$images}{list of arrays, each one containing the votes per raster cell for one iterated radius}
#' \item{$circles}{list of xy coordinates of all peripheral circles tested, each compartment contains a matrix with xy point coordinates from circles with same radius (the ones used in 'images')
#' }}
#' @export
hough = function(raster, rad = c(.025,.5), pixel_size = .025, min.val = .1, Plot = F, img.prefix = '', ...){

  rads = seq(rad[1],rad[2], pixel_size)
  angs = seq(0, 2*pi, pixel_size / rad[2])

  nRad = length(rads)
  nAng = length(angs)

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
  centers = matrix(ncol=4,nrow=0)
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

    bullseye = which(base == max(base), arr.ind = T)
    votes = max(base)

    if(nrow(bullseye) == 1){
      x.cen = means.x[bullseye[1]]
      y.cen = means.y[bullseye[2]]
      xy.cen = c(x.cen, y.cen)
      radii = rads[i]
      append = c(xy.cen, radii, votes)
    }else{
      x.cen = means.x[bullseye[,1]]
      y.cen = means.y[bullseye[,2]]
      xy.cen = cbind(x.cen, y.cen)
      votes = rep(votes, nrow(xy.cen))
      radii = rep(rads[i], nrow(xy.cen))
      append = cbind(xy.cen, radii, votes)
    }

    centers = rbind(centers, append)

    split.base = split(main.base, base.index)
    split.base = lapply(split.base, function(u) matrix(u, ncol=ncol(main.base), byrow = F))

    split.circs = split(circs, rad.index)
    split.circs = lapply(split.circs, function(u) matrix(u, ncol=2, byrow=F))
  }
  #main.base = split(main.base, base.index)
  #circs = split(circs, rad.index)
  return(list(centers = centers, images = split.base, circles = split.circs))
}

#' Create raster from a point cloud
#' @description extracts a raster (\emph{xy}) containing density or frequency information from a xyz a point cloud
#' @param XYZmat point cloud - \emph{xyz} matrix
#' @param cell.size pixel size
#' @param density if \code{TRUE}, each pixel will display point density values, if \code{FALSE}, frequency values will be used
#' @param image if \code{TRUE}, plots the raster as an image
#' @return object of clas 'list' with 4 compartments: \describe{\item{$x}{cell breaks in x}
#' \item{$y}{cell breaks in y}
#' \item{$z}{matrix representing each individual raster cell with their respective density or frequency values}
#' \item{$classes}{x and y intervals explicitly identified according to the values in the z compartment}}
#' @export
makeRaster = function(XYZmat, cell.size = .01, density = T, image = T){

  rgX = range(XYZmat[,1])
  rgY = range(XYZmat[,2])

  lx = rgX[2] - rgX[1]
  ly = rgY[2] - rgY[1]

  clx = seq(rgX[1]-(cell.size/2), rgX[2]+cell.size, by=cell.size)
  cly = seq(rgY[1]-(cell.size/2), rgY[2]+cell.size, by=cell.size)

  cutx = cut(x = XYZmat[,1], breaks = clx)
  cuty = cut(x = XYZmat[,2], breaks = cly)

  counts = table(cutx, cuty)
  densities = counts/max(counts)

  if(density){
    raster = matrix(densities, length(clx)-1, length(cly)-1)
  }else{
    raster = matrix(counts, length(clx)-1, length(cly)-1)
  }

  lst = list(x = clx, y = cly, z = as.matrix(raster), classes = if(density) densities else counts)

  if(image) image(lst, col=grey.colors(length(unique(as.vector(raster)))))

  return(lst)

}


###Circle fit

#' Least squares circle fit
#' @description Fits a circle to a set of points - adapted from the \code{pracma} package.
#' @param xp x coordinates
#' @param yp y coordinates
#' @param fast if \code{TRUE} skips the optimization step
#' @param c0 if \code{TRUE}, centers x and y to zero
#' @return vector object containing x, y center coordinates, circle radius and sum of squared errors
#' @export
circlefit = function (xp, yp, fast = FALSE, c0=T){

  if (!is.vector(xp, mode = "numeric") || !is.vector(yp, mode = "numeric"))
    stop("Arguments 'xp' and 'yp' must be numeric vectors.")
  if (length(xp) != length(yp))
    stop("Vectors 'xp' and 'yp' must be of the same length.")

  if(c0){
    cen = c(x=mean(xp), y=mean(yp))
    xp = xp-cen[1]
    yp = yp-cen[2]
  }

  n <- length(xp)
  p <- qr.solve(cbind(xp, yp, 1), matrix(xp^2 + yp^2, ncol = 1))
  r <- c(p[1]/2, p[2]/2, sqrt((p[1]^2 + p[2]^2)/4 + p[3]))
  cerr <- function(v) sqrt(sum((sqrt((xp - v[1])^2 + (yp -
                                                        v[2])^2) - v[3])^2)/n)
  if (fast) {
    cat("RMS error:", cerr(r), "\n")
  }
  else {
    q <- optim(p, cerr)
    #cat("RMS error:", q$value, "\n")
    r <- q$par
  }

  out = unlist(q[1:2])
  if(c0) out[1:2] = out[1:2]+cen

  return(out)
}

#' RAndom SAmple Consensus circle fit
#' @description returns the best circle parameters using the RANSAC algorithm
#' @param stem.sec stem section, \emph{xyz} matrix
#' @param n number of points to sample on every RANSAC iteration
#' @param p estimated proportion of inliers in the dataset
#' @param P level of confidence desired
#' @return vector object containing x, y center coordinates, circle radius and sum of squared errors
#' @export
RANSAC.circle = function(stem.sec, n=15, p=.8, P=.99){

  slc = stem.sec

  if(nrow(stem.sec) < n) n = nrow(stem.sec)

  N = log(1 - P) / log(1 - p^n)

  data = matrix(ncol=4, nrow=0)
  for(j in 1:(5*N)){
    a = sample(1:nrow(slc), size = n)

    b = tryCatch(circlefit(slc[a,1], slc[a,2]),
                 error = function(con){ return('next') },
                 warning = function(con) return('next'))

	if(b == 'next') next

    #if(class(try(circlefit(slc[a,1], slc[a,2]), silent = T)) == "try-error") next
    #b = circlefit(slc[a,1], slc[a,2])
    data = rbind(data, b)
  }

  if(nrow(data) == 0){ dt = NULL }else{

    c = which(data[,4] == min(data[,4]))
    dt = if(length(c) > 1) data[sample(c, size = 1),] else data[c,]

  }

  return(dt)

}


###Spectral decomposition filter

#' Spectral Decomposition point cloud filter
#' @description Removes points from the provided point cloud whose neighborhoods don't follow the specified criteria
#' @param sec tree section, \emph{xyz} matrix
#' @param k number of closest points in a neighborhood over which spectral decomposition is performed
#' @param flat.min minimum flatness accepted to keep points in the dataset
#' @param ang.tol angle tolerance (in degrees) between \code{sec}'s normal vector and \emph{z}
#' @return filtered tree section
#' @export
SD.prefilt = function(sec, k, flat.min, ang.tol){

  if(nrow(sec) < k) k = nrow(sec)

  dists = as.matrix(dist(sec))
  srr = apply(dists, 2, function(u){ a = sort(u, index.return=T)[[2]] ; b = a[1:k] ; return(b) })
  cld = apply(srr, 2, function(u){sec[u,]})

  fln = sapply(cld, FL)
  angs = sapply(cld, ang.deg)

  out = sec[fln > flat.min & (abs(angs-90) < ang.tol | abs(angs-270) < ang.tol),]

  return(out)
}


###Rough noise removal

#' Rough noise removal based on 3D spheres
#' @description removes isolated points from a point cloud
#' @param xyz.tree tree point cloud, \emph{xyz} matrix
#' @param ball.rad spheres radii for first cleaning
#' @param np.min minimum number of points per sphere - spheres with less points will be considered noise and thus removed
#' @param sec.filt apply noise filter a second time?
#' @param lball.rad spheres radii for second cleaning (should be larger than ball.rad) - only used if \code{sec.filt == TRUE}
#' @param min.ncov minimum number of points per sphere for the second cleaning - only used if \code{sec.filt == TRUE}
#' @return tree point cloud without rough noise
#' @export
balls = function(xyz.tree, ball.rad = .025, np.min = 2, sec.filt = T, lball.rad = .05, min.ncov = 3){

  tree.list = Vsections(xyz.tree, l.int = 3*ball.rad, Plot = F)

  #first filtering step
  for(i in 1:length(tree.list)){

    chk = tree.list[[i]]

    xz = chk[,1]
    zx = chk[,3]

    chk[,1] = zx
    chk[,3] = xz

    tmp.ls = Vsections(chk, l.int = ball.rad*3, overlap = 1/3, Plot = F)

    dists = lapply(tmp.ls, function(u){ as.matrix(dist(u)) })

    counts = lapply(dists, function(u) apply(u,1, function(v) length(v[v < ball.rad])) )

    chk = do.call(rbind, tmp.ls)
    counts = unlist(counts)

    chk = chk[counts > np.min,]
    chk = unique(chk)

    zx = chk[,3]
    xz = chk[,1]

    chk[,1] = zx
    chk[,3] = xz

    tree.list[[i]] = chk
  }

  #first tree output
  tree = do.call(rbind,tree.list)
  tree.list = Vsections(tree, l.int = 3*lball.rad, Plot=F)

  #second filtering step
  if(sec.filt){

    for(i in 1:length(tree.list)){

      chk = tree.list[[i]]

      xz = chk[,1]
      zx = chk[,3]

      chk[,1] = zx
      chk[,3] = xz

      tmp.ls = Vsections(chk, l.int = lball.rad*3, overlap = 1/3, Plot = F)

      dists = lapply(tmp.ls, function(u){ as.matrix(dist(u)) })

      counts = lapply(dists, function(u) apply(u,1, function(v) length(v[v < (ball.rad + lball.rad)])) )

      chk = do.call(rbind, tmp.ls)
      counts = unlist(counts)

      chk = chk[counts > min.ncov,]
      chk = unique(chk)

      zx = chk[,3]
      xz = chk[,1]

      chk[,1] = zx
      chk[,3] = xz

      tree.list[[i]] = chk
    }

  }

  #second tree output
  tree = do.call(rbind, tree.list)

  return(tree)
}


###Voxel space and 3D neighborhoods

#' 3D sample points cloud construction
#' @description creates sample points in the 3D space randomly or systematically distributed
#' @param xyz.tree tree point cloud, \emph{xyz} matrix
#' @param l minimum distance between sample points
#' @param systematic distribute the points systematically as a voxel grid? TRUE or FALSE
#' @param n,iterate only used if \code{systematic == FALSE}. The higher those values are, the "denser" is the resuling sample points point cloud
#' @return 3D point cloud of sample points - \emph{xyz} matrix
#' @export
cube.space = function(xyz.tree, l=.03, systematic = T, n=3 , iterate=1){

  #tree cloud processing

  rx = range(xyz.tree[,1])
  ry = range(xyz.tree[,2])
  rz = range(xyz.tree[,3])

  intx = seq(rx[1]-l, rx[2]+l, l)
  inty = seq(ry[1]-l, ry[2]+l, l)
  intz = seq(rz[1]-l, rz[2]+l, l)

  if(systematic){

    #systematic 3D sample points

    xx = (intx[-length(intx)] + intx[-1]) / 2
    yy = (inty[-length(inty)] + inty[-1]) / 2
    zz = (intz[-length(intz)] + intz[-1]) / 2

    dummy = cbind(NA,NA,zz)
    unsplit= xyz.tree#do.call(rbind,split)
    colnames(dummy) = colnames(unsplit)

    comb = rbind(unsplit,dummy)

    split2 = Vsections(comb, l.int = 10*l, Plot = F)
    z.smps = lapply(split2, function(u) u[is.na(u[,1]),3])

    zobs = sapply(z.smps, length)

    split2 = split2[zobs > 0]
    z.smps = z.smps[zobs > 0]

    smps = lapply(1:length(z.smps), function(u){
      a = split2[[u]]
      a = a[!is.na(a[,1]),]

      if(nrow(a) == 0) smps = NULL else{

        cax = cut(a[,1], breaks = intx)
        cay = cut(a[,2], breaks = inty)
        caz = cut(a[,3], breaks = intz)

        cl.a = cbind(cax,cay,caz)
        cl.a = unique(cl.a)

        rx = range(a[,1])
        ry = range(a[,2])

        sx = xx[xx >= rx[1] & xx <= rx[2]]
        sy = yy[yy >= ry[1] & yy <= ry[2]]
        sz = z.smps[[u]]

        xy = expand.grid(sx,sy)
        z = rep(sz, each=nrow(xy))

        smps = cbind(xy,z)

        clsx = cut(smps[,1], breaks = intx)
        clsy = cut(smps[,2], breaks = inty)
        clsz = cut(smps[,3], breaks = intz)

        cl.smps = cbind(clsx,clsy,clsz)

        vec.a = apply(cl.a, 1, paste, collapse=':')
        vec.smps = apply(cl.smps, 1, paste, collapse=':')

        log = vec.smps %in% vec.a

        smps = smps[log,]
      }
      return(smps)

    })

    samples = do.call(rbind,smps)

  }else{

    #random 3D sample points

    rd.gen = function(u, n){
      x = runif(n, intx[u[1]], intx[u[1]+1])
      y = runif(n, inty[u[2]], inty[u[2]+1])
      z = runif(n, intz[u[3]], intz[u[3]+1])
      out = c(x,y,z)
      return(out)
    }
    d.rm = function(chk){
      dst = as.matrix(dist(chk))
      dst[upper.tri(dst, diag = T)] = 100

      a = all(dst >= l)

      while(a == F){
        wch = which(dst < l, arr.ind = T)
        wch = unique(wch[,2])
        wch = wch[1:ceiling(length(wch)/2)]

        chk = chk[-wch,]
        dst = as.matrix(dist(chk))
        dst[upper.tri(dst, diag = T)] = 100

        a = all(dst >= l)
      }

      return(chk)

    }

    Ssamples = matrix(ncol=3,nrow=0)
    for(i in 1:iterate){
      spp = apply(u.ind, 1, rd.gen, n=n)
      nr =nrow(spp)/3
      sppx = c(spp[1:nr,])
      sppy = c(spp[(nr+1):(2*nr),])
      sppz = c(spp[(2*nr+1):(3*nr),])
      spp = cbind(sppx,sppy,sppz)

      spp.list = Vsections(spp, l.int = l*3, overlap = 1/3, Plot = F)
      spp.list = lapply(spp.list, d.rm)

      samples = do.call(rbind,spp.list)
      samples = as.matrix(unique(samples))
      Ssamples = rbind(Ssamples,samples)
    }

    tr.sp = apply(Ssamples, 1, function(u){ x = u[1] + runif(1,l/2,l)
    y = u[2] + runif(1,l/2,l)
    z = u[3] + l
    return(c(x,y,z)) })
    tr.sp2 = apply(Ssamples, 1, function(u){ x = u[1] - runif(1,l/2,l)
    y = u[2] - runif(1,l/2,l)
    z = u[3] - l
    return(c(x,y,z)) })

    Ssamples = rbind(Ssamples, t(tr.sp), t(tr.sp2))

    S.lst = Vsections(Ssamples, l.int = 3*l, overlap = 1/3, Plot = F)
    S.lst = lapply(S.lst, d.rm)
    samples = do.call(rbind, S.lst)
    samples = as.matrix(unique(samples))

  }

  return(samples)

}

#' Cover sets neghborhoods assignment
#' @description Defines cover set spherical neighborhoods for a 3D point cloud based on euclidian distances
#' @param xyz.tree tree point cloud, \emph{xyz} matrix
#' @param samples3d output point cloud from \code{\link{cube.space}} built for \code{xyz.tree}
#' @param d radius of spherical neighborhoods
#' @param neighborhood order of neighborhoods to include - 1 means only the direct neighbours, 2 includes the neighbours of first neighbours and so on
#' @return list object in which every compartment contains a cover set neighborhood
neighbours = function(xyz.tree, samples3d, d=.03, neighborhood=2){

  requireNamespace('foreach', quietly = T)

  xyz.tree = xyz.tree[order(xyz.tree[,3],xyz.tree[,2],xyz.tree[,1]),]
  samples3d = samples3d[order(samples3d[,3],samples3d[,2],samples3d[,1]),]

  #indexation of SAMPLE points
  ids = 1:nrow(samples3d)
  idt = rep(0,nrow(samples3d))
  smp = cbind(samples3d,ids,idt)

  #indexation of TREE points
  idt = 1:nrow(xyz.tree)
  ids = rep(0, nrow(xyz.tree))
  trr = cbind(xyz.tree, ids, idt)

  names(trr) = names(smp)
  #all points in a single dataset
  tr.sp = rbind(smp,trr)

  #split and reorder
  ptlist = Vsections(tr.sp, l.int = 10*d, overlap = 1/3, Plot = F)
  ptord = lapply(ptlist, function(u) u[order(u[,4],-u[,5]),])
  ptord = ptord[sapply(ptord, nrow)!=0]

  #indices assignment
  ptord = lapply(1:length(ptord), function(i){ #foreach(i = 1:length(ptord)) %do% {
    temp = ptord[[i]]
    trp = which(temp[,4]==0)
    spp = which(temp[,5]==0)

    dst = as.matrix(dist(temp[,1:3]))
    dst = dst[trp , spp]

    ngs = which(dst <= d, arr.ind = T)
    if(length(ngs) == 1) ngs = matrix(c(ngs,ngs),ncol=2)

    urw = unique(ngs[,1])
    vw=sapply(urw,function(u){
      a = ngs[ngs[,1] == u,]
      if(!is.matrix(a)) a = matrix(a, ncol=2)
      ln = temp[a[,2]+max(trp),4]
      id = paste(ln,collapse = ':')
      return(id)
    })

    urw2 = unique(ngs[,2])
    vw2=sapply(urw2,function(u){
      a = ngs[ngs[,2] == u,]
      if(!is.matrix(a)) a = matrix(a, ncol=2)
      ln = temp[a[,1],5]
      id = paste(ln,collapse = ':')
      return(id)
    })

    vw.na = rep(NA, nrow(temp))
    vw.na[urw] = vw
    vw2.na = rep(NA, nrow(temp))
    vw2.na[urw2+max(trp)] = vw2

    temp = cbind(temp,vw.na,vw2.na)
    return(temp)
  })

  #rejoining tree and sample matrices
  out = do.call(rbind, ptord)
  out1 = unique(out[!is.na(out[,6]),c(1:3,6)])
  out2 = unique(out[!is.na(out[,7]),c(1:3,7)])

  d.out1 = which(duplicated(out1[,1:3]) | duplicated(out1[,1:3], fromLast = T))
  d.out2 = which(duplicated(out2[,1:3]) | duplicated(out2[,1:3], fromLast = T))


  #sample points containing EACH TREE POINT

  temp1 = out1[d.out1,]
  indx1 = apply(temp1[,1:3], 1, paste, collapse=':')

  cb1 = sapply(indx1, function(u){ a = which(indx1 == u)
  b = paste(temp1[a,4], collapse=':')
  return(b) })
  temp1 = unique(cbind(temp1[,1:3],cb1))
  names(temp1) = names(out1)
  out1 = rbind(out1[-d.out1,], temp1)
  out1 = out1[order(out1[,3],out1[,2],out1[,1]),]

  #tree points comprised by EACH SAMPLE POINT

  temp2 = out2[d.out2,]
  indx2 = apply(temp2[,1:3], 1, paste, collapse=':')

  cb2 = sapply(indx2, function(u){ a = which(indx2 == u)
  b = paste(temp2[a,4], collapse=':')
  return(b) })
  temp2 = unique(cbind(temp2[,1:3],cb2))
  names(temp2) = names(out2)
  out2 = rbind(out2[-d.out2,], temp2)
  out2 = out2[order(out2[,3],out2[,2],out2[,1]),]

  row.names(out1) = NULL
  row.names(out2) = NULL

  ## merging 1st neighbour covers ##

  pts2 = out2[,4]
  pts2 = lapply(pts2, function(u) as.integer(unlist(unique(strsplit(as.character(u), split = ':')))))
  pts2 = lapply(pts2, unique)

  covers = lapply(pts2, function(u) xyz.tree[u,])

  pts1 = out1[,4]
  pts1 = lapply(pts1, function(u) as.integer(unlist(unique(strsplit(as.character(u), split = ':')))))
  pts1 = lapply(pts1, unique)

  cv2 = list()
  foreach(j = 1:length(pts1)) %do% {
    jn = pts1[[j]]
    cv2[jn] = lapply(jn, function(u) unlist(c(cv2[u],jn)))
  }

  cv2 = lapply(cv2, unique)

  foreach(1:neighborhood) %do%{ covers = lapply(cv2, function(u) do.call(rbind, covers[u])) }
  covers = lapply(covers, unique)
  covers = covers[!unlist(lapply(covers, function(u) is.null(u)))]
  covers = lapply(covers, function(u) u[order(u[,3],u[,2],u[,1]),])
  covers = unique(covers)
  covers = lapply(covers, function(u){row.names(u)=NULL;return(u)})

  return(covers)

}

#' Cover sets neghborhoods assignment
#' @description Defines cover set cubic neighborhoods based on a voxel space
#' @param xyz.tree tree point cloud, \emph{xyz} matrix
#' @param samples3d output point cloud from \code{\link{cube.space}} built for \emph{xyz.tree} - the points here must be systematically distributed
#' @param neighborhood order of voxel neighborhoods to include - 1 means only the direct neighbour voxels, 2 includes the neighbours of direct neighbours and so on
#' @return list object in which every compartment contains a cover set neighborhood
#' @export
neighbours2 = function(xyz.tree, samples3d, neighborhood=3){

  xvals = sort(unique(samples3d[,1]))
  yvals = sort(unique(samples3d[,2]))
  zvals = sort(unique(samples3d[,3]))

  d = xvals[2] - xvals[1]

  xvals = c(xvals[1]-d/2, xvals, rev(xvals)[1]+d/2)
  yvals = c(yvals[1]-d/2, yvals, rev(yvals)[1]+d/2)
  zvals = c(zvals[1]-d/2, zvals, rev(zvals)[1]+d/2)

  k = 1:neighborhood

  kn = lapply(k, function(i){

    kx = seq(i, length(xvals) - (neighborhood-i), by = neighborhood) #; if(rev(kx)[1] != length(xvals)) kx = c(kx,length(xvals))
    ky = seq(i, length(yvals) - (neighborhood-i), by = neighborhood) #; if(rev(ky)[1] != length(yvals)) kx = c(ky,length(yvals))
    kz = seq(i, length(zvals) - (neighborhood-i), by = neighborhood) #; if(rev(kz)[1] != length(zvals)) kz = c(kz,length(zvals))

    clx = cut(xyz.tree[,1], breaks = xvals[kx])
    cly = cut(xyz.tree[,2], breaks = yvals[ky])
    clz = cut(xyz.tree[,3], breaks = zvals[kz])

    classes = cbind(clx,cly,clz)
    uni.cls = unique(classes)

    cl.vec = apply(classes, 1, paste, collapse=':')
    uni.vec = apply(uni.cls, 1, paste, collapse=':')

    covs = lapply(uni.vec, function(u) which(cl.vec == u))

    groups = lapply(1:length(covs), function(u){
      a = covs[[u]]
      grp = if(length(a) < 3) NULL else xyz.tree[a,]
      return(grp)
    })

    groups = groups[!sapply(groups, is.null)]

    return(groups)

  })

  covers = do.call(c, kn)

  return(covers)
}


###Cylinder fit

#' Distance from the cylinder surface
#' @description Calculates the distances from every point in a point cloud to its model cylinder
#' @param ang.rad vector containing the 5 cylinder parameters: \emph{rho, theta, phi, alpha and radius} - as in \code{\link{cyl.parameters}} output
#' @param P cylinder point cloud - \emph{xyz} matrix
#' @return vector of distances from every point to the cylinder's surface
#' @export
cyl.dists = function(ang.rad, P){

  rho = ang.rad[1]
  theta = ang.rad[2]
  phi = ang.rad[3]
  alpha = ang.rad[4]
  r = ang.rad[5]

  n = c( cos(phi) * sin(theta) , sin(phi) * sin(theta) , cos(theta) )
  n.theta = c( cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta) )
  n.phi = c(-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0 )
  n.phi.bar = n.phi / sin(theta)

  a = n.theta * cos(alpha) + n.phi.bar * sin(alpha)

  Q = crossprod(rho + r,n)

  dists = apply(P, 1,
                #function(u) (inv.r/2) * ( crossprod(u) - 2*rho*crossprod(u,n) - crossprod(u,a)^2 + rho^2 ) + rho - crossprod(u,n) )
                function(u) sqrt(sum( Xprod((u - Q),a)^2)) - r )

  return(dists)

}

#' Sum of squared distances from a cylinder surface
#' @description Calculates the sum of squared distances from a cylinder surface
#' @param ang.rad vector containing the 5 cylinder parameters: \emph{rho, theta, phi, alpha and radius} - as in \code{\link{cyl.parameters}} output
#' @param P cylinder point cloud - \emph{xyz} matrix
#' @return sum of squared distances from \emph{P} to the cylinder surface with \emph{ang.rad} parameters
#' @export
cyl.fit = function(ang.rad, P){

  #ang.rad = c(rho, theta, phi, alpha, r)

  rho = ang.rad[1]
  theta = ang.rad[2]
  phi = ang.rad[3]
  alpha = ang.rad[4]
  r = ang.rad[5]

  n = c( cos(phi) * sin(theta) , sin(phi) * sin(theta) , cos(theta) )
  n.theta = c( cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta) )
  n.phi = c(-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0 )
  n.phi.bar = n.phi / sin(theta)

  a = n.theta * cos(alpha) + n.phi.bar * sin(alpha)

  Q = crossprod(rho + r,n)

  dists = apply(P, 1,
                #function(u) (inv.r/2) * ( crossprod(u) - 2*rho*crossprod(u,n) - crossprod(u,a)^2 + rho^2 ) + rho - crossprod(u,n) )
                function(u) sqrt(sum( Xprod((u - Q),a)^2)) - r )

  sumsq = sum(dists^2)

  return(sumsq)

}

#' Estimate cylinder parameters from a point cloud
#' @description Estimates the cylinder parameters of a point cloud
#' @param P cylinder point cloud - \emph{xyz} matrix
#' @param init (optional) initial guess of cylinder parameters
#' @param opt.method optimization method - passed on to \code{\link[stats]{optim}}
#' @return \code{\link[stats]{optim}} output containing the optimal cylinder parameters of \emph{P}, respectively: \emph{rho, theta, phi, alpha and radius}
#' @export
cyl.parameters = function(P, init = NULL, opt.method = 'Nelder-Mead'){

  #require(optimx)

  if(is.null(init)){

    center = apply(P, 2, mean)
    center[3] = 0
    rho = sqrt(sum(center^2))

    n = center / rho
    n[3] = 0
    theta = acos(n[3])
    phi = asin(n[2]/sin(theta))

    n.theta = c( cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta) )
    n.phi = c(-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0 )
    n.phi.bar = n.phi / sin(theta)

    a = function(alpha){
      a1 = n.theta * cos(alpha) + n.phi.bar * sin(alpha)
      pr = crossprod(n, a1)
      return(abs(pr))
    }
    alpha = optimize(f = a, interval = c(0,2*pi))[[1]]

    init = c(rho, theta, phi, alpha, 0)
    init = c(rho, pi/2, 0, 0, 0)
  }

  out = optim(fn = cyl.fit, par = init, P=P, method = opt.method)

  return(out)
}

#' Cylinder direction vectors extraction
#' @description Extracts the vectors that define a cylinder's direction based on its five main parameters - as in \code{\link{cyl.parameters}} output
#' @param ang.rad vector containing the 5 cylinder parameters: \emph{rho, theta, phi, alpha and radius} - as in \code{\link{cyl.parameters}}
#' @return list containing the vectors, respectively: \emph{n, a, Q, n^theta and n^phi bar}
#' @export
cyl.vectors = function(ang.rad){

  #ang.rad = c(rho, theta, phi, alpha, r)

  rho = ang.rad[1]
  theta = ang.rad[2]
  phi = ang.rad[3]
  alpha = ang.rad[4]
  r = ang.rad[5]

  n = c( cos(phi) * sin(theta) , sin(phi) * sin(theta) , cos(theta) )
  n.theta = c( cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta) )
  n.phi = c(-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0 )
  n.phi.bar = n.phi / sin(theta)

  a = n.theta * cos(alpha) + n.phi.bar * sin(alpha)

  Q = crossprod(rho + r,n)

  out = list(n = n, a = a, Q=c(Q), n.theta = n.theta, n.phi.bar = n.phi.bar)

  return(out)

}

dist.filt = function(xyz, bprs, inf = .1, abs = NULL){
  dst = cyl.dists(ang.rad = bprs, P = xyz)
  if(is.null(abs)) dst = abs(dst/bprs[5])
  fac = ifelse(is.null(abs), inf, abs)
  keep = xyz[dst < fac,]
  return(keep)
}

#' Iterated Reweighted Least Squares cylinder fit
#' @description Estimates the optimal cylinder parameters of a point cloud using Iterated Reweighted Least Squares
#' @param P cylinder point cloud - \emph{xyz} matrix
#' @param init (optional) initial guess of cylinder parameters
#' @param max.iter maximum number of iterations
#' @param Plot if TRUE plots the x and y cylinder coordinates highlighting the 25\% lowest weight points
#' @param speed.up for large point clouds (> 500 points) - if TRUE, performs optimization over 500 randomly selected points only
#' @param opt.method optimization method - passed on to \code{\link[stats]{optim}}
#' @return list containing the cylinder optimal parameters, its sum of squared distances from the point cloud, number of iterations and point weights, respectively
#' @export
IRLS = function(P, init=NULL, max.iter = 10, Plot=F, speed.up=F, opt.method = 'Nelder-Mead'){

  #require(optimx)

  if(is.null(init)){

    center = apply(P, 2, mean)
    center[3] = 0
    rho = sqrt(sum(center^2))

    n = center / rho
    theta = acos(n[3])
    phi = asin(n[2]/sin(theta))

    n.theta = c( cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta) )
    n.phi = c(-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0 )
    n.phi.bar = n.phi / sin(theta)

    a = function(alpha){
      a1 = n.theta * cos(alpha) + n.phi.bar * sin(alpha)
      pr = crossprod(n, a1)
      return(abs(pr))
    }
    alpha = optimize(f = a, interval = c(0,2*pi))[[1]]

    init = c(rho, theta, phi, alpha, 0)
    init = c(rho, pi/2, 0, 0, 0)
  }

  if(speed.up){ if(nrow(P) > 500) P = P[sample(1:nrow(P), 500, replace = F),]}

  w = rep(1,nrow(P))

  crit = F
  i = 1
  ss = 0
  while(crit == F){
    ss.prev = ss
    opt = optim(par = init, fn = function(u){ a=cyl.dists(u, P = P) ; ssq = sum(w*a^2) ; return(ssq) }, method = opt.method)
    #,control = list(all.methods=T))
    init = opt[[1]]
    ss = opt$value

    dst = cyl.dists(ang.rad = init, P = P)
    wh = tukey.estimator(dst)

    Y = wh$Y
    w = wh$weights #; print(opt)
    i = i+1
    crit = abs(ss - ss.prev) < .001 || i == max.iter
  }

  #print(opt)

  if(Plot){
    qt = quantile(w)
    plot(P[,2]~P[,1], xlab='x', ylab='y', pch=20, col = ifelse(w<qt[2],'red','black'), main=paste(round(range(P[,3]),2),collapse = '-'))
  }

  return(list(pars = init, sqd = ss, iter=i, weights = w))
}

#' RAndom SAmple Consensus cylinder fit
#' @param stem.sec cylinder point cloud - \emph{xyz} matrix
#' @param n number of points sampled in every RANSAC interation
#' @param p proportion of inliers in \emph{stem.sec}
#' @param P level of confidence desired
#' @param timesN inflation factor to multiply the number of RANSAC iterations by
#' @param init initial cylinder parameter estimates, if any
#' @param opt.method optimization method - passed on to \code{\link[stats]{optim}}
#' @return optimal cylinder parameters of \emph{stem.sec} and sum of squared distances, respectively: \emph{rho, theta, phi, alpha, radius, ssq}
#' @export
RANSAC.cylinder = function(stem.sec, n=20, p=.9, P=.99, timesN = 5, init = NULL, opt.method = 'Nelder-Mead'){

  slc = stem.sec

  if(nrow(stem.sec) < n){ dt = NULL }else{

    N = log(1 - P) / log(1 - p^n)

    data = matrix(ncol=6, nrow=0)
    for(j in 1:(timesN*N)){
      a = sample(1:nrow(slc), size = n)
      temp = slc[a,]

      #if(class(try(cyl.parameters(temp), silent = T)) == "try-error") next

      b = cyl.parameters(temp, init = init, opt.method = opt.method)
      bp = unlist(b[1:2])
      data = rbind(data, bp)
    }

    #if(nrow(data) == 0){ dt = NULL }else{

    c = which(data[,6] == min(data[,6]))
    dt = if(length(c) > 1) data[sample(c, size = 1),] else data[c,]

  }

  return(dt)

}


###Diagnostics tools

#' Isolation of cylinder from stem point cloud
#' @description Takes out a stem section in a specified interval and gets its circle or cylinder parameters
#' @param stem stem point cloud - \emph{xyz} matrix
#' @param H cylinder's/circle's mean height
#' @param bound height interval from H, i.e. H +/- bound
#' @param method \emph{"circle"} or \emph{"cylinder"}
#' @param RANSAC if \code{TRUE} estimates using RANSAC, otherwise applies regular fitting
#' @param fast if \code{TRUE}, uses a reduced samples from the point cloud to speed up processing
#' @param resample sample size for \emph{fast}
#' @return vector of circle or cylinder parameters and its residual sum of squares related to the extracted \emph{stem}'s section
#' @export
cyl.iso = function(stem, H=1.3, bound=.5, method='circle', RANSAC=F, fast=T, resample = 500){
  chunk = stem[stem[,3]<H+bound & stem[,3]>H-bound , ]

  if(method == 'circle'){
    params = if(nrow(chunk) < 3) rep(100,4) else{ if(RANSAC) RANSAC.circle(chunk, p = .99) else circlefit(chunk[,1], chunk[,2]) }
    if(is.null(params)) params = rep(100,4)
    names(params) = c('x','y','r','ssq')
  }
  if(method == 'cylinder'){
    if(fast & nrow(chunk) > resample & !RANSAC) chunk = chunk[sample(1:nrow(chunk), size = resample, replace = F),]
    params = if(nrow(chunk) < 3) rep(100,6) else{ if(RANSAC) RANSAC.cylinder(chunk, p = .99, timesN = 5) else unlist(cyl.parameters(chunk)[1:2]) }
    if(is.null(params)) params = rep(100,6)
    names(params) = c('rho', 'theta', 'phi', 'alpha', 'r', 'ssq')
  }
  return(params)
}

#' Performance check of estimated vs. measured diameters
#' @description Extracts precision and accuracy measures and plots estimated vs. measured diameters for a single tree
#' @param tree tree point cloud - \emph{xyz} matrix
#' @param stem stem point cloud - \emph{xyz} matrix
#' @param d measured vector of diameters along the stem
#' @param H vector of heights on which \emph{d} were taken
#' @param model output from \code{\link{taper.mod}}
#' @param rad.max maximum radius tolerated from automatic extraction
#' @param b shape factor for tree taper, passed on to \code{\link{taper.mod}}
#' @param cyl.len cylinder length, interval from which diameter estimates are taken
#' @param plot.dd plot \emph{estimated vs measured} diameters? \code{TRUE} or \emph{FALSE}
#' @param method 'circle' or 'cylinder'
#' @param ... further arguments passed on to \code{\link[graphics]{plot}}
#' @return list containing summarized stem section-wise results and RMSE
#' @export
real.vs.est = function(tree, stem, d, H, model, rad.max=.5, b=2, cyl.len = .5, plot.dd = T, method = 'circle', ...){
  #Sys.time()
  a = t(sapply(H, function(u) cyl.iso(stem, H=u, bound = cyl.len/2, method = method)))
  #Sys.time()

  Theight = max(tree[,3]) - min(tree[,3])
  rads = solid(Theight, H, rad.max, b)

  a[a[,'r'] > rads,] = NA

  b = cbind(a, real.diam = d, height = H)
  b = b[order(b[,'height']),]

  est.rad = predict(object = model, newdata = list(x = b[is.na(b[,1]),'height']))
  b[is.na(b[,1]),'r'] = est.rad

  rmse = sqrt(sum((b[,'real.diam'] - b[,'r']*2)^2)/nrow(b))

  if(plot.dd){

    plot(b[,'r']*2 ~ b[,'real.diam'], xlab = 'Measured D (m)', ylab = 'Estimated D (m)', cex=2, pch=20,
         col = ifelse(is.na(b[,1]), 'red', 'black'), ...)

    lines(-1:10,-1:10, lty=2, lwd=2, col='blue')
  }

  return(list(result = b, D_RMSE = rmse))

}

#' Precision and accuracy of estimated diameters
#' @description extracts RMSE, bias and Pearson's correlation from estimated diameters in relation to their respective field measured values
#' @param rad vector of estimated radii along a tree's stem
#' @param d correspondent diameters to \emph{rad}, measured from the field
#' @return vector containing RMSE, bias and correlation, respectively
#' @export
res.stats = function(rad, d){
  a = sqrt(sum((2*rad - d)^2 / length(rad)))
  b = sum((2*rad - d)) / length(rad)
  c = cor(2*rad, d)
  return(c(RMSE = a, bias = b, r = c))
}

#' Build 3D reference solid
#' @description builds a 3D reference solid wiath a tree's height and arbitrary shape and base diameter
#' @param tree.height a tree's total height
#' @param z.stem height values from which diameters were taken on the stem
#' @param rad.max maximum base radius tolerated from estimations
#' @param b shape factor for tree taper - 0 for cylinder, 1 for paraboloid, 2 for cone, 3 for neiloid
#' @return radii on \emph{z.stem} heights
#' @export
solid = function(tree.height, z.stem, rad.max, b=2){
  sol = sqrt( ((tree.height - z.stem)/tree.height)^b )
  rads = sol*rad.max
  return(rads)
}

#' Tree taper model
#' @description builds a simple model estimating stem section radii from height above ground
#' @param stem stem point cloud - \emph{xyz} matrix OR stem summary table from a fitting algorithm
#' @param tree tree point cloud - \emph{xyz} matrix
#' @param l.int interval of estimation along the stem (cylinder length)
#' @param method 'circle' or 'cylinder'
#' @param rad.max maximum radius tolerated from automatic extraction
#' @param b shape factor for tree taper, passed on to \code{\link{taper.mod}}
#' @return \code{lm} function output
#' @export
taper.mod = function(stem, tree, l.int=.5, method = 'circle', rad.max=.5, b=2){

  tree.height = max(tree[,3]) - min(tree[,3])

  if(ncol(stem) == 3){

    z.stem = seq(min(stem[,3]), max(stem[,3]), by = l.int)+l.int
    a = t(sapply(z.stem, function(u) cyl.iso(stem, H=u, bound = l.int/2, method = method)))

  } else {

    z.stem = rowSums(stem[,1:2])/2
    a = stem

  }

  solid = sqrt( ((tree.height - z.stem)/tree.height)^b )
  rads = solid*rad.max

  log = a[,'r'] > rads

  a[log,] = NA

  y = a[,'r']
  x = z.stem

  mod = lm(y ~x)

  return(mod)

}
