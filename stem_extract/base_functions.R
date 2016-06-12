###Data simulation and analysis

ang.deg = function(XYZplane){
  e = eigen(cov(XYZplane))
  ang = ( e$vectors[,3] %*% c(0,0,1) ) / ( sqrt(sum(e$vectors[,3]^2)) * sqrt(sum(c(0,0,1)^2)) )
  ang = ang[,,drop=T]
  degs = acos(ang)*180/pi
  return(degs)
}

angle = function(a,b){
  prod = a %*% b
  lprod = sqrt(sum(a^2)) * sqrt(sum(b^2))
  ang = prod/lprod
  cang = acos(ang) * 180/pi
  return(cang[,,drop=T])
}

change.coords = function(xyz.points, rot.mat, shift=c(0,0,0)){
  xyz.points = as.matrix(xyz.points)
  rot = xyz.points %*% rot.mat
  mrot = t(apply(rot, 1, function(u) u + shift))
  return(mrot)  
}

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

FL = function(XYZplane){
  e = eigen(cov(XYZplane))
  flat = 1 - ( e$values[3] / sum(e$values) )
  return(flat)
}

pizza = function(cylinderXYZ, n.slices = 8){
  cX = mean(range(cylinderXYZ[,1]))
  cY = mean(range(cylinderXYZ[,2]))
  
  center = c(cX,cY)
  ang.int = seq(0, 2*pi, length.out = n.slices+1) - pi/n.slices
  
  plot(cylinderXYZ[,2] ~ cylinderXYZ[,1], xlab='x', ylab='y')
  apply(cbind(cX, tan(ang.int)), 1, function(x) abline(a=x[1], b=x[2], col='red'))
  abline(h = cX, lty = 3)
  abline(v = cY, lty = 3)
  
  angles = atan( (cylinderXYZ[,2] - cY) / (cylinderXYZ[,1] - cX) )
  upCyl = cbind(cylinderXYZ, angles)
  
  upCyl[(upCyl[,1] < cX),4] = upCyl[(upCyl[,1] < cX),4] + pi
  upCyl[(upCyl[,1] > cX & upCyl[,2] < cY),4] = upCyl[(upCyl[,1] > cX & upCyl[,2] < cY),4] + 2*pi
  
  fac = cut(upCyl[,4], breaks = c(0,ang.int[-1],pi*2))
  
  slices = split(x = as.data.frame(upCyl), f = fac)
  first = rbind(slices[[1]], slices[[length(slices)]])
  slices[[1]] = first
  slices = slices[-length(slices)]
  
  return(slices)
}

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

Xprod = function(a,b){
  x = c(a[2]*b[3] - a[3]*b[2] ,
        a[3]*b[1] - a[1]*b[3] ,
        a[1]*b[2] - a[2]*b[1])
  return(x)
}

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

center.zero = function(xyz.cloud){
  mn = colMeans(xyz.cloud)
  names(mn) = c('x','y','z')
  xyz.cloud[,1] = xyz.cloud[,1] - mn[1]
  xyz.cloud[,2] = xyz.cloud[,2] - mn[2]
  return(list(cloud = xyz.cloud, center = mn[1:2]))
}

clip.XY = function(cloud, rad = 1.5, center = c(0,0)){
  
  dists = apply(cloud, 1, function(u) sqrt(sum((u[1:2]-center)^2)))
  out = cloud[dists <= rad,]
  return(out)
  
}

rglAXES = function(xyz = c(1,1,1), cols = c('red','green','blue'), ...){
  rgl.lines(c(0,xyz[1]), c(0,0), c(0,0), col=cols[1], ...)
  rgl.lines(c(0,0), c(0,xyz[2]), c(0,0), col=cols[2], ...)
  rgl.lines(c(0,0), c(0,0), c(0,xyz[3]), col=cols[3], ...)
}

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

filter = function(XYZmat, z.int = NULL, rad.inf = 2, cell.size = .025, min.val=.3, Plot=F){
  
  #Description: identification of reference cylinder for filtering outliers 
  
  #XYZmat == matrix or data frame with 3 columns containing x, y and z coordinates, respectively
  #z.int == lower and upper height limits for extracting the reference cylinder
  #rad.inf == inflation factor to apply over the cylinder radius
  #Plot == if TRUE plots the xy coordinates of the extracted section of the input point cloud
  
  #output == vector of length 5, containing the upper and lower height limits of the reference cylinder,
  #xy coordinates of circle center and its radius
  
  if(is.null(z.int)) z.int = min(XYZmat[,3]) + c(.05,.1)*(max(XYZmat[,3]) - min(XYZmat[,3])) else z.int = min(XYZmat[,3]) + z.int
  
  chunk = XYZmat[XYZmat[,3] >= z.int[1] & XYZmat[,3] <= z.int[2],]
  
  ras = makeRaster(chunk, cell.size = cell.size, image = F)
  rad = hough(ras, nAng=120, nRad = 50, Plot = F, min.val = min.val)
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

hough = function(raster, rad = c(.025,.5), nRad = 100, nAng = 360, min.val = .1, Plot = F, img.prefix = '', ...){
  
  #Description: hough transformation for circles in xyz point clouds
  
  #raster == output from makeRaster
  #rad == lower and upper limits of radiuses to test
  #nRad == how many radiuses to test within the interval specified by 'rad'
  #nAng == number of sample points for each tested circle
  #min.val == minimum density or frequency value to test for centers
  #Plot == if TRUE saves a .png file showing all sample points generated for the input raster
  #img.prefix == file name for 'Plot'
  #... == further plot arguments
  
  #output == object of class 'list' with 3 compartments:
  #centers == matrix with 4 columns, each row contains the x and y center coordinates, 
  #circle radius and number of 'votes' (overlapping peripheral circles)
  #images == list of arrays, each one containing the votes per raster cell for one iterated radius
  #circles == list of xy coordinates of all peripheral circles tested,
  #each compartment contains a matrix with xy point coordinates from circles with same radius
  #(the ones used in 'images')
  
  
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

makeRaster = function(XYZmat, cell.size = .01, density = T, image = T){
  
  #Description: extracts a raster of xy density or frequency from a xyz a point cloud
  
  #XYZmat == matrix or data frame with 3 columns containing x, y and z coordinates, respectively
  #cell.size == size of the raster cells (side length)
  #density == if TRUE each raster cell registers point density information, if FALSE frequency is recorded
  #image == if TRUE outputs a plot of the generated raster
  
  #output == object of clas 'list' with 4 compartments:
  #x == cell breaks in x
  #y == cell breaks in y
  #z == matrix representing each individual raster cell with their respective density or frequency values
  #classes == x and y intervals explicitly identified according to the values in the z compartment
  
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

circlefit = function (xp, yp, fast = FALSE, c0=T){
  
  #Description: adapted from the 'pracma' package. Fits a circle to xy coordinates
  
  #xp == x points coordinates
  #yp == y points coordinates
  #fast == if TRUE doesn't perform optimization
  
  #output == x, y, radius and sum of squared errors
  
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

RANSAC.circle = function(stem.sec, n=15, p=.8, P=.99){
  
  slc = stem.sec
  
  if(nrow(stem.sec) < n) n = nrow(stem.sec)
  
  N = log(1 - P) / log(1 - p^n)
  
  data = matrix(ncol=4, nrow=0)
  for(j in 1:(5*N)){
    a = sample(1:nrow(slc), size = n)
    
    if(class(try(circlefit(slc[a,1], slc[a,2]), silent = T)) == "try-error") next
    
    b = circlefit(slc[a,1], slc[a,2])
    data = rbind(data, b)
  }
  
  if(nrow(data) == 0){ dt = NULL }else{
    
    c = which(data[,4] == min(data[,4]))
    dt = if(length(c) > 1) data[sample(c, size = 1),] else data[c,]
    
  }
  
  return(dt)
  
}


###Spectral decomposition filter

SD.prefilt = function(sec, k, flat.min, ang.tol){
  
  if(nrow(sec) < k) k = nrow(sec)
  
  dists = as.matrix(dist(sec))
  srr = apply(dists, 2, function(u){ a = sort(u, index.return=T)[[2]] ; b = a[1:k] ; return(b) })
  cld = apply(srr, 2, function(u){sec[u,]})
  
  fln = sapply(cld, FL)
  angs = sapply(cld, ang.deg)
  
  out = sec[fln > flat.min & abs(angs-90) < ang.tol,]
  
  return(out)
}


###Rough noise removal

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
    
    rd.gen = function(u, n=3){
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

neighbours = function(xyz.tree, samples3d, d=.03, neighborhood=2){
  
  require(foreach, quietly = T)
  
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

dist.filt = function(xyz, bprs, inf = .1, abs = NULL){
  dst = cyl.dists(ang.rad = bprs, P = xyz)
  if(is.null(abs)) dst = abs(dst/bprs[5])
  fac = ifelse(is.null(abs), inf, abs)
  keep = xyz[dst < fac,]
  return(keep)
}

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

cyl.iso = function(stem, H=1.3, bound=.5, method='circle'){
  chunk = stem[stem[,3]<H+bound & stem[,3]>H-bound , ]
  
  if(method == 'circle'){
    params = if(nrow(chunk) < 3) rep(100,4) else RANSAC.circle(chunk, p = .99)
    if(is.null(params)) params = rep(100,4)
    names(params) = c('x','y','r','ssq')
  }
  if(method == 'cylinder'){
    params = if(nrow(chunk) < 3) rep(100,6) else RANSAC.cylinder(chunk, p = .99, timesN = 5)
    if(is.null(params)) params = rep(100,6)
    names(params) = c('rho', 'theta', 'phi', 'alpha', 'r', 'ssq')
  }
  return(params)
}

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

resid.eval = function(rad, d){
  a = sqrt(sum((2*rad - d)^2 / length(rad)))
  b = sum((2*rad - d)) / length(rad)
  c = cor(2*rad, d)
  return(c(RMSE = a, bias = b, r = c))
}

solid = function(tree.height, z.stem, rad.max, b=2){
  sol = sqrt( ((tree.height - z.stem)/tree.height)^b )
  rads = sol*rad.max
  return(rads)
}

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
