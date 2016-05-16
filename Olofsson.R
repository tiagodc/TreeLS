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

circlefit = function (xp, yp, fast = FALSE){
  
  #Description: adapted from the 'pracma' package. Fits a circle to xy coordinates
  
  #xp == x points coordinates
  #yp == y points coordinates
  #fast == if TRUE doesn't perform optimization
  
  #output == x, y, radius and sum of squared errors
  
  
  if (!is.vector(xp, mode = "numeric") || !is.vector(yp, mode = "numeric")) 
    stop("Arguments 'xp' and 'yp' must be numeric vectors.")
  if (length(xp) != length(yp)) 
    stop("Vectors 'xp' and 'yp' must be of the same length.")
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
  return(unlist(q[1:2]))
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

#presets

Olofsson.pref = function(tree, l.int = .5, cell.size = .025, min.den = .3){
  
  baseline = filter(tree, Plot = F, cell.size = cell.size, rad.inf = 2, min.val = min.den)
  dists = apply(tree, 1, function(u) sum( (u[1:2] - baseline[3:4])^2 )^(1/2) )
  noise = which(dists > baseline[5] & tree[,3] < baseline[2])
  tree.cloud = if(length(noise) == 0) tree else tree[-noise,]
  
  slices = Vsections(tree.cloud, l.int = l.int, Plot = F)
  tops = sapply(slices, function(u) max(u[,3]))
  
  take = which(tops <= baseline[2])
  keep = which(tops > baseline[2])
  
  angs = seq(0,2*pi,length.out = 180+1) 
  out = data.frame()
  
  rd.bff = .05
  
  for(i in take){
    ras = makeRaster(slices[[i]], cell.size = cell.size, image = F)
    ctr = hough(ras, nRad = 50, rad = c(.025,baseline[5]*.75), nAng = 120, min.val = min.den, Plot = F)
    row = which(ctr[[1]][,4] == max(ctr[[1]][,4]))
    
    if(length(row) > 1) sec.info = ctr[[1]][sample(row,1),] else sec.info = ctr[[1]][row,]
    
    dis = apply(slices[[i]], 1, function(u) sqrt(sum((u[-3]-sec.info[1:2])^2)) )
    slices[[i]] = slices[[i]][dis<sec.info[3]+rd.bff,]
  }
  
  for(i in keep){
    
    if(i == keep[1]){
      
      dists = apply(slices[[i]], 1, function(u) sum( (u[1:2] - baseline[3:4])^2 )^(1/2) )
      noise = which(dists > baseline[5])
      
    } else {
      
      dists = apply(slices[[i]], 1, function(u) sum( (u[1:2] - sec.info[1:2])^2 )^(1/2) )
      noise = which(dists > sec.info[3]+rd.bff)
      
    }
    
    if(length(noise) > 0) slices[[i]] = slices[[i]][-noise,]
    if(nrow(slices[[i]]) == 0){ next }
    
    ras = makeRaster(slices[[i]], cell.size = cell.size, image = F)
    if(length(ras$z) <= 1) next
    ctr = hough(ras, rad = c(cell.size, (4/3)*sec.info[3]), nRad = 30 ,nAng = 120, min.val = min.den, Plot = F)
    row = which(ctr[[1]][,4] == max(ctr[[1]][,4]))
    
    if(length(row) > 1) sec.info = ctr[[1]][sample(row, 1),] else sec.info = ctr[[1]][row,]
    
    dists = apply(slices[[i]], 1, function(u) sum( (u[1:2] - sec.info[1:2])^2 )^(1/2) )
    noise = which(dists > sec.info[3]+rd.bff/2)
    if(length(noise) > 0) slices[[i]] = slices[[i]][-noise,]
    
  }
  
  slices = slices[sapply(slices,nrow) > 0]
  unsliced = do.call(rbind, slices)
  
  return(unsliced)
  
}

Olofsson.stem = function(trunk, l.int = .5, cut.rad = .01, n=15, p=.8, P=.99){
  
  slices = Vsections(trunk, l.int = l.int, Plot = F)
  
  N = log(1 - P) / log(1 - p^n)
  
  slices2 = list()
  best = matrix(ncol=4,nrow=0)
  dt=NULL
  for(i in 1:length(slices)){
    
    slc = slices[[i]]
    
    if(nrow(slc) < n){ best = rbind(best, rep(0,4)) ; next}
    
    data = matrix(ncol=4, nrow=0)
    for(j in 1:(5*N)){
      a = sample(1:nrow(slc), size = n)
      
      if(class(try(circlefit(slc[a,1], slc[a,2]), silent = T)) == "try-error") next
      
      b = circlefit(slc[a,1], slc[a,2])
      data = rbind(data, b)
    }
    
    if(nrow(data) == 0) next
    
    c = which(data[,4] == min(data[,4]))
    if(length(c) > 1) dt = data[sample(c, size = 1),] else dt = data[c,]
    
    dst = apply(slc, 1, function(u) sqrt( sum( (u[1:2] - dt[1:2])^2)) )
    slices2[[i]] = slc[dst < cut.rad+dt[3],]
    
    best = rbind(best, dt)
    
  }
  
  colnames(best) = c('x','y','radius','ssq')
  best = data.frame(best, h=names(slices), row.names = NULL)
  stem = do.call(rbind, slices2)
  
  return(stem)
}