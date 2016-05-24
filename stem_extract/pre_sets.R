##Pre-filtering phase

# Circle hough transformation pre-filter (based on Olofsson et al. 2014)
pref_HT = function(tree, l.int = .5, cell.size = .025, min.den = .3){
  
  ### INSTRUCTIONS ###
  
  # tree = single tree point cloud. Matrix with 3 columns: x, y and z coordinates, respectively
  # l.int = length of segments to split the tree into, from bottom to top (in meters)
  # cell.size = pixel size for creating raster layers (in meters)
  # min.den = minimum point density of pixels to test circles using the hough transformation
  
    # OUTPUT = rough trunk point cloud, matrix with 3 columns, x, y and z, respectively
  
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


# Spectral decomposition pre-filter (based on Liang et. al 2012)
pref_SD = function(tree, k=30, flat.min=.9, ang.tol=10, l.int=.5, freq.ratio = .25){
  
  ### INSTRUCTIONS ###
  
  # tree = single tree point cloud. Matrix with 3 columns: x, y and z coordinates, respectively
  # k = number of points belonging to a surface to be tested for flatness
  # flat.min = minimum tolerated flatness (between 0 and 1)
  # ang.tol = maximum deviation from a perpendicular angle with the z axis (90 +/- ang.tol, in degrees)
  # l.int = length of segments to split the tree into, from bottom to top (in meters) (for this method it only relates to code speed and to recuce RAM usage)
  # freq.ratio = minimum proportion of number of points in a group, in relation to the largest found group, in order to keep it in the output point cloud (between 0 and 1)
  
    # OUTPUT = rough trunk point cloud, matrix with 3 columns, x, y and z, respectively

  
  comps = Vsections(tree, l.int = l.int, overlap = 1/3 ,Plot = F)
  comps = comps[sapply(comps, nrow) > 3]
  comps = lapply(comps, SD.prefilt, k=k, flat.min=flat.min, ang.tol=ang.tol)
  
  retree = do.call(rbind, comps)
  retree = unique(retree)
  colnames(retree) = c('V1','V2','V3')
  retree = retree[with(retree, order(V3,V2,V1)),]
  
  smps = cube.space(retree, l = .1, systematic = T)
  ngs = neighbours2(retree, smps, neighborhood = 3)
  ng.means = t(sapply(ngs, function(u) apply(u, 2, mean)))
  XY.dists = as.matrix(dist(ng.means[,-3]))
  
  max.d = .2
  groups = list()
  pts.check = c()
  for(i in 1:nrow(ng.means)){
    temp = XY.dists[,i]
    pts = which(temp < max.d)
    pts.filt = pts %in% pts.check
    pts = pts[!pts.filt]
    
    groups[[i]] = do.call(rbind, ngs[pts])
    pts.check = c(pts.check, pts)
  }
  
  groups = groups[!sapply(groups, is.null)]
  groups = lapply(groups, unique)
  
  #for(i in 1:length(groups)) rgl.points(groups[[i]], col=sample(rainbow(100), 1))  
  
  group.size = sapply(groups, nrow)
  group.ratio = group.size / max(group.size)
  
  trunk = groups[group.ratio >= freq.ratio]
  trunk = do.call(rbind, trunk)
  trunk = unique(trunk)
  
  return(trunk)
}


# Voxel neighborhoods pre-filter (based on Raumonen et al. 2013)
pref_VN = function(tree, noise1.rad = .05, noise2.rad=.1 , flat.min = .9, ang.tol=10, neighborhood = 4, largest.cov=NULL, axis.dist = .5){
  
  ### INSTRUCTIONS ###
  
  # tree = single tree point cloud. Matrix with 3 columns: x, y and z coordinates, respectively
  # noise1.rad = sphere radii to use for the first rough noise filtering step (in meters)
  # noise2.rad = sphere radii to use for the second rough noise filtering step (in meters)
  # flat.min = minimum tolerated flatness (between 0 and 1)
  # ang.tol = maximum deviation from a perpendicular angle with the z axis (90 +/- ang.tol, in degrees)
  # neighborhood = order of voxel neighborhood to merge (for voxels of 5 cm)
  # largest.cov = minimum accepeted proportion of points in a cover set, in relation to the largest cover set, to keep it in the output cloud. If =NULL an automated maximum curve detection method is applied.
  # axis.dist = maximum distance from an estimated z axis tolerated to keep a cover set in the output cloud
  
    # OUTPUT = rough trunk point cloud, matrix with 3 columns, x, y and z, respectively
  
  tree2 = balls(tree, ball.rad = noise1.rad, sec.filt = noise2.rad)
  sps = cube.space(tree2)
  
  ngass = neighbours2(tree2, sps, neighborhood = 3)
  
  flness = sapply(ngass, FL)
  norms = sapply(ngass, ang.deg)
  pars = abs(norms - 90)
  
  #flat.min = .95
  #ang.tol = 10
  
  log = flness >= flat.min & pars < ang.tol
  trunk.list = ngass[log]
  
  trunk = do.call(rbind, trunk.list)
  trunk = unique(trunk)
  
  smp.trunk = cube.space(trunk, l = .05)
  retrunk = neighbours2(trunk, smp.trunk, neighborhood)
  
  obs = sapply(retrunk, nrow)
  obs.rat = obs/max(obs)
  
  index = sort(obs, decreasing = T, index.return=T)[[2]]
  i=1
  crit = F
  while(crit == F){
    chunk = do.call(rbind, retrunk[index[1:i]])
    cen = circlefit(chunk[,1], chunk[,2])
    
    i = i+1
    crit = cen[3] >= .05 && cen[3] <= .5 && cen[4] < 1 || i == 10
  }
  
  cen = cen[1:2]
  
  if(is.null(largest.cov)){
    x = length(obs):1
    y = sort(obs.rat)
    
    a = lm(rev(range(y))~range(x))
    
    ang = -atan(a$coefficients[2])
    
    rtm = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2, byrow = T)
    rot = cbind(x,y) %*% rtm
    
    curv = which(rot[,2] == min(rot[,2]))
    
    largest.cov = y[curv]
  }
  
  trunk2 = retrunk[obs.rat > largest.cov]
  
  meandis = sapply(trunk2, function(u){
    xy = apply(u[,1:2], 2, mean)
    euc = sqrt(sum((xy-cen)^2))
    return(euc)
  })
  
  if(i < 10) trunk2 = trunk2[meandis < axis.dist]
  trunk2 = do.call(rbind,trunk2)
  trunk2 = unique(trunk2)
  
  return(trunk2)
}


###Fitting phase

# RANSAC circle fit
fit_RANSAC_circle = function(trunk, l.int = .5, cut.rad = .01, n=15, p=.8, P=.99){

  ### INSTRUCTIONS ###
  
  # trunk = output point cloud from a pre-filtering method. Matrix with 3 columns: x, y and z coordinates, respectively
  # l.int = length of stem segments to fit (in meters)
  # cut.rad = value to add to the radius, removing all points outside a range of radius+cut.rad from the point cloud (in meters)
  # n = number of points to take in every RANSAC iteration
  # p = estimated proportion of inliers (stem points)
  # P = confidence level
  
    # OUTPUT = single stem point cloud, matrix with 3 columns, x, y and z, respectively

  
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


# Iterated reqeighted total least squares cylinder fit
fit_IRTLS = function(trunk, c.len = .5, max.rad=.5, s.height = 1, speed.up = T, opt.method = 'Nelder-Mead'){
  
  ### INSTRUCTIONS ###
  
  # trunk = output point cloud from a pre-filtering method. Matrix with 3 columns: x, y and z coordinates, respectively
  # c.len = length of stem segments to fit (in meters)
  # max.rad = maximum radius accepted as estimate (in meters)
  # s.height = starting height for acquiring baseline cylinder parameters (in meters)
  # speed.up = TRUE or FALSE, if TRUE takes no more than 500 points to fit a cylinder, randomly selected (applicable for large point clouds)
  # opt.method = optimization method, passed to optim (R base function)
  
    # OUTPUT = single stem point cloud, matrix with 3 columns, x, y and z, respectively

  
  zbound = c((s.height-c.len/2) , (s.height+c.len/2)) + min(trunk[,3])
  
  #max.rad = .5
  rrad = 100
  thr=.05
  
  while(max.rad <= rrad){
    
    base = trunk[trunk[,3]>zbound[1] & trunk[,3]<zbound[2],]
    if(nrow(base) < 100*c.len){ zbound = zbound + c.len ; next }
    #plot(base[,-3])
    rreg = IRLS(base, speed.up = speed.up, opt.method = opt.method)
    
    
    rrad = rreg[[1]][5]
    zbound = zbound + c.len
  }
  
  below = trunk[trunk[,3]<zbound[1],]
  u.dist = cyl.dists(rreg$pars, below)
  below = below[abs(u.dist) < thr,]
  
  above = trunk[trunk[,3]>=zbound[1],]
  ab.list = Vsections(above, l.int = c.len, Plot = F)
  
  par = rreg$pars
  for(i in 1:length(ab.list)){
    temp = ab.list[[i]]
    a.dist = cyl.dists(par, temp)
    temp = temp[abs(a.dist)<thr,]
    
    if(nrow(temp) < 20){ab.list[[i]] = temp ; next}
    
    a.rreg = IRLS(temp, init = par, speed.up = speed.up, opt.method = opt.method)
    if(a.rreg[[1]][5] > max.rad) next
    par = a.rreg$pars
    
    a.dist = cyl.dists(par, temp)
    temp = temp[abs(a.dist)<.01,]
    
    ab.list[[i]] = temp
  }
  above = do.call(rbind,ab.list)
  
  stem = rbind(below,above)
  
  return(stem)
}


# RANSAC cylinder fit
fit_RANSAC_cylinder = function(trunk, c.len = .5, h.init = 1, max.rad = .5, timesN = 2, opt.method = 'Nelder-Mead'){
  
  ### INSTRUCTIONS ###
  
  # trunk = output point cloud from a pre-filtering method. Matrix with 3 columns: x, y and z coordinates, respectively
  # c.len = length of stem segments to fit (in meters)
  # max.rad = maximum radius accepted as estimate (in meters)
  # timesN = factor that multiplies N, for the number of iterations of the RANSAC, N=36 as default
  # opt.method = optimization method, passed to optim (R base function)
  
    # OUTPUT = single stem point cloud, matrix with 3 columns, x, y and z, respectively

  
  zbound = c((h.init-c.len/2) , (h.init+c.len/2)) + min(trunk[,3])
  
  rrad = 100
  thr=.05
  
  while(max.rad <= rrad){
    
    base = trunk[trunk[,3]>zbound[1] & trunk[,3]<zbound[2],]
    if(nrow(base) < 100*c.len){ zbound = zbound + c.len ; next }
    #plot(base[,-3])
    rreg = RANSAC.cylinder(base, opt.method = opt.method)
    
    rrad = rreg[5]
    zbound = zbound + c.len
  }
  
  below = trunk[trunk[,3]<zbound[1],]
  u.dist = cyl.dists(rreg[-6], below)
  below = below[abs(u.dist) < thr,]
  
  above = trunk[trunk[,3]>=zbound[1],]
  ab.list = Vsections(above, l.int = c.len, Plot = F)
  
  par = rreg[-6]
  for(i in 1:length(ab.list)){
    #print(Sys.time())
    temp = ab.list[[i]]
    a.dist = cyl.dists(par, temp)
    temp = temp[abs(a.dist)<thr,]
    
    if(nrow(temp) < 20){ab.list[[i]] = temp ; next}
    
    a.rreg = RANSAC.cylinder(temp, timesN = timesN, init = par, p = .95, opt.method = opt.method)
    if(a.rreg[5] > max.rad || a.rreg[5] < .025) next
    par = a.rreg[-6]
    
    a.dist = cyl.dists(par, temp)
    temp = temp[abs(a.dist)<.01,]
    
    #plot(temp[,-3])
    
    ab.list[[i]] = temp
  }
  
  above = do.call(rbind,ab.list)
  stem = rbind(below,above)
  
  return(stem)
}
