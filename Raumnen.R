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

rglAXES = function(xyz = c(1,1,1), cols = c('red','green','blue')){
  rgl.lines(c(0,xyz[1]), c(0,0), c(0,0), col=cols[1])
  rgl.lines(c(0,0), c(0,xyz[2]), c(0,0), col=cols[2])
  rgl.lines(c(0,0), c(0,0), c(0,xyz[3]), col=cols[3])
}

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

FL = function(XYZplane){
  e = eigen(cov(XYZplane))
  flat = 1 - ( e$values[3] / sum(e$values) )
  return(flat)
}

ang.deg = function(XYZplane){
  e = eigen(cov(XYZplane))
  ang = ( e$vectors[,3] %*% c(0,0,1) ) / ( sqrt(sum(e$vectors[,3]^2)) * sqrt(sum(c(0,0,1)^2)) )
  ang = ang[,,drop=T]
  degs = acos(ang)*180/pi
  return(degs)
}

Xprod = function(a,b){
  x = c(a[2]*b[3] - a[3]*b[2] ,
        a[3]*b[1] - a[1]*b[3] ,
        a[1]*b[2] - a[2]*b[1])
  return(x)
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

#presets

Raumonen.trunk = function(tree, noise1.rad = .05, noise2.rad=.1 , flat.min = .9, ang.tol=10, neighborhood = 4, largest.cov=NULL, axis.dist = .5){

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

func = lm(log(y)~log(x)-1) #nls(y~x^i), start = list(i=1,j=1))
cf = coef(func)

nlfunc = nls(y~x^a, start = list(a=cf))
nlcf = coef(nlfunc)

dk = function(x, n = nlcf){#a=nlcf[2], b=nlcf[1]){
  
  #dk=abs(n-1)*abs(n)*abs(x^(n-2))/(n^2*x^(2*(n-1))+1)^(3/2)
  dk=-abs(n-1)*abs(n)*x^(2*n-7)*((2*n^3-n^2)*x^(2*n)+(2-n)*x^2)/((n^2*x^(2*(n-1))+1)^(5/2)*abs(x^(n-2)))
  
  return(abs(dk))
}
curvmax = nlm(dk, p = 1)[[2]]
largest.cov = curvmax^nlcf

#png('curve.png', width = 15, height = 15, units = 'cm', res = 300)
#plot(y~x, pch=20, xlab = 'Order of sorted relative point frequency', ylab = 'Relative point frequency')
#lines(x, predict(nlfunc), col='red', lwd=2)
#points(curvmax, largest.cov, col='red', cex=2, pch=18)
#lines(c(-100000, curvmax, curvmax), c(largest.cov, largest.cov, -1), col='red', lty=2)
#legend('topright', pch=c(20, 18, -1), pt.cex=1, lty=c(0,0,1), lwd=2, col=c('black', 'red', 'red'), 
#       legend=c('cover set point frequency', 'estimated point of maximum curvature', expression(paste('rough power function ', y==x^beta))))
#dev.off()

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

Raumonen.stem = function(trunk, c.len = .5, h.init = 1, max.rad = .5, timesN = 2, opt.method = 'Nelder-Mead'){
  
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