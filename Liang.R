#functions

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
  a = Vsections(XYZtree, l.int = l.int, Plot = F)
  a = lapply(a , function(u){ if(nrow(u)>thr) u = u[sample(1:nrow(u), size = thr),] ; return(u) })
  a = do.call(rbind, a)
  return(a)
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
  
  #all points in a single dataset
  tr.sp = rbind(smp,as.matrix(trr))
  
  #split and reorder
  ptlist = Vsections(tr.sp, l.int = 10*d, overlap = 1/3, Plot = F)
  ptord = lapply(ptlist, function(u) u[order(u[,4],-u[,5]),])
  ptord = ptord[sapply(ptord, nrow)!=0]
  
  #indices assignment
  foreach(i = 1:length(ptord)) %do% {
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
    
    ptord[[i]] = cbind(temp,vw.na,vw2.na)
    
  }
  
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

circlefit = function (xp, yp, fast = FALSE){
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

change.coords = function(xyz.points, rot.mat, shift=c(0,0,0)){
  xyz.points = as.matrix(xyz.points)
  rot = xyz.points %*% rot.mat
  mrot = t(apply(rot, 1, function(u) u + shift))
  return(mrot)  
}

rglAXES = function(xyz = c(1,1,1), cols = c('red','green','blue'), ...){
  rgl.lines(c(0,xyz[1]), c(0,0), c(0,0), col=cols[1], ...)
  rgl.lines(c(0,0), c(0,xyz[2]), c(0,0), col=cols[2], ...)
  rgl.lines(c(0,0), c(0,0), c(0,xyz[3]), col=cols[3], ...)
}

angle = function(a,b){
  prod = a %*% b
  lprod = sqrt(sum(a^2)) * sqrt(sum(b^2))
  ang = prod/lprod
  cang = acos(ang) * 180/pi
  return(cang[,,drop=T])
}

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

Xprod = function(a,b){
  x = c(a[2]*b[3] - a[3]*b[2] ,
        a[3]*b[1] - a[1]*b[3] ,
        a[1]*b[2] - a[2]*b[1])
  return(x)
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

dist.filt = function(xyz, bprs, inf = .1, abs = NULL){
  dst = cyl.dists(ang.rad = bprs, P = xyz)
  if(is.null(abs)) dst = abs(dst/bprs[5])
  fac = ifelse(is.null(abs), inf, abs)
  keep = xyz[dst < fac,]
  return(keep)
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

#presets

Liang.pref = function(tree, k=30, flat.min=.9, ang.tol=10, l.int=.5, freq.ratio = .25){

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

Liang.stem = function(trunk, c.len = .5, max.rad=.5, s.height = 1, speed.up = T, opt.method = 'Nelder-Mead'){

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