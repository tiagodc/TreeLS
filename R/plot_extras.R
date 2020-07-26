xprod = function(a, b){
  x = c(
    a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1]
  )
  return(x);
}

tlsPlot.dh = function(las, pars, clear=T, wired=T, col='white'){
  if('ax' %in% names(pars)){
    tlsPlot.dh.bf.cylinder(las, pars$x, pars$y, pars$radius, pars$ax, pars$ay, clear, wired, col)
  }else if('rho' %in% names(pars)){
    tlsPlot.dh.cylinder(las, pars$rho, pars$theta, pars$phi, pars$alpha, pars$radius, clear, wired, col)
  }else{
    tlsPlot.dh.circle(las, pars$X, pars$Y, pars$radius, clear, wired, col)
  }
}

tlsPlot.dh.bf.cylinder = function(las, x, y, r, ax, ay, clear=F, wired=T, col='white'){
  pt3d = las@data[,.(X,Y,Z)]

  rmat = rotationMatrix(ax*pi/180, ay*pi/180, 0)
  temp = as.matrix(pt3d) %*% rmat %>% as.data.table
  names(temp) = names(pt3d)

  rmat = rotationMatrix(-ax*pi/180, 0, 0) %*% rotationMatrix(0, -ay*pi/180, 0)
  cbase = c(x,y,min(temp$Z)) %*% rmat %>% as.double
  ctop  = c(x,y,max(temp$Z)) %*% rmat %>% as.double

  ccen  = c(x,y,range(temp$Z) %>% mean) %*% rmat %>% as.double
  crad  = c(x+r,y,range(temp$Z) %>% mean) %*% rmat %>% as.double

  tlsPlot.dh.3d(las, cbind(cbase, ctop), cbind(ccen, crad), r, clear, wired, col)

  ptring = cbind(
    X = x + cos(seq(0,2*pi,length.out = 36)) * r,
    Y = y + sin(seq(0,2*pi,length.out = 36)) * r,
    Z = range(temp$Z) %>% mean
  ) %*% rmat

  if(!is.null(las)) tlsPlot.dh.2d(las, ptring, r)
}

tlsPlot.dh.cylinder = function(las, rho, theta, phi, alpha, r, clear=F, wired=T, col='white'){
  n = c(cos(phi) * sin(theta) , sin(phi) * sin(theta) , cos(theta))
  ntheta = c(cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta))
  nphi = c(-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0);
  nphibar = nphi/sin(theta)

  a = ntheta * cos(alpha) + nphibar * sin(alpha)
  q = n*(r+rho)

  meds = apply(las@data[,.(X,Y,Z)], 2, median) %>% as.double

  pt3d = las@data[,.(X,Y,Z)]
  height = las$Z %>% range %>% diff %>% abs
  height = height/2

  tlsPlot.dh.3d(las, cbind(meds-a*height, meds+a*height), cbind(meds, meds+n*r), r, clear, wired, col)

  # cols = if(hasField(las, 'gpstime')) lidR:::set.colors(las$gpstime, las$gpstime %>% unique %>% length %>% height.colors) else 'darkgrey'
  # if(clear) clear3d() # else rgl.open()
  # bg3d('black') ; axes3d(col='white')
  # lines3d(data.frame(meds+q-a*height, meds+q+a*height) %>% t, color='darkred', lwd=3)
  # lines3d(data.frame(meds, meds+q) %>% t, color='blue', lwd=3)
  # rgl.points(pt3d, color=cols)
  #
  # cyl = cylinder3d(rbind(meds+q-a*height, meds+q+a*height), radius=r, sides=36)
  # wire3d(cyl, color='white', lwd=3)
  # pan3d(2)

  # ptsx = cos(seq(0,2*pi,length.out = 36)) * r
  # ptsy = sin(seq(0,2*pi,length.out = 36)) * r
  # ptsz = height
  # ptring = lapply(ptsz, function(x) cbind(ptsx,ptsy,x)) %>% do.call(what=rbind)

  ptring = cbind(
    cos(seq(0,2*pi,length.out = 36)) * r,
    sin(seq(0,2*pi,length.out = 36)) * r,
    0
  )

  v = xprod(a,c(0,0,1))
  vx = matrix(c(0,-v[3],v[2],v[3],0,-v[1],-v[2],v[1],0),ncol=3,byrow = T)
  rmat = diag(3) + vx + (vx %*% vx)/(1+ as.double(a %*% c(0,0,1)))

  ptring = ptring %*% rmat
  ptring[,1] = ptring[,1] + meds[1]
  ptring[,2] = ptring[,2] + meds[2]
  ptring[,3] = ptring[,3] + meds[3]

  if(!is.null(las)) tlsPlot.dh.2d(las, ptring, r)

  # tid = if(hasField(las, 'TreeID')) paste('tree', las$TreeID[1], '-') else ''
  # pt3d[,1:2] %>% plot(pch=20, cex=1, asp=1, main=paste(tid ,'d =',round(r*200,2),'cm'), col=cols)
  # abline(v=seq(min(las$X),max(las$X),.02),col='grey',lty=2)
  # abline(h=seq(min(las$Y),max(las$Y),.02),col='grey',lty=2)
  # ptring[,1:2] %>% lines(col='darkred', lwd=2)
  # points((meds+q)[1], (meds+q)[2], col='darkred', cex=2, pch=3)
}

tlsPlot.dh.circle = function(las, x, y, r, clear=F, wired=T, col='white'){
  pt3d = las@data[,.(X,Y,Z)]
  cbase = c(x,y,min(las$Z))
  ctop  = c(x,y,max(las$Z))

  ccen  = c(x,y,range(las$Z) %>% mean)
  crad  = c(x+r,y,range(las$Z) %>% mean)

  tlsPlot.dh.3d(las, cbind(cbase, ctop), cbind(ccen, crad), r, clear, wired, col)

  # cols = if(hasField(las, 'gpstime')) lidR:::set.colors(las$gpstime, las$gpstime %>% unique %>% length %>% height.colors) else 'darkgrey'
  # if(clear) clear3d() # else rgl.open()
  # bg3d('black') ; axes3d(col='white')
  # lines3d(data.frame(cbase, ctop) %>% t, color='darkred', lwd=3)
  # lines3d(data.frame(ccen, crad) %>% t, color='blue', lwd=3)
  # rgl.points(pt3d, color=cols)
  #
  # cyl = cylinder3d(rbind(cbase, ctop), radius=r, sides=36)
  # wire3d(cyl, color='white', lwd=3)
  # pan3d(2)

  ptring = data.frame(x + cos(seq(0,2*pi,length.out = 36)) * r,
                      y + sin(seq(0,2*pi,length.out = 36)) * r)

  if(!is.null(las)) tlsPlot.dh.2d(las, ptring, r)

  # tid = if(hasField(las, 'TreeID')) paste('tree', las$TreeID[1], '-') else ''
  # pt3d[,1:2] %>% plot(pch=20, cex=1, asp=1, main=paste(tid, 'd =',round(r*200,2),'cm'), col=cols)
  # abline(v=seq(min(las$X),max(las$X),.02),col='grey',lty=2)
  # abline(h=seq(min(las$Y),max(las$Y),.02),col='grey',lty=2)
  # ptring %>% lines(col='darkred', lwd=2)
  # points(x,y, col='darkred', cex=2, pch=3)
}

tlsPlot.dh.3d = function(las, rings, rVec, r, clear=T, wired=T, col='white'){
  cols = if(hasField(las, 'gpstime')) lidR:::set.colors(las$gpstime, las$gpstime %>% unique %>% length %>% height.colors) else 'darkgrey'

  if(clear) clear3d() # else rgl.open()
  bg3d('black') ; axes3d(col='white')
  lines3d(t(rings), color='darkred', lwd=3)
  lines3d(t(rVec), color='blue', lwd=3)
  # rgl.points(las %>% las2xyz, color=cols)

  cyl = cylinder3d(t(rings), radius=r, sides=36)
  if(wired) wire3d(cyl, color=col, lwd=3) else shade3d(cyl, color=col)
  # pan3d(2)
}

tlsPlot.dh.2d = function(las, ring, r, col='darkred'){
  cols = if(hasField(las, 'gpstime')) lidR:::set.colors(las$gpstime, las$gpstime %>% unique %>% length %>% height.colors) else 'darkgrey'

  tid = if(hasField(las, 'TreeID')) paste('tree', las$TreeID[1], '-') else ''
  las@data[,.(X,Y)] %>% plot(pch=20, cex=1, asp=1, main=paste(tid, 'd =',round(r*200,2),'cm'), col=cols)
  abline(v=seq(min(las$X),max(las$X),.02),col='grey',lty=2)
  abline(h=seq(min(las$Y),max(las$Y),.02),col='grey',lty=2)
  ring[,1:2] %>% lines(col=col, lwd=2)
  points(mean(ring[,1]), mean(ring[,2]), col=col, cex=2, pch=3)
}

bringToOrigin = function(las, x){
  if(!(length(x) == 2 && is.numeric(x)[1])) return(las)

  if(class(las)[1] == 'LAS'){
    las@data$X = las@data$X - x[1]
    las@data$Y = las@data$Y - x[2]
  }else{
    las$X = las$X - x[1]
    las$Y = las$Y - x[2]
  }

  return(las)
}

add_treeIDs = function(x, las, ...){
  las = bringToOrigin(las, x)
  minz = min(las$Z) - 0.5
  if(class(las)[1]){
    las = las@data[.(TreeID,X,Y)]
  }
  las = las[TreeID > 0,.(X=mean(X), Y=mean(Y)), by='TreeID']
  las %$% text3d(X,Y,minz,TreeID, ...)
}

add_treeMap = function(){}

add_treePoints = function(x, las, color_func=pastel.colors, ...){
  isLAS(las)
  if(!hasField(las, 'TreeID')){
    stop('TreeID field not found.')
  }
  las = filter_poi(las, TreeID > 0)
  las = bringToOrigin(las, x)
  colors = las$TreeID %>% unique %>% length %>% color_func
  colors = lidR:::set.colors(las$TreeID, colors)
  las@data %$% rgl.points(X,Y,Z,color=colors,...)
}

add_stemPoints = function(x, las, ...){
  isLAS(las)
  if(!hasField(las, 'Stem')){
    stop('Stem field not found.')
  }
  las = filter_poi(las, Stem)
  las = bringToOrigin(las, x)
  las@data %$% rgl.points(X,Y,Z,...)
}

add_stemSegments = function(x, las, stems_data_table, cylinders=FALSE, color='white'){

  isLAS(las)

  if(!hasField(las, 'Segment')){
    stop('las must be the output from stemPoints')
  }

  if(!(hasAttribute(stems_data_table, 'single_stem_dt') || hasAttribute(stems_data_table, 'multiple_stems_dt'))){
    stop('stems_data_table must be the output from stemSegmentation')
  }

  stems_data_table = bringToOrigin(stems_data_table, x)

  nms = colnames(stems_data_table)
  if('DX' %in% nms){
    vals = stems_data_table[,c('X','Y','Radius', 'DX', 'DY')]
    colnames(vals) = c('x', 'y', 'radius', 'ax', 'ay')
  }else if('rho' %in% nms){
    vals = stems_data_table[,c('rho', 'theta', 'phi', 'alpha', 'Radius')]
    colnames(vals) %<>% tolower
  }else{
    vals = stems_data_table[,c('X', 'Y', 'Radius')]
    colnames(vals)[3] %<>% tolower
  }

  las = filter_poi(las, Stem)
  las@data

  if('TreeID' %in% colnames(stems_data_table)){

  }else{
    for(i in 1:nrow(stems_data_table)){
      seg = stems_data_table$Segment[i]
      temp = filter_poi(las, Segment == seg)
      tlsPlot.dh(temp, vals[i,], F, T, 'white')
    }
  }
}

add_tlsInventory = function(x){}
