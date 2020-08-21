xprod = function(a, b){
  x = c(
    a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1]
  )
  return(x);
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

tfBruteForceCoordinates = function(dt, ax, ay){
  nm = names(dt)
  for(i in 1:nrow(dt)){
    rmat = rotationMatrix(-ax[i]*pi/180, 0, 0) %*% rotationMatrix(0, -ay[i]*pi/180, 0)
    vals = as.matrix(dt[i,]) %*% rmat %>% as.double
    for(j in 1:2) dt[i,j] = vals[j]
  }
  names(dt) = nm
  return(dt)
}

.pan3d = function(button=2){
  start <- list()

  begin <- function(x, y) {
    start$userMatrix <<- par3d("userMatrix")
    start$viewport <<- par3d("viewport")
    start$scale <<- par3d("scale")
    start$projection <<- rgl.projection()
    start$pos <<- rgl.window2user( x/start$viewport[3], 1 - y/start$viewport[4], 0.5,
                                   projection=start$projection)
  }

  update <- function(x, y) {
    xlat <- (rgl.window2user( x/start$viewport[3], 1 - y/start$viewport[4], 0.5,
                              projection = start$projection) - start$pos)*start$scale
    mouseMatrix <- translationMatrix(xlat[1], xlat[2], xlat[3])
    par3d(userMatrix = start$userMatrix %*% t(mouseMatrix) )
  }
  rgl.setMouseCallbacks(button, begin, update)
  # cat("Pan set on button", button, "of rgl device",rgl.cur(),"\n")
}

# --- Hidden function copied from lidR's source code (reference: https://github.com/Jean-Romain/lidR):
# From rgl.setMouseCallbacks man page
# nocov start
.pan3d <- function(button=2, dev = rgl::rgl.cur(), subscene = rgl::currentSubscene3d(dev))
{
  start <- list()

  begin <- function(x, y)
  {
    activeSubscene <- rgl::par3d("activeSubscene", dev = dev)
    start$listeners <<- rgl::par3d("listeners", dev = dev, subscene = activeSubscene)

    for (sub in start$listeners)
    {
      init <- rgl::par3d(c("userProjection","viewport"), dev = dev, subscene = sub)
      init$pos <- c(x/init$viewport[3], 1 - y/init$viewport[4], 0.5)
      start[[as.character(sub)]] <<- init
    }
  }

  update <- function(x, y)
  {
    for (sub in start$listeners)
    {
      init <- start[[as.character(sub)]]
      xlat <- 2*(c(x/init$viewport[3], 1 - y/init$viewport[4], 0.5) - init$pos)
      mouseMatrix <- rgl::translationMatrix(xlat[1], xlat[2], xlat[3])
      rgl::par3d(userProjection = mouseMatrix %*% init$userProjection, dev = dev, subscene = sub )
    }
  }
  rgl::rgl.setMouseCallbacks(button, begin, update, dev = dev, subscene = subscene)
}
# nocov end


# --- Hidden function copied from lidR's source code (reference: https://github.com/Jean-Romain/lidR):
set.colors = function (x, palette, trim = Inf, value_index = FALSE){
  if (all(is.na(x)))
    return()
  if (value_index) {
    x[x >= length(palette)] <- length(palette) - 1
    return(palette[x + 1])
  }
  ncolors <- length(palette)
  if (!is.infinite(trim))
    x[x > trim] <- trim
  minx <- min(x, na.rm = T)
  maxx <- max(x, na.rm = T)
  if (maxx - minx == 0) {
    if (!anyNA(x)) {
      colors <- palette[1]
    }
    else {
      colors <- rep(NA_character_, length(x))
      colors[!is.na(x)] <- palette[1]
    }
  }
  else {
    idx <- findInterval(x, seq(minx, maxx, length.out = ncolors))
    colors <- palette[idx]
  }
  return(colors)
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
  is_las = class(las)[1] == 'LAS'
  pt3d = if(is_las) las@data[,.(X,Y,Z)] else las

  rmat = rotationMatrix(ax*pi/180, ay*pi/180, 0)
  temp = as.matrix(pt3d) %*% rmat %>% as.data.table
  names(temp) = names(pt3d)

  rmat = rotationMatrix(-ax*pi/180, 0, 0) %*% rotationMatrix(0, -ay*pi/180, 0)
  cbase = c(x,y,min(temp$Z)) %*% rmat %>% as.double
  ctop  = c(x,y,max(temp$Z)) %*% rmat %>% as.double

  ccen  = c(x,y,range(temp$Z) %>% mean) %*% rmat %>% as.double
  crad  = c(x+r,y,range(temp$Z) %>% mean) %*% rmat %>% as.double

  tlsPlot.dh.3d(las, cbind(cbase, ctop), cbind(ccen, crad), r, clear, wired, col)
  if(!is_las) return(NULL)

  ptring = cbind(
    X = x + cos(seq(0,2*pi,length.out = 36)) * r,
    Y = y + sin(seq(0,2*pi,length.out = 36)) * r,
    Z = range(temp$Z) %>% mean
  ) %*% rmat

  tlsPlot.dh.2d(las, ptring, r)
}

tlsPlot.dh.cylinder = function(las, rho, theta, phi, alpha, r, clear=F, wired=T, col='white'){
  n = c(cos(phi) * sin(theta) , sin(phi) * sin(theta) , cos(theta))
  ntheta = c(cos(phi) * cos(theta) , sin(phi) * cos(theta) , -sin(theta))
  nphi = c(-sin(phi) * sin(theta) , cos(phi) * sin(theta) , 0);
  nphibar = nphi/sin(theta)

  a = ntheta * cos(alpha) + nphibar * sin(alpha)
  q = n*(r+rho)

  is_las = class(las)[1] == 'LAS'
  if(is_las){
    meds = apply(las@data[,.(X,Y,Z)], 2, function(x) sum(range(x))/2) %>% as.double
    pt3d = las@data[,.(X,Y,Z)]
  }else{
    meds = apply(las, 2, mean) %>% as.double
  }

  height = las$Z %>% range %>% diff %>% abs
  height = height/2

  tlsPlot.dh.3d(las, cbind(meds-a*height, meds+a*height), cbind(meds, meds+n*r), r, clear, wired, col)
  if(!is_las) return(NULL)

  # cols = if(hasField(las, 'gpstime')) set.colors(las$gpstime, las$gpstime %>% unique %>% length %>% height.colors) else 'darkgrey'
  # if(clear) clear3d() # else rgl.open()
  # bg3d('black') ; axes3d(col='white')
  # lines3d(data.frame(meds+q-a*height, meds+q+a*height) %>% t, color='darkred', lwd=3)
  # lines3d(data.frame(meds, meds+q) %>% t, color='blue', lwd=3)
  # rgl.points(pt3d, color=cols)
  #
  # cyl = cylinder3d(rbind(meds+q-a*height, meds+q+a*height), radius=r, sides=36)
  # wire3d(cyl, color='white', lwd=3)
  # .pan3d(2)

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

  tlsPlot.dh.2d(las, ptring, r)

  # tid = if(hasField(las, 'TreeID')) paste('tree', las$TreeID[1], '-') else ''
  # pt3d[,1:2] %>% plot(pch=20, cex=1, asp=1, main=paste(tid ,'d =',round(r*200,2),'cm'), col=cols)
  # abline(v=seq(min(las$X),max(las$X),.02),col='grey',lty=2)
  # abline(h=seq(min(las$Y),max(las$Y),.02),col='grey',lty=2)
  # ptring[,1:2] %>% lines(col='darkred', lwd=2)
  # points((meds+q)[1], (meds+q)[2], col='darkred', cex=2, pch=3)
}

tlsPlot.dh.circle = function(las, x, y, r, clear=F, wired=T, col='white'){
  is_las = class(las)[1] == 'LAS'
  if(is_las) pt3d = las@data[,.(X,Y,Z)]
  cbase = c(x,y,min(las$Z))
  ctop  = c(x,y,max(las$Z))

  ccen  = c(x,y,range(las$Z) %>% mean)
  crad  = c(x+r,y,range(las$Z) %>% mean)

  tlsPlot.dh.3d(las, cbind(cbase, ctop), cbind(ccen, crad), r, clear, wired, col)
  if(!is_las) return(NULL)

  # cols = if(hasField(las, 'gpstime')) set.colors(las$gpstime, las$gpstime %>% unique %>% length %>% height.colors) else 'darkgrey'
  # if(clear) clear3d() # else rgl.open()
  # bg3d('black') ; axes3d(col='white')
  # lines3d(data.frame(cbase, ctop) %>% t, color='darkred', lwd=3)
  # lines3d(data.frame(ccen, crad) %>% t, color='blue', lwd=3)
  # rgl.points(pt3d, color=cols)
  #
  # cyl = cylinder3d(rbind(cbase, ctop), radius=r, sides=36)
  # wire3d(cyl, color='white', lwd=3)
  # .pan3d(2)

  ptring = data.frame(x + cos(seq(0,2*pi,length.out = 36)) * r,
                      y + sin(seq(0,2*pi,length.out = 36)) * r)

  tlsPlot.dh.2d(las, ptring, r)

  # tid = if(hasField(las, 'TreeID')) paste('tree', las$TreeID[1], '-') else ''
  # pt3d[,1:2] %>% plot(pch=20, cex=1, asp=1, main=paste(tid, 'd =',round(r*200,2),'cm'), col=cols)
  # abline(v=seq(min(las$X),max(las$X),.02),col='grey',lty=2)
  # abline(h=seq(min(las$Y),max(las$Y),.02),col='grey',lty=2)
  # ptring %>% lines(col='darkred', lwd=2)
  # points(x,y, col='darkred', cex=2, pch=3)
}

tlsPlot.dh.3d = function(las, rings, rVec, r, clear=T, wired=T, col='white'){
  is_las = class(las)[1] == 'LAS'

  if(clear) clear3d() # else rgl.open()
  bg3d('black') ; axes3d(col='white')
  lines3d(t(rings), color='darkred', lwd=3)
  lines3d(t(rVec), color='blue', lwd=3)

  if(is_las){
    cols = if(hasField(las, 'gpstime')) set.colors(las$gpstime, las$gpstime %>% unique %>% length %>% height.colors) else 'darkgrey'
    rgl.points(las %>% las2xyz, color=cols)
  }

  cyl = cylinder3d(t(rings), radius=r, sides=36)
  if(wired) wire3d(cyl, color=col, lwd=3) else shade3d(cyl, color=col)
  # .pan3d(2)
}

#' @importFrom graphics abline lines
tlsPlot.dh.2d = function(las, ring, r, col='darkred'){
  cols = if(hasField(las, 'gpstime')) set.colors(las$gpstime, las$gpstime %>% unique %>% length %>% height.colors) else 'darkgrey'

  tid = if(hasField(las, 'TreeID')) paste('tree', las$TreeID[1], '-') else ''
  las@data[,.(X,Y)] %>% plot(pch=20, cex=1, asp=1, main=paste(tid, 'd =',round(r*200,2),'cm'), col=cols)
  abline(v=seq(min(las$X),max(las$X),.02),col='grey',lty=2)
  abline(h=seq(min(las$Y),max(las$Y),.02),col='grey',lty=2)
  ring[,1:2] %>% lines(col=col, lwd=2)
  points(mean(ring[,1]), mean(ring[,2]), col=col, cex=2, pch=3)
}

#' @rdname tlsPlot
#' @export
add_segmentIDs = function(x, las, ...){
  if(class(las)[1] == 'LAS'){
    if(!hasField(las, 'Segment')) stop('Segment field not found.')
    las = filter_poi(las, Stem == T)
    las = las@data
  }else if(!('Segment' %in% colnames(las))){
    stop('Segment field not found.')
  }

  by = 'Segment'
  if('TreeID' %in% colnames(las)){
    by = c('TreeID', by)
    las = las[TreeID > 0]
  }

  las = las[,.(X=mean(X),Y=mean(Y),Z=mean(Z)),by=by]
  las = bringToOrigin(las, x)
  las %$% text3d(X,Y,Z,Segment, ...)
}

#' @rdname tlsPlot
#' @export
add_treeIDs = function(x, las, ...){
  if(class(las)[1] == 'LAS'){
    if(!hasField(las, 'TreeID')) stop('TreeID field not found.')
    las = las@data
  }else{
    if(!all(c('X','Y','TreeID') %in% colnames(las))){
      stop('X, Y and TreeID fields must be present.')
    }else if(!hasField(las, 'Z')){
      las$Z = 0
    }
  }

  las = bringToOrigin(las, x)
  z = min(las$Z) - .5
  las = las[TreeID > 0,.(X=mean(X), Y=mean(Y)), by='TreeID']
  las %$% text3d(X,Y,z,TreeID, ...)
}

#' @rdname tlsPlot
#' @export
add_treeMap = function(x, las, ...){
  if(!hasAttribute(las, 'tree_map') && !hasAttribute(las, 'tree_map_dt'))
    stop('las is not a tree_map object: check ?treeMap')

  h = 1.3
  if(hasAttribute(las, 'tree_map')){
    if(hasField(las, 'TreePosition')){
      las %>% bringToOrigin(x) %>% las2xyz %>% rgl.points(...)
      return(NULL)
    }
    h = mean(las$Z)
    las %<>% treeMap.positions(F)
  }
  las = bringToOrigin(las,x)
  spheres3d(las$X, las$Y, h, radius=.25, ...)
}


#' @rdname tlsPlot
#' @export
add_treePoints = function(x, las, color_func=pastel.colors, ...){
  isLAS(las)
  if(!hasField(las, 'TreeID')){
    stop('TreeID field not found.')
  }
  las = filter_poi(las, TreeID > 0)
  las = bringToOrigin(las, x)
  colors = las$TreeID %>% unique %>% length %>% color_func
  colors = set.colors(las$TreeID, colors)
  las@data %$% rgl.points(X,Y,Z,color=colors,...)
}

#' @rdname tlsPlot
#' @export
add_stemPoints = function(x, las, ...){
  isLAS(las)
  if(!hasField(las, 'Stem')){
    stop('Stem field not found.')
  }
  las = filter_poi(las, Stem)
  las = bringToOrigin(las, x)
  las@data %$% rgl.points(X,Y,Z,...)
}

#' @rdname tlsPlot
#' @importFrom stats median
#' @export
add_stemSegments = function(x, stems_data_table, color='white', fast=FALSE){

  if(!(hasAttribute(stems_data_table, 'single_stem_dt') || hasAttribute(stems_data_table, 'multiple_stems_dt'))){
    stop('stems_data_table must be the output from stemSegmentation')
  }

  nms = colnames(stems_data_table)
  if('DX' %in% nms){
    stems_data_table = bringToOrigin(stems_data_table, x)
    vals = stems_data_table[,c('X','Y','Radius', 'DX', 'DY')]
    colnames(vals) = c('x', 'y', 'radius', 'ax', 'ay')
    positions = stems_data_table[,.(X,Y,Z=AvgHeight)]
    if(fast) positions = tfBruteForceCoordinates(positions, vals$ax, vals$ay)
  }else if('rho' %in% nms){
    vals = stems_data_table[,c('rho', 'theta', 'phi', 'alpha', 'Radius')]
    colnames(vals) %<>% tolower
    positions = stems_data_table[,.(X=PX,Y=PY,Z=PZ)] %>% bringToOrigin(x)
  }else{
    stems_data_table = bringToOrigin(stems_data_table, x)
    vals = stems_data_table[,c('X', 'Y', 'Radius')]
    colnames(vals)[3] %<>% tolower
    positions = stems_data_table[,.(X,Y,Z=AvgHeight)]
  }

  if(fast){
    spheres3d(positions, radius = stems_data_table$Radius, color=color)
  }else{
    len = ifelse(nrow(stems_data_table) == 1, .15, median(positions$Z[-1] - positions$Z[-nrow(positions)])/2)
    for(i in 1:nrow(stems_data_table)){
      temp = positions[i,]
      temp = rbind(temp,temp)
      temp$Z[1] = temp$Z[1] - len
      temp$Z[2] = temp$Z[2] + len
      tlsPlot.dh(temp, vals[i,], F, T, color)
    }
  }
}

#' @rdname tlsPlot
#' @export
add_tlsInventory = function(x, inventory_data_table, color='white', fast=FALSE){

  if(!hasAttribute(inventory_data_table, 'tls_inventory_dt')){
    stop('inventory_data_table must be the output from tlsInventory')
  }

  nms = colnames(inventory_data_table)
  if('DX' %in% nms){
    inventory_data_table = bringToOrigin(inventory_data_table, x)
    vals = inventory_data_table[,c('X','Y','Radius', 'DX', 'DY')]
    colnames(vals) = c('x', 'y', 'radius', 'ax', 'ay')
    positions = inventory_data_table[,.(X,Y,Z=h_radius)]
    if(fast) positions = tfBruteForceCoordinates(positions, vals$ax, vals$ay)
  }else if('rho' %in% nms){
    vals = inventory_data_table[,c('rho', 'theta', 'phi', 'alpha', 'Radius')]
    colnames(vals) %<>% tolower
    positions = inventory_data_table[,.(X=PX,Y=PY,Z=PZ)] %>% bringToOrigin(x)
  }else{
    inventory_data_table = bringToOrigin(inventory_data_table, x)
    vals = inventory_data_table[,c('X', 'Y', 'Radius')]
    colnames(vals)[3] %<>% tolower
    positions = inventory_data_table[,.(X,Y,Z=h_radius)]
  }

  if(fast){
    spheres3d(positions, radius = inventory_data_table$Radius, color=color)
  }else{
    len = .15
    for(i in 1:nrow(inventory_data_table)){
      temp = positions[i,]
      temp = rbind(temp,temp)
      temp$Z[1] = temp$Z[1] - len
      temp$Z[2] = temp$Z[2] + len
      tlsPlot.dh(temp, vals[i,], F, T, color)
    }
  }

  for(i in 1:nrow(positions)){
    positions[i,] %$% arrow3d(c(X,Y,0), c(X,Y,inventory_data_table$H[i]), color=color, width = .1, thickness = .1, s=.05, type='rotation', n=18)
  }

}
