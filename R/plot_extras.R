xprod = function(a, b){
  x = c(
    a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1]
  )
  return(x);
}

tlsPlot.dh = function(las, pars, clear=T){
  pars %<>% as.double
  if(length(pars) > 5)
    tlsPlot.dh.cylinder(las, pars[1], pars[2], pars[3], pars[4], pars[5]/200, clear)
  else
    tlsPlot.dh.circle(las, pars[1], pars[2], pars[3]/200, clear)
}

tlsPlot.dh.cylinder = function(las, rho, theta, phi, alpha, r, clear=F){
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

  cols = if(hasField(las, 'gpstime')) lidR:::set.colors(las$gpstime, las$gpstime %>% unique %>% length %>% height.colors) else 'darkgrey'

  if(clear) clear3d() else rgl.open()
  bg3d('black') ; axes3d(col='white')
  lines3d(data.frame(meds+q-a*height, meds+q+a*height) %>% t, color='red', lwd=3)
  lines3d(data.frame(meds, meds+q) %>% t, color='blue', lwd=3)
  rgl.points(pt3d, color=cols)

  cyl = cylinder3d(rbind(meds+q-a*height, meds+q+a*height), radius=r, sides=36)
  wire3d(cyl, color='white', lwd=3)
  pan3d(2)

  ptsx = cos(seq(0,2*pi,length.out = 36)) * r
  ptsy = sin(seq(0,2*pi,length.out = 36)) * r
  ptsz = height

  ptring = lapply(ptsz, function(x) cbind(ptsx,ptsy,x)) %>% do.call(what=rbind)

  v = xprod(a,c(0,0,1))
  vx = matrix(c(0,-v[3],v[2],v[3],0,-v[1],-v[2],v[1],0),ncol=3,byrow = T)

  rmat = diag(c(1,1,1)) + vx + (vx %*% vx)/(1+ as.double(a %*% c(0,0,1)))
  ptring = ptring %*% rmat
  ptring[,1] = ptring[,1] + q[1] + meds[1]
  ptring[,2] = ptring[,2] + q[2] + meds[2]
  ptring[,3] = ptring[,3] + q[3] + meds[3]

  tid = if(hasField(las, 'TreeID')) paste('tree', las$TreeID[1], '-') else ''
  pt3d[,1:2] %>% plot(pch=20, cex=.5, asp=1, main=paste(tid ,'d =',round(r*200,2),'cm'), col=cols)
  abline(v=seq(min(las$X),max(las$X),.02),col='grey',lty=2)
  abline(h=seq(min(las$Y),max(las$Y),.02),col='grey',lty=2)
  ptring[,1:2] %>% lines(col='red', lwd=2)
  points((meds+q)[1], (meds+q)[2], col='red', cex=2, pch=3)
}

tlsPlot.dh.circle = function(las, x, y, r, clear=F){
  pt3d = las@data[,.(X,Y,Z)]
  cbase = c(x,y,min(las$Z))
  ctop  = c(x,y,max(las$Z))

  ccen  = c(x,y,range(las$Z) %>% mean)
  crad  = c(x+r,y,range(las$Z) %>% mean)

  cols = if(hasField(las, 'gpstime')) lidR:::set.colors(las$gpstime, las$gpstime %>% unique %>% length %>% height.colors) else 'darkgrey'

  if(clear) clear3d() else rgl.open()
  bg3d('black') ; axes3d(col='white')
  lines3d(data.frame(cbase, ctop) %>% t, color='red', lwd=3)
  lines3d(data.frame(ccen, crad) %>% t, color='blue', lwd=3)
  rgl.points(pt3d, color=cols)

  cyl = cylinder3d(rbind(cbase, ctop), radius=r, sides=36)
  wire3d(cyl, color='white', lwd=3)
  pan3d(2)

  ptring = data.frame(x + cos(seq(0,2*pi,length.out = 36)) * r,
                      y + sin(seq(0,2*pi,length.out = 36)) * r)

  tid = if(hasField(las, 'TreeID')) paste('tree', las$TreeID[1], '-') else ''
  pt3d[,1:2] %>% plot(pch=20, cex=.5, asp=1, main=paste(tid, 'd =',round(r*200,2),'cm'), col=cols)
  abline(v=seq(min(las$X),max(las$X),.02),col='grey',lty=2)
  abline(h=seq(min(las$Y),max(las$Y),.02),col='grey',lty=2)
  ptring %>% lines(col='red', lwd=2)
  points(x,y, col='red', cex=2, pch=3)
}
