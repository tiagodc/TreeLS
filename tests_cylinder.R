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

a = TreeLS:::tlsCylinder(500, rad = .08, h = .5, dev = .015)
a = TreeLS:::las2xyz(a)

# tls = readTLS('inst/extdata/spruce.laz')
# a = lasfilter(tls, Z > 1 & Z < 1.5)

cyl.parameters(a)[[1]]

temp(a)

require(optimx)

init = c(0,pi/2,0,0,0)
optimx(init, cyl.fit, P=a, method = c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin'))

circlefit(a[,1], a[,2])
temp(a)
