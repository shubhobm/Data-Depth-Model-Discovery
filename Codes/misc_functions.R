# misc functions
## function to generate from multivariate normal
my.mvrnorm = function(n, mu, Sigma){
  p = length(mu)
  # compute square root of covariance matrix
  eo=eigen(Sigma, symmetric=TRUE)
  sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
  
  # generate random normals from runif by box-muller transform
  rnorm.vec = sqrt(-2*log(runif(n*p)))*cos(2*pi*runif(n*p))
  
  # generate sample matrix
  sample.matrix = matrix(rep(mu, n), nrow=n, byrow=T) +
    matrix(rnorm.vec, nrow=n, ncol=p)%*%sigma.sqrt
  return(sample.matrix)
}

ecmeml1 = function (y, subj, pred, xcol, zcol, vmax, occ, start, maxits = 1000, 
          eps = 1e-04) 
{
  tmp <- table(subj)
  m <- length(tmp)
  nmax <- max(tmp)
  ntot <- length(y)
  pcol <- ncol(pred)
  q <- length(zcol)
  p <- length(xcol)
{
    if (missing(vmax)) {
      vmax <- diag(rep(1, nmax))
      occ <- integer(ntot)
      iflag <- as.integer(1)
    }
    else iflag <- as.integer(0)
  }
storage.mode(vmax) <- "double"
{
  if (!missing(start)) {
    beta <- start$beta
    sigma2 <- start$sigma2
    xi <- start$psi/start$sigma2
    storage.mode(beta) <- "double"
    storage.mode(xi) <- "double"
    storage.mode(sigma2) <- "double"
    sflag <- as.integer(1)
  }
  else {
    beta <- numeric(p)
    xi <- matrix(0, q, q)
    sigma2 <- 0
    sflag <- as.integer(0)
  }
}
# if(trace==F){
#   cat("Performing ECME...")
# }
now <- proc.time()
err <- 0
tmp <- .Fortran("ecmeml", ntot, as.integer(subj), m, ist = integer(m), 
                ifin = integer(m), as.integer(occ), nmax, vmax, w = array(0, 
                                                                          c(nmax, nmax, m)), vinv = array(0, c(nmax, nmax, 
                                                                                                               m)), pcol, as.double(pred), q, as.integer(zcol), 
                ztvinv = array(0, c(q, nmax, m)), ztvinvz = array(0, 
                                                                  c(q, q, m)), iflag = iflag, err = as.integer(err), 
                msg = integer(1), u = array(0, c(q, q, m)), iter = integer(1), 
                sflag, sigma2 = sigma2, p, as.integer(xcol), beta = beta, 
                as.double(y), delta = rep(0, ntot), xtw = matrix(0, p, 
                                                                 nmax), xtwx = matrix(0, p, p), xtwy = numeric(p), 
                xtwxinv = matrix(0, p, p), wkqq1 = matrix(0, q, q), wkqq2 = matrix(0, 
                                                                                   q, q), xi = xi, wkqnm = array(0, c(q, nmax, m)), 
                b = matrix(0, q, m), cvgd = integer(1), obeta = rep(0, 
                                                                    p), oxi = matrix(0, q, q), maxits = as.integer(maxits), 
                llvec = numeric(as.integer(maxits)), eps = as.double(eps), 
                PACKAGE = "lmm")
clock <- proc.time() - now
# if(trace==F){
#   cat("\n")
# }
iter <- tmp$iter
msg <- tmp$msg
{
  if (msg == 1) 
    warning("Supplied V <- i matrix is not positive definite")
  else if (msg == 2) 
    warning("GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank")
  else if (msg == 3) 
    warning("Inadequate information to obtain starting value of psi")
  else if (msg == 4) 
    warning("Value of psi became non-pos.def. during iterations")
  else if (msg == 5) 
    warning("t(X)%*%W%*%X became non-pos.def. during iterations")
}
llvec <- tmp$llvec[1:iter]
converged <- tmp$cvgd == as.integer(1)
cov.beta <- tmp$xtwxinv * tmp$sigma2
b.hat <- tmp$b
cov.b <- tmp$u * tmp$sigma2
psi <- tmp$xi * tmp$sigma2
beta <- tmp$beta
if (!is.null(dimnames(pred)[[2]])) {
  colnames <- dimnames(pred)[[2]]
  names(beta) <- colnames[xcol]
  dimnames(psi) <- list(colnames[zcol], colnames[zcol])
}
list(beta = beta, sigma2 = tmp$sigma2, psi = psi, converged = converged, 
     iter = iter, loglik = llvec, cov.beta = cov.beta, b.hat = b.hat, 
     cov.b = cov.b)
}
