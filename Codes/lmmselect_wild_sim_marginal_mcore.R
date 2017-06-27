## Sample analysis
## Adult income data from UCI repository
rm(list=ls())
# setwd('d:/Study/My projects/Depth-model-selection/Codes')
source('misc_functions.R')

library(R.utils)
library(fda.usc)
library(parallel)
library(lmm)

## overhead function
get.err = function(beta, sig, D, ni, m, sdm, nboot, n.iter, ncores){
  
  library(R.utils)
  library(fda.usc)
  library(lmm)
  
  ## set up data
  n = ni*m
  # sdm = seq(1,10, by=.1)
  sdm = 1:24
  p = length(beta)
  pZ = 4
  X = cbind(1,matrix(runif(n*(p-1), -2, 2), ncol=p-1))
  Z = X[,1:pZ]
  subj = rep(1:m, rep(ni,m))
  
  ## calculate score matrices and true value of Cn
  Cn.true = rep(0, p+1)
  
  W.true = matrix(0, ncol=n, nrow=n)
  for(im in 1:m){
    iind = ((im-1)*ni+1):(im*ni)
    Zi = Z[iind,]
    W.true[iind,iind] = sig * diag(ni) + Zi %*% D %*% t(Zi)
  }
  W.true = solve(W.true)
  beta.sim = my.mvrnorm(1e4, beta, solve(t(X) %*% W.true %*% X))
  Cn.true[p+1] = mean(mdepth.TD(beta.sim, beta.sim)$dep)
  
  for(j in 1:p){
    Xj = X[,-j]
    jbeta.sim = matrix(0, ncol=p, nrow=1e4)
    jbeta.sim[,-j] = my.mvrnorm(1e4, beta[-j], solve(t(Xj) %*% W.true %*% Xj))
    Cn.true[j] = mean(mdepth.TD(jbeta.sim, beta.sim)$dep)
  }
    
  ## Set up coef vector
  beta.names = "x2x3"
  X.names = c('x1','x2','x3','x4','x5','x6','x7','x8','x9','x10')
      
  set.seed(11092015)
  ## depth model selection function
  loopfun = function(sd){
    require(fda.usc)
    n.Cn = rep(0,5)
    Cn.mat = matrix(0, ncol=p+1, nrow=n.iter)
    
    set.seed(11092015)
    pb = txtProgressBar(0,n.iter)
    
    for(iter in 1:n.iter){
      
      ## Generate samples
      y = rep(0,n)
      for(im in 1:m){
        iind = ((im-1)*ni+1):(im*ni)
        y[iind] = my.mvrnorm(1, X[iind,] %*% beta, W.true[iind,iind])
      }
      
      ## full model constants
      lmmfull = ecmeml1(y=y, subj=subj, pred=X, xcol=1:p, zcol=1:pZ)
      
      W = matrix(0, nrow=n, ncol=n)
      for(im in 1:m){
        iind = ((im-1)*ni+1):(im*ni)
        Zi = Z[iind,]
        W[iind,iind] = with(lmmfull, solve(sigma2 * diag(ni) + Zi %*% psi %*% t(Zi)))
      }
      beta.hat = lmmfull$beta
      H = lmmfull$cov.beta %*% t(X) %*% W
      r = y - X %*% beta.hat; r = r - mean(r)
      
      # matrix of full model estimates: Fn is approx by Fn^b1
      SSPmat.d = matrix(0, nrow=nboot, ncol=p+1)
      
      # depth model selection
      beta.mat = matrix(0 ,nrow=nboot, ncol=p)
      for(i in 1:nboot){
        iresid = as.matrix(sd * rep(rnorm(m), rep(ni,m)) * r, ncol=1)
        beta.mat[i,] = as.numeric(beta.hat) + as.numeric(H %*% iresid)
      }
      
      beta.mat1 = matrix(0 ,nrow=nboot, ncol=p)
      for(i in 1:nboot){
        iresid = as.matrix(sd * rep(rnorm(m), rep(ni,m)) * r, ncol=1)
        beta.mat1[i,] = as.numeric(beta.hat) + as.numeric(H %*% iresid)
      }
      SSPmat.d[,p+1] = mdepth.TD(beta.mat1, beta.mat)$dep
      
      ## now marginal models
      for(j in 1:p){
        
        jbeta.mat = matrix(0 , nrow=nboot, ncol=p)
        for(i in 1:nboot){
          iresid = sd * rep(rnorm(m), rep(ni,m))* (y - X[,-j] %*% beta.hat[-j])
          iresid = iresid - mean(iresid)
          jbeta.mat[i,-j] = beta.hat[-j] + H[-j,] %*% iresid
        }

        SSPmat.d[,j] = mdepth.TD(jbeta.mat, beta.mat)$dep
      }
      
      Cn.vec = apply(SSPmat.d, 2, mean)
      which.sel = which(Cn.vec[-(p+1)] < Cn.vec[p+1])
      name.Cn = paste(X.names[which.sel], collapse="")

      if(name.Cn==beta.names){
        n.Cn[1] = n.Cn[1] + 1
      }
      
      if(grepl(beta.names, name.Cn)){
        n.Cn[2] = n.Cn[2] + 1
      }
      
      n.Cn[3] = n.Cn[3] + length(which.sel) # model size
      n.Cn[4] = n.Cn[4] + sum(c(1,4:10) %in% which.sel)/length(which.sel) # false positive rate
      n.Cn[5] = n.Cn[5] + sum((2:3) %in% (1:10)[-which.sel])/(10-length(which.sel)) # false negative
      
      Cn.mat[iter,] = Cn.vec
      setTxtProgressBar(pb,iter)    
    }
    close(pb)
    c(n.Cn, abs(apply(Cn.mat - Cn.true,2,mean)), apply(Cn.mat - Cn.true,2,sd))
  }
  
  ## Do it multicore!
  system.time(n.Cn.mat <- mclapply(sdm, loopfun, mc.cores=ncores))
  n.Cn.mat = matrix(unlist(n.Cn.mat), ncol=2*p+7, byrow=T) 
  
  n.Cn.mat[,1:5] = n.Cn.mat[,1:5]/n.iter
  cbind(sdm, n.Cn.mat)
}

## Run it!!
set.seed(11092015)
p = 10
sig = 1
nboot = 1e2
n.iter = 1e2
beta = c(0,rep(1,2),rep(0,p-3)) # change this so that length 9
beta.names = "x2x3"
X.names = c('x1','x2','x3','x4','x5','x6','x7','x8','x9','x10')

D = matrix(c(9,   4.8, 0.6, 0,
             4.8, 4,   1,   0,
             0.6, 1,   1,   0,
             0,   0,   0,   0),
           ncol=4, byrow=T)

out150c = get.err(beta=beta, sig=sig, D=D, ni=5, m=30, nboot=nboot, n.iter=n.iter, ncores=8)
save(out150c, file='lmmout150marginal.rda')

out600c = get.err(beta=beta, sig=sig, D=D, ni=10, m=60, nboot=nboot, n.iter=n.iter, ncores=8)
save(out600c, file='lmmout600marginal.rda')

# library(mailR)
# sender <- "bstat1@gmail.com"
# recipients <- c("zoom.subha@gmail.com")
# send.mail(from = sender,
#           to = recipients,
#           subject="SIMULATION DONE!!!",
#           body = "SIMULATION DONE!!!",
#           smtp = list(host.name = "smtp.gmail.com", port=465,
#                       user.name="bstat1@gmail.com", passwd="jhaatjalasna", ssl=TRUE),
#           authenticate = TRUE,
#           send = TRUE)