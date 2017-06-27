## Sample analysis
## Adult income data from UCI repository
rm(list=ls())
# setwd('C:/Study/My projects/Depth-model selection/Codes')
source('misc_functions.R')

library(R.utils)
library(fda.usc)
library(parallel)
library(leaps)

## overhead function
get.err = function(n.sig, X, sig, mseq, nboot, n.iter, ncores, corr=F){
  
  n = nrow(X)
  p = ncol(X)
  
  ## Set up coef vector
  beta = c(rep(1,n.sig),rep(0,p - n.sig))
  X.names = paste0(rep("x",p), 1:p)
  beta.names = paste0(X.names[1:n.sig], collapse="")
  form = as.formula(paste0("y~", paste(X.names, collapse="+")))
  
  ## calculate score matrices
  P = solve(t(X) %*% X)
  H = P %*% t(X)
  
  ## Calculate true value of Cn
  Cn.true = rep(0, p+1)
  beta.sim = my.mvrnorm(1e4, beta, P)
  Cn.true[p+1] = mean(mdepth.TD(beta.sim, beta.sim)$dep)
  
  for(j in 1:p){
    jbeta.sim = beta.sim
    jbeta.sim[,j] = 0
    Cn.true[j] = mean(mdepth.TD(jbeta.sim, beta.sim)$dep)
  }
  
  ## AIC and BIC: backward deletion and all subset
  n.AIC.bk = rep(0,2)
  n.AIC.all = rep(0,2)
  n.BIC.bk = rep(0,2)
  n.BIC.all = rep(0,2)
  
  set.seed(10052015)
  for(iter in 1:n.iter){
    y = X %*% beta + sig*rnorm(n)
    
    lmod = lm(form, data=data.frame(X))
    selmod = step(lmod, direction="backward", k=log(n), trace=0)
    selnames = names(selmod$coef)[-1]
    best.ind = sort(which(selnames %in% X.names))
    name.paste = paste0(X.names[best.ind], collapse="")
    
    if(name.paste==beta.names){
      n.BIC.bk[1] = n.BIC.bk[1]+ 1
    }
    if(grepl(beta.names, name.paste)){
      n.BIC.bk[2] = n.BIC.bk[2]+ 1
    }
    
    selmod = step(lmod, direction="backward", trace=0)
    selnames = names(selmod$coef)[-1]
    best.ind = sort(which(selnames %in% X.names))
    name.paste = paste0(X.names[best.ind], collapse="")
    
    if(name.paste==beta.names){
      n.AIC.bk[1] = n.AIC.bk[1]+ 1
    }
    if(grepl(beta.names, name.paste)){
      n.AIC.bk[2] = n.AIC.bk[2]+ 1
    }
    
    # all subsets regression BIC
    subsetObj = summary(regsubsets(X,y,nvmax=p))
    bicvals = subsetObj$bic
    best.ind = which(subsetObj$which[which.min(bicvals),-1])
    name.paste = paste0(X.names[sort(best.ind)], collapse="")
    
    if(name.paste==beta.names){
      n.BIC.all[1] = n.BIC.all[1]+ 1
    }
    if(grepl(beta.names, name.paste)){
      n.BIC.all[2] = n.BIC.all[2]+ 1
    }
    
    # all subsets regression AIC
    aicvals = subsetObj$cp
    best.ind = which(subsetObj$which[which.min(aicvals),-1])
    name.paste = paste0(X.names[sort(best.ind)], collapse="")
    
    if(name.paste==beta.names){
      n.AIC.all[1] = n.AIC.all[1]+ 1
    }
    if(grepl(beta.names, name.paste)){
      n.AIC.all[2] = n.AIC.all[2]+ 1
    }
  }
  
  ## depth model selection function
  loopfun = function(m){
    require(fda.usc)
    
    n.Cn = rep(0,2)
    sd = 1
    Cn.mat = matrix(0, ncol=p+1, nrow=n.iter)
    
    pb = txtProgressBar(0,n.iter)
    set.seed(10052015)
    for(iter in 1:n.iter){
      y = X %*% beta + sig*rnorm(n)
      form = as.formula(paste0("y~", paste(X.names, collapse="+")))
      
      ## full model constants
      beta.hat = H %*% y
      r = y - X %*% beta.hat

      # Cn for full model estimates: Fn is approx by Fn^b1
      SSPmat.d = matrix(0, nrow=nboot, ncol=p+1)
      beta.mat.m = matrix(0 ,nrow=nboot, ncol=p)
      for(i in 1:nboot){
        samp = 1
        while(length(unique(samp)) < p){
          samp = sample(1:n, m, replace=T)
        }
        H.samp = solve(t(X[samp,]) %*% X[samp,]) %*% t(X[samp,])
        beta.mat.m[i,] = H.samp %*% y[samp]
      }
      SSPmat.d[,p+1] = mdepth.TD(beta.mat.m, beta.mat.m)$dep
      
      # marginal models
      for(j in 1:p){
        jbeta.mat.m = beta.mat.m
        jbeta.mat.m[,j] = 0
        SSPmat.d[,j] = mdepth.TD(jbeta.mat.m, beta.mat.m)$dep
      }
      
      (Cn.vec = apply(SSPmat.d, 2, mean))
      which.sel = which(Cn.vec[-(p+1)] < Cn.vec[p+1])
      name.Cn = paste(X.names[which.sel], collapse="")
      if(name.Cn==beta.names){
        n.Cn[1] = n.Cn[1] + 1
      }
      
      if(grepl(beta.names, name.Cn)){
        n.Cn[2] = n.Cn[2] + 1
      }
      Cn.mat[iter,] = Cn.vec
      setTxtProgressBar(pb,iter)
    }
    
    close(pb)
    c(n.Cn,
      abs(apply(Cn.mat - matrix(Cn.true, nrow=n.iter, ncol=p+1, byrow=T),2,mean)),
      abs(apply(Cn.mat - matrix(Cn.true, nrow=n.iter, ncol=p+1, byrow=T),2,sd)))
  }
  
  ## Do it multicore!
  system.time(n.Cn.mat <- mclapply(mseq, loopfun, mc.cores=ncores))
  # system.time(n.Cn.mat <- lapply(mseq, loopfun))
  n.Cn.mat = matrix(unlist(n.Cn.mat), ncol=2*p+4, byrow=T)
  
  n.Cn.mat[,1:2] = n.Cn.mat[,1:2]/n.iter
  list(IC.vec = rbind(n.AIC.bk, n.AIC.all,n.BIC.bk, n.BIC.all)/n.iter,
       Cn.mat = cbind(mseq, n.Cn.mat))
}

## Run it!!
set.seed(10052015)
n = 1e2
mseq = seq(11,50,by=1)
p = 10
sig = 1
nboot = 1e2
n.iter = 1e2

d = read.table("../data/charliesim.txt", header=T)
X = as.matrix(d[sample(1:nrow(d),100,replace=F),1:p])

out2 = get.err(n.sig=2, X=X, sig=sig, mseq=mseq, nboot=nboot, n.iter=n.iter, ncores=15)
out4 = get.err(n.sig=4, X=X, sig=sig, mseq=mseq, nboot=nboot, n.iter=n.iter, ncores=15)
out6 = get.err(n.sig=6, X=X, sig=sig, mseq=mseq, nboot=nboot, n.iter=n.iter, ncores=15)
out8 = get.err(n.sig=8, X=X, sig=sig, mseq=mseq, nboot=nboot, n.iter=n.iter, ncores=15)
lmOutMarginalMooN = list(out2,out4,out6,out8)
save(lmOutMarginalMooN, file='lmOutMarginalMooN1.rda')

# out2c = get.err(n.sig=2, X=X, sig=sig, mseq=mseq, nboot=nboot, n.iter=n.iter, ncores=15, corr=T)
# out4c = get.err(n.sig=4, X=X, sig=sig, mseq=mseq, nboot=nboot, n.iter=n.iter, ncores=15, corr=T)
# out6c = get.err(n.sig=6, X=X, sig=sig, mseq=mseq, nboot=nboot, n.iter=n.iter, ncores=15, corr=T)
# out8c = get.err(n.sig=8, X=X, sig=sig, mseq=mseq, nboot=nboot, n.iter=n.iter, ncores=15, corr=T)
# lmOutCorrMarginalMooN = list(out2c,out4c,out6c,out8c)
# save(lmOutCorrMarginalMooN, file='lmOutCorrMarginalMooN.rda')

# sends email when it's done!
# library(mailR)
# sender <- "bstat1@gmail.com"
# recipients <- c("zoom.subha@gmail.com")
# send.mail(from = sender,
#           to = recipients,
#           subject="Subject of the email",
#           body = "Body of the email",
#           smtp = list(host.name = "smtp.gmail.com", port=465,
#                       user.name="bstat1@gmail.com", passwd="jhaatjalasna", ssl=TRUE),
#           authenticate = TRUE,
#           send = TRUE)