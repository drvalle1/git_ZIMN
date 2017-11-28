tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  #qnorm can give some imprecise results
  cond=z<lo;    z[cond] = lo[cond]
  cond=z==-Inf; z[cond] = lo[cond]
  cond=z>hi;    z[cond] = hi[cond]
  cond=z==Inf;  z[cond] = hi[cond]
  z
}
#----------------------------------------------------------------------------------------------
sample.w=function(y,xmat.alpha,alpha){
  media=xmat.alpha%*%alpha
  cond=y==0
  w=rep(NA,length(y))
  w[cond]=tnorm(sum(cond),lo=-Inf,hi=0, mu=media[cond],sig=1)
  w[!cond]=tnorm(sum(!cond),lo=0,hi=Inf,mu=media[!cond],sig=1)
  w
}
#----------------------------------------------------------------------------------------------
sample.alpha=function(w,xmat.alpha,var1.alpha){
  pmedia=var1.alpha%*%t(xmat.alpha)%*%w
  t(rmvnorm(1,pmedia,var1.alpha))
}
#----------------------------------------------------------------------------------------------
sample.z=function(y,xmat.betas,betas,b,class1,nobs,nclass){
  z=rep(NA,nobs)
  media=xmat.betas%*%betas
  
  #------------------------------
  #get z
  cond=y==class1[1] & y!=0
  z[cond]=tnorm(sum(cond),lo=-100,hi=b[1],mu=media[cond],sig=1)
  
  for (i in 2:(nclass-1)){
    cond=y==class1[i] & y!=0
    z[cond]=tnorm(sum(cond),lo=b[i-1],hi=b[i],mu=media[cond],sig=1)
  }

  cond=y==class1[nclass] & y!=0
  z[cond]=tnorm(sum(cond),lo=b[nclass-1],hi=100,mu=media[cond],sig=1)
  z
}
#----------------------------------------------------------------------------------------------
sample.betas=function(z,xmat.betas,var1.betas,y){
  cond=y!=0
  pmedia=var1.betas%*%t(xmat.betas[cond,])%*%z[cond]
  t(rmvnorm(1,pmedia,var1.betas))
}
#----------------------------------------------------------------------------------------------
sample.b=function(z,y,nclass,class1){
  b=rep(NA,nclass-1)
  for (i in 2:nclass){
    cond1=y==class1[i-1] & y!=0
    cond2=y==class1[i]   & y!=0
    lo=max(z[cond1])
    hi=min(z[cond2]);
    # if (lo==hi) hi=hi+0.001
    b[i-1]=runif(1,lo,hi)
  }
  b
}