ZIMN.gibbs=function(dat,ngibbs,covs,burnin,prior.var){
  nobs=nrow(dat)
  
  tmp=sort(unique(dat$y))
  tmp1=tmp[tmp!=0]
  nclass=length(tmp1)
  class1=1:nclass
  if (sum(tmp==0)>0)  uni=data.frame(y=c(0,tmp1),y1=c(0,class1)) #if zero exists
  if (sum(tmp==0)==0) uni=data.frame(y=tmp1,y1=class1)           #if zero does not exist
  
  dat$ind=1:nrow(dat)
  dat1=merge(dat,uni,all=T); dim(dat); dim(dat1)
  dat2=dat1[order(dat1$ind),c('y1',covs)]
  
  #get initial values
  z=(dat2$y1-mean(dat2$y1))/sd(dat2$y1)
  w=rep(NA,nobs); cond=dat2$y1==0
  w[cond]=-1; w[!cond]=1
  
  b=class1[1:(nclass-1)]+0.1
  b=(b-mean(dat2$y1,na.rm=T))/sd(dat2$y1,na.rm=T)

  #get useful stuff to sample betas
  cond=dat2$y1!=0
  xmat.tmp=data.matrix(dat2[cond,covs])
  npar=ncol(xmat.tmp)
  betas=rep(0,npar)
  xtx=t(xmat.tmp)%*%xmat.tmp
  prec=xtx+diag(1/prior.var,npar)
  var1.betas=solve(prec)
  xmat.betas=data.matrix(dat2[,covs])

  #get useful stuff to sample alpha
  xmat.alpha=data.matrix(cbind(1,dat2[,covs]))
  alpha=rep(0,npar+1)
  xtx=t(xmat.alpha)%*%xmat.alpha
  prec=xtx+diag(c(1/10,rep(1/prior.var,npar)),npar+1)
  var1.alpha=solve(prec)
  
  #gibbs sampler
  store.alpha=matrix(NA,ngibbs,npar+1) #includes intercept
  store.b=matrix(NA,ngibbs,nclass-1)
  store.betas=matrix(NA,ngibbs,npar) #does not include intercept
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)
    z=sample.z(y=dat2$y1,xmat.betas=xmat.betas,betas=betas,b=b,class1=class1,nobs=nobs,nclass=nclass)
    betas=sample.betas(z=z,xmat.betas=xmat.betas,var1.betas=var1.betas,y=dat2$y1)
    b=sample.b(z=z,y=dat2$y1,nclass=nclass,class1=class1)
    w=sample.w(y=dat2$y1,xmat.alpha=xmat.alpha,alpha=alpha)
    alpha=sample.alpha(w=w,xmat.alpha=xmat.alpha,var1.alpha=var1.alpha)
    
    #store results
    store.b[i,]=b
    store.betas[i,]=betas
    store.alpha[i,]=alpha
  }
  
  after.burn=burnin:ngibbs
  list(b=store.b[after.burn,],
       betas=store.betas[after.burn,],
       alpha=store.alpha[after.burn,])
}
