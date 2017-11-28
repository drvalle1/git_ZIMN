rm(list=ls(all=TRUE))
set.seed(1)

n=500
x=seq(from=-1,to=1,length.out=n)

true.a0=a0=1
true.a1=a1=-2
true.w=w=rnorm(n,mean=a0+a1*x,sd=1)
plot(x,w)

true.b1=b1=-1
true.b2=b2=0
med=b1*x+b2*(x^2)
true.z=z=rnorm(n,mean=med,sd=1)
plot(x,z)

nclass=20
true.breaks=breaks=sort(rnorm(nclass-1,mean=0,sd=1))

y=rep(NA,n)
cond=w<0
y[cond]=0
cond=z<breaks[1] & w>0; y[cond]=1
for (i in 2:(nclass-1)){
  cond=z<breaks[i] & z>breaks[i-1] & w>0
  y[cond]=i
}
cond=z>breaks[nclass-1] & w>0
y[cond]=nclass

table(y)
hist(y)
plot(x,y)
dat=data.frame(x=x,y=y)

setwd('U:\\modeling abundance\\git_ZIMN')
write.csv(dat,'fake data.csv',row.names=F)