# rm(list=ls(all=TRUE))
library('mvtnorm')
set.seed(1)

setwd('U:\\modeling abundance\\git_ZIMN')
source('multinom functions.R')
source('multinom main func.R')
dat=read.csv('fake data.csv',as.is=T)
dat$x1=dat$x; dat$x2=dat$x^2;

ngibbs=50000
covs=c('x1','x2')
prior.var=1

res=multinom.gibbs(dat=dat,ngibbs=ngibbs,covs=covs,burnin=5000,prior.var=prior.var)

par(mfrow=c(1,1))
ebreaks=apply(res$b,2,mean)
rango=range(c(true.breaks,ebreaks))
plot(true.breaks,ebreaks,xlim=rango,ylim=rango)
lines(rango,rango)

par(mfrow=c(2,2))
plot(res$betas[,1],type='l')
abline(h=true.b1,col='red')

plot(res$betas[,2],type='l')
abline(h=true.b2,col='red')

plot(res$alpha[,1],type='l')
abline(h=true.a0,col='red')

plot(res$alpha[,2],type='l')
abline(h=true.a1,col='red')
