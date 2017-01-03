#Diversity indices - performance with marine biodiversity
#bias in shannon
library(entropy)
library(vegan)

make.div<-function(x){
  d<-list()
  N<-sapply(x,function(n)sum(n))
  d$shannon<-sapply(x,function(n)diversity(n))
  d$shannon.bc<-sapply(x,function(n)div(n))
  d$expshannon<-exp(d$shannon)
  d$expshannon.bc<-exp(d$shannon.bc)
  d$simpson<-sapply(x,function(n)diversity(n, "simpson"))
  d$invsimpson<-sapply(x,function(n)diversity(n, "invsimpson"))
  
  d$fishers.alpha<-sapply(x,function(n)fisher.alpha(n))
  d$ES100<-sapply(x,function(n)rarefy(n, 100))
  
  d$ES100[N<100]<-NA
  d$ES50<-sapply(x,function(n)rarefy(n, 50))
  d$ES50[N<50]<-NA

  d$margalef<-sapply(x,function(n)log(length(n))/log(sum(n)))
  d$SN<-sapply(x,function(n)log(length(n))/log(log(sum(n))))
#  d$SNc<-sapply(x,function(n)log(length(n))/log(log(sum(n)))*sum(n)/(sum(n)+5))
  d$S<-sapply(x,function(n)length(n))
  d$logS<-sapply(x,function(n)log(length(n)))
  
 # d$shannonE<-d$shannon/log(d$S)
 # d$shannon.bcE<-d$shannon.bc/log(d$S)
 # d$simpsonE<-d$simpson/d$S

  as.data.frame(d)
}
get.div<-function(x, xx){
  d<-make.div(x)
  summ<-lapply(d,function(e){
    x<-tapply(e,xx,function(f)c(mean(f,na.rm=T),sd(f,na.rm=T),quantile(f,probs=c(.025,.975),na.rm=T)))
    x<-as.data.frame(matrix(unlist(x), ncol=4, byrow=T))
    names(x)<-c("mean","sd","x.025","x.975")
    x
  })
  list(summ=summ,divesities=d)
}  


plot.div<-function(mul,org){
x11();par(mfrow=c(4,4),mar=c(3,3,1,1), mgp=c(1.5,.5,0))
mapply(function(mul,n,org){
  matplot(unique(xs),mul[,-2],main=n, type="l", col=c(2,4,4), xlim=c(0,500), xlab="N", ylab="diversity")
  abline(h=org[[1]]*c(1,.95,.9), col=c(1,2,4), lty=2)
  }
  ,mul=mul, n=names(mul), org=org)
}

plot.coefv<-function(mul){
  matplot(unique(xs),sapply(mul,function(x)x[,2]/x[,1]*100), type="l", ylab="Coefficient of error", xlab="N", xlim=c(0,500))
  legend("topright", legend=names(mul),col=1:6, lty=1:5)
}


runstuff<-function(sam,xs){
  resampled<-sapply(xs,function(n)t(table(sample(rep(1:length(sam),sam), size=n, replace=T))))
  
  resampled1000<-sapply(rep(1000,100),function(n)t(table(sample(rep(1:length(sam),sam), size=n, replace=T))))
  
  re1000.div<-get.div(resampled1000, rep(1000,100))
  re.div<-get.div(resampled, xs)
  list(resampled=resampled,resampled1000=resampled1000,re1000.div=re1000.div,re.div=re.div)
}


xs<-rep(seq(20,1000,10),100)

  tmp<-sapply(xs,function(n)t(table(sample(rep(1:ncol(foram8),foram8[1,,drop=F]), size=n, replace=T))))
  
  tmp1000<-sapply(rep(1000,100),function(n)t(table(sample(rep(1:ncol(foram8),foram8[1,,drop=F]), size=n, replace=T))))
  
  tmp1000.div<-get.div(tmp1000, rep(1000,100))
  tmp1.div<-get.div(tmp, xs)

f1<-runstuff(foram8[1,], xs)
plot.div(mul=f1$re.div[[1]], org=f1$re1000.div[[1]])
x11();plot.coefv(mul=f1$re.div[[1]][1:13])

f6<-runstuff(foram8[6,], xs)
plot.div(mul=f6$re.div[[1]], org=f6$re1000.div[[1]])
x11();plot.coefv(mul=f6$re.div[[1]][1:13])


f16<-runstuff(foram8[16,], xs)
plot.div(mul=f16$re.div[[1]], org=f16$re1000.div[[1]])
x11();plot.coefv(mul=f16$re.div[[1]][1:13])

x11();plot.coefv(mul=f16$re.div[[1]][1:13])

all.div<-lapply(1:nrow(fall),function(n)runstuff(fall[n,], xs)$re.div$summ)

x11(4,4);par(mar=c(3,3,1,1), mgp=c(1.5,.5,0))
matplot(unique(xs),sapply(all.div,function(a){x=a$expshannon.bc;x[,2]/x[,1]*100}), type="l", ylab="Coefficient of error", xlab="Number individuals", xlim=c(0,500), col=1, lty=1)



########plot bias, variance for different samples - check constant pattern
##### pca of diversity indices based on constant count


#for all obs with counts >1000 calculate diverity @1000 & @50 calculate bias50-1000 and % of @1000
#boxplots of bias50-1000 and plots of bias vs @1000
x50<-rep(c(50,1000),each=100)
f8<-as.data.frame(unclass(t(foram8)))
bias<-lapply(f8[,colSums(f8)>1000],function(sam){
  resampled<-sapply(x50,function(n)t(table(sample(rep(1:length(sam),sam), size=n, replace=T))))
  re.div<-get.div(resampled, x50)$summ
  sapply(re.div,function(d)c(bias=(d[1,1]-d[2,1])/d[2,1],div=d[,1],cv=d[,2]/d[,1]*100))
})



x11();par(mar=c(7,3,1,1),mgp=c(2,.5,0), tcl=-.4)
boxplot(lapply(colnames(bias[[1]])[-8],function(n)sapply(bias,function(S)S[1,n])),names=colnames(bias[[1]])[-8], las=2, ylab=expression(bias[50-1000]))
abline(h=0, col="grey80", lty=2)

x11();par(mfrow=c(4,4),mar=c(3,3,1,1),mgp=c(1.5,.5,0))
lapply(colnames(bias[[1]])[-8],function(n){
  x<-t(sapply(bias,function(S)(S[2:1,n])))
  plot(x, main=n, xlab=expression(Diversity[1000]), ylab=expression(bias[50-1000]), ylim=range(sapply(bias,function(S)S[1,]), na.rm=T))
  abline(h=0, col="grey80", lty=2)
})


x11();par(mfrow=c(4,4),mar=c(3,3,1,1),mgp=c(1.5,.5,0))
lapply(colnames(bias[[1]])[-8],function(n){
  x<-t(sapply(bias,function(S)(S[c(2,4),n])))
  plot(x, main=n, xlab=expression(Diversity[1000]), ylab="cv50", ylim=range(sapply(bias,function(S)S[4,]), na.rm=T))
})

x11();par(mfrow=c(4,4),mar=c(3,3,1,1),mgp=c(1.5,.5,0))
lapply(colnames(bias[[1]])[-8],function(n){
  x<-t(sapply(bias,function(S)(S[c(2,5),n])))
  plot(x, main=n, xlab=expression(Diversity[1000]), ylab="cv1000", ylim=range(sapply(bias,function(S)S[5,]), na.rm=T))
})

div50<-sapply(colnames(bias[[1]])[-8],function(n){
  sapply(bias,function(S)(S[2,n]))
  })

round(cor(div50),2)
pca50<-prcomp(div50,scale =T)
x11()
biplot(pca50)
x11()
screeplot(pca50)



div1000<-sapply(colnames(bias[[1]])[-8],function(n){
  sapply(bias,function(S)(S[3,n]))
  })
round(cor(div1000),2)
pca1000<-prcomp(div1000,scale =T)
x11()
biplot(pca1000)
x11()
screeplot(pca1000)

apply(div50,2,cor,scores(rda(div1000,scale =T),disp="sites", choice=1))

mapply(cor,as.data.frame(div50),as.data.frame(div1000))

##################
#macro

m8<-as.data.frame(unclass(t(macro8g)))
mbias<-lapply(m8[,colSums(m8)>500],function(sam){
  resampled<-lapply(x50,function(n)t(table(sample(rep(1:length(sam),sam), size=n, replace=T))))
  re.div<-get.div(resampled, x50)$summ
  sapply(re.div,function(d)c(bias=(d[1,1]-d[2,1])/d[2,1],div=d[,1],cv=d[,2]/d[,1]*100))
})

x11();par(mfrow=c(4,4),mar=c(3,3,1,1),mgp=c(1.5,.5,0))
lapply(colnames(mbias[[1]])[-8],function(n){
  x<-t(sapply(mbias,function(S)(S[2:1,n])))
  plot(x, main=n, xlab=expression(Diversity[1000]), ylab=expression(bias[50-1000]), ylim=range(sapply(bias,function(S)S[1,]), na.rm=T))
  abline(h=0, col="grey80", lty=2)
})


x11();par(mfrow=c(4,4),mar=c(3,3,1,1),mgp=c(1.5,.5,0))
lapply(colnames(mbias[[1]])[-8],function(n){
  x<-t(sapply(mbias,function(S)(S[c(2,4),n])))
  plot(x, main=n, xlab=expression(Diversity[1000]), ylab="cv50", ylim=range(sapply(bias,function(S)S[4,]), na.rm=T))
})

x11();par(mfrow=c(4,4),mar=c(3,3,1,1),mgp=c(1.5,.5,0))
lapply(colnames(mbias[[1]])[-8],function(n){
  x<-t(sapply(mbias,function(S)(S[c(2,5),n])))
  plot(x, main=n, xlab=expression(Diversity[1000]), ylab="cv1000", ylim=range(sapply(bias,function(S)S[5,]), na.rm=T))
})

