#import foram data
#import chem data
#import ntubes data
#harmonise data - may be easier to do heavy lifting in SQL

#regression ntubes vs o2
#nforams vs ntubes
#nforams vs O2 - forams at different depths
#Any species associated with tubes/deep ? use PRC
#ordination of 

#size fraction by depth/species




library(vegan)
library(mgcv)
library(RODBC)
library(rioja)

con<-odbcConnectAccess("PES_DB_3Feb2011.mdb")
foramF<-sqlQuery(con,"select * from PES_DB_foraminifera_species_data where DAT>20050000")
levels(foramF$Size)[levels(foramF$Size)== "063-125"]<-"63-125"
names(foramF)<-make.names(names(foramF))
dim(foramF)

close(con)

foram<-xtabs(Number.split~paste(Station_code,Replicate,Slice_numeric, sep=":")+CodeForam, data=foramF, subset=!Station_code%in%c("RC5","RC8","RC9"),drop.unused.levels = TRUE)
dim(foram)
foram<-as.data.frame(unclass(foram))

foram.meta<-as.data.frame(matrix(unlist(strsplit(rownames(foram),":")), ncol=3, byrow=T))
names(foram.meta)<-c("station", "replicate", "depth")
table(foram.meta[,c(1,3)])
foram.meta$mO2<-chem$mO2[sapply(foram.meta$station,function(g){n=which(rownames(chem)==g);ifelse(length(n)==0,NA,n)})]

#################
plot.depth.grad<-function(x,y,depth, n, lines=T, legend=F){
  xnew<-data.frame(x1=seq(-1,5,.1))
#browser()
  x11(4.75,4.75)
  ly<-layout(cbind(1:2),heights=c(.25,.75))
  layout.show(ly)
  par(mar=c(0,3,2,1),mgp=c(1.5,.5,0))
  if(any(y[depth==-1]>0)){
    plot(x[depth==-1],y[depth==-1], col=1, log="y", main=n, xlim=range(x, na.rm=T), ylab="Forams per polytube", xlab="", xaxt="n")
    if(lines){
      x1<-x[depth==-1]                      
      mod<-glm(y[depth==-1]~x1+I(x1^2)+I(x1^3), family=quasipoisson)
      lines(xnew[,1],predict(mod,newdata=xnew,type="response"), col=1)
    }
  }else{
     plot(x[depth==-1],y[depth==-1], col=1, log="y", main=n, ylim=c(1,10), ylab="Forams per polytube", xlab="", xaxt="n")
  }

  
  par(mar=c(3,3,0,1),mgp=c(1.5,.5,0))
  plot(x[depth!=-1],y[depth!=-1], col=depth[depth!=-1], log="y", xlim=range(x, na.rm=T), xlab=expression(minO[2]), ylab="Forams")
  
  if(lines){
    x1<-x[depth==0.5]
    mod<-glm(y[depth==0.5]~x1+I(x1^2)+I(x1^3), family=quasipoisson)
    lines(xnew[,1],predict(mod,newdata=xnew,type="response"), col=2)
  
    x1<-x[depth==1.5]
    mod<-glm(y[depth==1.5]~x1+I(x1^2)+I(x1^3), family=quasipoisson)
    lines(xnew[,1],predict(mod,newdata=xnew,type="response"), col=3)
  
    x1<-x[depth==3]
    mod<-glm(y[depth==3]~x1+I(x1^2)+I(x1^3), family=quasipoisson)
    lines(xnew[,1],predict(mod,newdata=xnew,type="response"), col=4)
  
    x1<-x[depth==5]
    mod<-glm(y[depth==5]~x1+I(x1^2)+I(x1^3), family=quasipoisson)
    lines(xnew[,1],predict(mod,newdata=xnew,type="response"), col=5)
  }
  if(legend)legend("bottomright", legend=levels(depth),pch=1, col=1:length(levels(depth)), cex=.7)
}


plot.depth.grad(foram.meta$mO2,rowSums(foram),foram.meta$depth, "No Forams", legend=T)

plot.depth.grad(foram.meta$mO2,foram$STAIFUS,foram.meta$depth, "STAIFUS")
plot.depth.grad(foram.meta$mO2,foram$NONILAB,foram.meta$depth, "NONILAB")
plot.depth.grad(foram.meta$mO2,foram$CASSLAE,foram.meta$depth, "CASSLAE")
plot.depth.grad(foram.meta$mO2,foram$TEXTEAR,foram.meta$depth, "TEXTEAR")
plot.depth.grad(foram.meta$mO2,foram$BULIMAR,foram.meta$depth, "BULIMAR")
plot.depth.grad(foram.meta$mO2,foram$LEPTSCO,foram.meta$depth, "LEPTSCO")
plot.depth.grad(foram.meta$mO2,foram$MICRHYA,foram.meta$depth, "MICRHYA")
plot.depth.grad(foram.meta$mO2,foram$GLOBAUR,foram.meta$depth, "GLOBAUR")
plot.depth.grad(foram.meta$mO2,foram$NONITUR,foram.meta$depth, "NONITUR")

    
#####
xforam63<-xtabs(Number.split~paste(Station_code,Replicate,Slice_numeric, sep=":")+CodeForam, data=foramF[foramF$Size=="63-125",], subset=!Station_code%in%c("RC5","RC8","RC9"),drop.unused.levels = F)
xforam125<-xtabs(Number.split~paste(Station_code,Replicate,Slice_numeric, sep=":")+CodeForam, data=foramF[foramF$Size=="125-500",], subset=!Station_code%in%c("RC5","RC8","RC9"),drop.unused.levels = F)


foram125<-matrix(0,ncol=ncol(xforam125),nrow=nrow(foram))
rownames(foram125)<-rownames(foram)
colnames(foram125)<-colnames(xforam125)
foram125<-as.data.frame(foram125)
foram63<-foram125

foram63[rownames(foram63)%in%rownames(xforam63),]<-xforam63
foram125[rownames(foram125)%in%rownames(xforam125),]<-xforam125
foram63$Sum<-rowSums(foram63)
foram125$Sum<-rowSums(foram125)


plot.depth.grad.pie<-function(x,sp,depth, n, lines=T, legend=F){
  y63<-foram63[[sp]]
  y125<-foram125[[sp]]
  y<-foram[[sp]]  
  depth<-jitter(as.numeric(as.character(depth)))
  plot(x,depth, ylim=c(5.5,0), type="n", main=n)
  sapply(1:length(y63),function(i){
    if(sum(c(y63[i],y125[i]))>0){
      if(y63[i]==0)y63[i]=0.000000001
      if(y125[i]==0)y125[i]=0.000000001
      floating.pie(xpos=x[i],ypos=depth[i], x=c(y63[i],y125[i]), col=c(0,1), radius=.1)
    }
  })
  
}



#plot.depth.grad.pie(foram.meta$mO2,"STAIFUS",foram.meta$depth, "STAIFUS")

lapply(c("STAIFUS","NONILAB","CASSLAE","TEXTEAR","BULIMAR","LEPTSCO","MICRHYA","GLOBAUR","NONITUR"),function(n){
  x11()
  plot.depth.grad.pie(foram.meta$mO2,n,foram.meta$depth, n)

})




library(plotrix)
surf<-function (x, depth, Y , col = "red", thinplate = TRUE, main, nlevels = 10, labcex = 0.6, cex = 1, ...) 
{
    miss<-is.na(x)|depth==-1
    x<-x[!miss]
    depth<-depth[!miss]
    Y<-Y[!miss,]
    GRID = 31
    require(mgcv) || stop("Requires package 'mgcv'")
    x1 <- x
    x2 <- depth
    y=Y[,1]/rowSums(Y)
#    mod <- gam(y ~ s(x1, x2, k = 10), family = quasibinomial, weights = rowSums(Y))
     mod <- gam(y ~ s(x1, k=4) + s(x2, k=4), family = quasibinomial, weights = rowSums(Y))
    xn1 <- seq(min(x1), max(x1), len = GRID)
    xn2 <- seq(min(x2), max(x2), len = GRID)
    newd <- expand.grid(x1 = xn1, x2 = xn2)
    fit <- predict(mod, type = "response", newdata = as.data.frame(newd))
    plot(x,depth, ylim=c(5.5,0), type="n")
    sapply(1:nrow(Y),function(i){
      if(sum(Y[i,])>0){
        if(Y[i,1]==0)Y[i,1]=0.000000001
        if(Y[i,2]==0)Y[i,2]=0.000000001
        floating.pie(xpos=x[i],ypos=jitter(depth)[i], x=Y[i,], col=c(0,1), radius=sqrt(sum(Y[i,]))/sqrt(max(rowSums(Y)))*.1)
      }
    })

    title(main = paste(main,paste(paste(c("O2", "Depth"),round(anova(mod)$s.pv,3)), collapse=" ")))
    levels <- pretty(range(fit, finite = TRUE), nlevels)
    contour(xn1, xn2, matrix(fit, nrow = GRID), col = col, add = TRUE, levels = levels, labcex = labcex, drawlabels = !is.null(labcex) && labcex > 0)
    anova(mod)
 }

#surf(foram.meta$mO2,as.numeric(as.character(foram.meta$depth)),Y=cbind(foram63$"NONITUR",foram125$"NONITUR"), main="NONITUR")


lapply(c("STAIFUS","NONILAB","CASSLAE","TEXTEAR","BULIMAR","LEPTSCO","MICRHYA","GLOBAUR","NONITUR", "Sum"),function(n){
  x11()
  surf(foram.meta$mO2,as.numeric(as.character(foram.meta$depth)),Y=cbind(foram63[[n]],foram125[[n]]), main=n)
  dev.print(png,file=paste(n,"size.png"), height=500,width=500)
})
#GlobAur has very few obs with both size fractions - numerical instability

keep<-foram.meta$depth!=-1&!is.na(foram.meta$O2)
prc.response<-decostand(foram[keep,], "total")
prc.response<-foram[keep,]
prc.response<-log(prc.response[,colSums(prc.response>0)>1]+1)
prc.treat<-foram.meta$depth[keep]
prc.var<-cut(foram.meta$O2[keep], breaks=seq(0,5,.5)-0.000001) 
#prc.var<-rep(cut(chem$O2, breaks=seq(0,5,.5))[keep],2)
levels(prc.var)<-seq(.25,4.75,.5)                                        
prc.mod<-prc(prc.response, prc.treat,prc.var)
prc.mod               
anova(prc.mod, strata = prc.var, first=TRUE)
x11();par(mar=c(3,3,1,3), mgp=c(1.5,.5,0))
plot(prc.mod, cex=.6, air=1.01, select=colSums(prc.response>0)>5, xlab="O2")


prc.mod<-prc(prc.response, prc.var,prc.treat)
prc.mod               
anova(prc.mod, strata = prc.var, first=TRUE)
x11();par(mar=c(3,3,1,3), mgp=c(1.5,.5,0))
plot(prc.mod, cex=.6, air=1.01, select=colSums(prc.response>0)>5, xlab="Depth")

#depth optima 
###for low intermediate and high O2 fit WA model and extract optima

wa<-WA(sqrt(foram/rowSums(foram))[,colSums(foram>0)>4], as.numeric(foram.meta$depth), tolDW = T)
opt<-wa$coef
opt<-opt[order(opt[,1]),]
x11();par(mar=c(3,8.5,1,.5), mgp=c(1.5,.5,0), tcl=-.4)
plot(opt[,1], 1:nrow(opt), yaxt="n", xlim=range(c(opt[,1]-opt[,2],opt[,1]+opt[,2])), pch=16, ylab="", ylim=c(0,nrow(opt)+1), yaxs="i", xlab="Depth")
arrows(opt[,1]-opt[,2], 1:nrow(opt), x1 = opt[,1]+opt[,2],  code = 0)
axis(2,at=1:nrow(opt),rownames(opt), las=1, cex.axis=.65)


#Thin to 1 replicate so unbiased. Use same species list for all Ox concentrations
low<-foram.meta$mO2<1&!is.na(foram.meta$mO2)
med<-foram.meta$mO2>=1&foram.meta$mO2<3&!is.na(foram.meta$mO2)
high<-foram.meta$mO2>=3&!is.na(foram.meta$mO2)

lowspp<-colSums(foram[low,]>0)>0&colSums(foram>0)>6
medspp<-colSums(foram[med,]>0)>0&colSums(foram>0)>6
highspp<-colSums(foram[high,]>0)>0&colSums(foram>0)>6

low.opt<-med.opt<-high.opt<-opt
low.opt[low.opt>0]<-NA
med.opt[med.opt>0]<-NA
high.opt[high.opt>0]<-NA

waL<-WA(sqrt(foram)[low,lowspp], as.numeric(foram.meta$depth)[low], tolDW = T)
low.opt<-waL$coef
waM<-WA(sqrt(foram)[med,medspp], as.numeric(foram.meta$depth)[med], tolDW = T)
med.opt<-waM$coef
waH<-WA(sqrt(foram)[high,highspp], as.numeric(foram.meta$depth)[high], tolDW = T)
high.opt<-waH$coef

allspp<-sort(unique(c(rownames(low.opt),rownames(med.opt),rownames(high.opt))))
opt<-matrix(NA,nrow=length(allspp), ncol=6)
rownames(opt)<-allspp
opt[allspp%in%rownames(low.opt),1:2]<-low.opt
opt[allspp%in%rownames(med.opt),3:4]<-med.opt
opt[allspp%in%rownames(high.opt),5:6]<-high.opt

opt<-opt[order(apply(opt,1,max, na.rm=T)),]
x11();par(mar=c(3,8.5,1,.5), mgp=c(1.5,.5,0), tcl=-.4)
plot(opt[,1], 1:nrow(opt)-.1, yaxt="n", xlim=c(0,5), pch=16, ylab="", ylim=c(0,nrow(opt)+1), yaxs="i", xlab="Depth", col=1)
arrows(opt[,1]-opt[,2], 1:nrow(opt)-.1, x1 = opt[,1]+opt[,2],  code = 0, col=1)
points(opt[,3], 1:nrow(opt), yaxt="n", pch=16, col=2)
arrows(opt[,3]-opt[,4], 1:nrow(opt), x1 = opt[,3]+opt[,4],  code = 0, col=2)
points(opt[,5], 1:nrow(opt)+.1, yaxt="n", pch=16, col=3)
arrows(opt[,5]-opt[,6], 1:nrow(opt)+.1, x1 = opt[,5]+opt[,6],  code = 0, col=3)
axis(2,at=1:nrow(opt),rownames(opt), las=1, cex.axis=.65)

  legend("bottomright",legend=c("low O2", "medium O2", "high O2"), pch=16, lty=1, col=1:3) 


       
dat<-read.table("o:/data/pes/polytube species.txt", header=T, sep="\t")
dat$Ntubes
dat<-dat[order(dat$site),]
mtubes<-as.vector(tapply(dat$Ntubes, dat$site, mean))
names(mtubes)<-levels(dat$site)

mtubes<-mtubes[names(mtubes)%in%rownames(chem)]
keep<-rownames(chem)%in%names(mtubes)
chem<-chem[keep,]
fpt<-fpt[keep,]
fp1<-fp1[keep,]
fp2<-fp2[keep,]

#ntubes vs O2
x11(6,6,14)
mod1.to<-glm(mtubes~O2, data=chem, family=quasipoisson)
anova(mod1.to, test="F")
mod2.to<-glm(mtubes~O2+I(O2^2), data=chem, family=quasipoisson)
anova(mod1.to,mod2.to, test="F")
mod3.to<-gam(mtubes~s(O2), data=chem, family=quasipoisson)
anova(mod2.to,mod3.to, test="F")
anova(mod1.to,mod3.to, test="F")
anova(mod3.to)
summary(mod3.to)

plot(mod1.to)
plot(chem$O2,mtubes)
lines(seq(0,5,.1), predict(mod1.to, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=2)
lines(seq(0,5,.1), predict(mod2.to, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=2)
lines(seq(0,5,.1), predict(mod3.to, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=3)



#number of forams vs number tubes
fpert<-rowSums(fpt)/mtubes
quantile(fpert, na.rm=T)
 x11(6,6,14)
plot(mtubes,fpert, xlab="number tubes", ylab="number forams")
mod1<-glm(fpert~mtubes, family=quasipoisson, na.action=na.exclude, subset=is.finite(fpert))
lines(mtubes[is.finite(fpert)],fitted(mod1))
summary(mod1)


#n forams vs O2
plot(fpert~O2, data=chem)
mod3<-glm(fpert ~O2, data=chem, family=quasipoisson, subset=is.finite(fpert))
summary(mod3)
mod4<-glm(fpert~O2+I(O2^2), data=chem, family=quasipoisson,subset=is.finite(fpert))
summary(mod4)
mod5<-gam(fpert~s(O2), data=chem, family=quasipoisson,subset=is.finite(fpert))
summary(mod5)
anova(mod3, test="F")
anova(mod3,mod4, test="F")
anova(mod4,mod5, test="F")

plot(fpert~O2, data=chem, xlab="O2", ylab="no forams per tube")
lines(seq(0,5,.1),predict(mod4, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=2)
lines(seq(0,5,.1),predict(mod5, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=4)

#shallow
plot(rowSums(fp1)~O2, data=chem)
mod3s<-glm(rowSums(fp1)~O2, data=chem, family=quasipoisson)
summary(mod3s)
mod4s<-glm(rowSums(fp1)~O2+I(O2^2), data=chem, family=quasipoisson)
summary(mod4s)
mod5s<-gam(rowSums(fp1)~s(O2), data=chem, family=quasipoisson)
summary(mod5s)
anova(mod3s, test="F")
anova(mod3s,mod4s, test="F")
anova(mod4s,mod5s, test="F")

plot(rowSums(fp1)~O2, data=chem, xlab="O2", ylab="no forams 0-1 cm")
lines(seq(0,5,.1),predict(mod4s, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=2)
lines(seq(0,5,.1),predict(mod5s, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=4)


#deep
plot(rowSums(fp2)~O2, data=chem)
mod3d<-glm(rowSums(fp2)~O2, data=chem, family=quasipoisson)
summary(mod3d)
mod4d<-glm(rowSums(fp2)~O2+I(O2^2), data=chem, family=quasipoisson)
summary(mod4d)
mod5d<-gam(rowSums(fp2)~s(O2), data=chem, family=quasipoisson)
summary(mod5d)
anova(mod3d, test="F")
anova(mod3d,mod4d, test="F")
anova(mod4d,mod5d, test="F")

plot(rowSums(fp2)~O2, data=chem, xlab="O2", ylab="no forams  1-2 cm")
lines(seq(0,5,.1),predict(mod4d, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=2)
lines(seq(0,5,.1),predict(mod5d, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=4)


plot(chem$O2,rowSums(fp1)/(rowSums(fp1)+rowSums(fp2)), ylab="Proportion of forams in upper cm")

x11();par(mfrow=c(3,1))
plot(fpert~O2, data=chem, xlab="O2", ylab="no forams per tube")
lines(seq(0,5,.1),predict(mod5, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=4)

plot(rowSums(fp1)~O2, data=chem, xlab="O2", ylab="no forams 0-1 cm")
lines(seq(0,5,.1),predict(mod5s, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=4)

plot(rowSums(fp2)~O2, data=chem, xlab="O2", ylab="no forams  1-2 cm")
lines(seq(0,5,.1),predict(mod5d, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=4)


#species richness vs n tubes 
quantile(rowSums(fpt))
quantile(rowSums(fp1))
quantile(rowSums(fp2))

min.count=8
rarep<-rarefy(fpt,min.count)
rarep[rowSums(fpt)<min.count]<-NA
rare1<-rarefy(fp1,min.count)
rare2<-rarefy(fp2,min.count)
rare2[rowSums(fp2)<min.count]<-NA



plot(mtubes,rarep, xlab="Number Tubes",ylab="rarefied number species")

modst1<-glm(rarep~mtubes,family=quasipoisson)
anova(modst1, test="F")
modst2<-glm(rarep~mtubes+I(mtubes^2),family=quasipoisson)
anova(modst1,modst2, test="F")
modst3<-gam(rarep~s(mtubes),family=quasipoisson)
anova(modst2,modst3, test="F")
anova(modst1,modst3, test="F")

plot(mtubes,rare, xlab="Number Tubes",ylab="rarefied number species")
lines(0:50,predict(modst1, newdata=data.frame(mtubes=0:50), type="response"), col=2)
lines(0:50,predict(modst2, newdata=data.frame(mtubes=0:50), type="response"), col=2)
lines(0:50,predict(modst3, newdata=data.frame(mtubes=0:50), type="response"), col=4)

#species richness vs O2

x11();par(mfrow=c(3,1), mar=c(0,3,1,1), oma=c(3,0,0,0), mgp=c(1.5,.5,0))
plot(chem$O2,rarep, xlab="",ylab="rarefied number species", xaxt="n")
plot(chem$O2,rare1, xlab="",ylab="rarefied number species", xaxt="n")
plot(chem$O2,rare2, xlab="",ylab="rarefied number species", xaxt="n")
axis(1, outer=T)
title(xlab="O2", outer=T)

modso1<-glm(rarep~O2,family=quasipoisson, data=chem)
anova(modso1, test="F")
modso2<-glm(rarep~O2+I(O2^2),family=quasipoisson, data=chem)
anova(modso1,modso2, test="F")
modso3<-gam(rarep~s(O2),family=quasipoisson, data=chem)
anova(modso2,modso3, test="F")
anova(modso1,modso3, test="F")

plot(chem$O2,rarep, xlab="O2",ylab="rarefied number species")
lines(seq(0,5,.1),predict(modso1, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=1)
lines(seq(0,5,.1),predict(modso2, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=4)
lines(seq(0,5,.1),predict(modso3, newdata=data.frame(O2=seq(0,5,.1)), type="response"), col=4)


mod<-cca(log(fp1[,colSums(fp1>0)>1]+1)~1, data=chem)
screeplot(mod, bstick=T)



spp2<-decostand(spp[rowSums(>3,colSums(spp>0)>1], "hellinger")
(DCA<-decorana(spp2))
plot(DCA)

mod.ca<-cca(spp2)
plot(mod.ca)
(e.fit<-envfit(mod.ca,nenv[nind>3,], permutations=999))                        
plot(e.fit, add=T, p.max=.05)
screeplot(mod.ca, bstick=T)

nmds<-metaMDS(spp[nind>3,colSums(spp>0)>1], zerodist="add")
plot(nmds)
(e.fit<-envfit(nmds,nenv[nind>3,], permutations=999))                        
plot(e.fit, add=T, p.max=.05)

mod1<-cca(spp2~O2, data=nenv[nind>3,])
mod2<-cca(spp2~tubes, data=nenv[nind>3,])
mod3<-cca(spp2~O2+Condition(tubes), data=nenv[nind>3,])
anova(mod1)
anova(mod2)
anova(mod3)
mod4<-cca(spp2~O2+tubes, data=nenv[nind>3,])
plot(mod4)


sapply((1:nrow(fp1))[order(chem$O2)],function(n){
  x11()
  x<-cbind(decostand(fp1,"total")[n,],decostand(fp2,"total")[n,])
  plot(x, main=chem$O2[n])
  text(x[rowSums(x)>0,], labels=colnames(fp1)[rowSums(x)>0])
})

sapply((1:nrow(fp1))[order(chem$O2)],function(n){
  x11()
  x<-cbind(decostand(fp1,"total")[n,],decostand(fpt,"total")[n,])
  plot(x, main=chem$O2[n])
  text(x[rowSums(x)>0,], labels=colnames(fp1)[rowSums(x)>0])
})


keep<-rowSums(fp2)>0
prc.response<-decostand(rbind(fp1,fp2), "total")[keep,]
prc.response<-sqrt(prc.response[,colSums(prc.response>0)>1])
prc.treat<-gl(2,nrow(prc.response)/2,labels=c("surface","deep"))
prc.var<-as.factor(rep(chem$O2,2)[keep]) 
prc.var<-rep(cut(chem$O2, breaks=seq(0,5,.5))[keep],2)
levels(prc.var)<-seq(.25,4.75,.5)                                        
prc.mod<-prc(prc.response, prc.treat,prc.var)
prc.mod
anova(prc.mod, strata = prc.var, first=TRUE)
x11();par(mar=c(3,3,1,3), mgp=c(1.5,.5,0))
plot(prc.mod, cex=.6, air=1.01, select=colSums(prc.response>0)>5, xlab="O2")




keep<-rowSums(fpt)>0
prc.response<-decostand(rbind(fp1,fpt), "total")[keep,]
prc.response<-sqrt(prc.response[,colSums(prc.response>0)>1])
prc.treat<-gl(2,nrow(prc.response)/2,labels=c("surface","tubes"))
prc.var<-as.factor(rep(chem$O2,2)[keep]) 
prc.var<-rep(cut(chem$O2, breaks=seq(0,5,.5))[keep],2)
levels(prc.var)<-seq(.25,4.75,.5)                                        

cca.mod<-cca(prc.response)
plot(cca.mod, type="n")
points(cca.mod, disp="sp", pch="+")
points(cca.mod, disp="sites", col=prc.treat)

cca.mod<-cca(prc.response~(rep(chem$O2,2)[keep])+prc.treat)
plot(cca.mod)
plot(cca.mod, type="n")
points(cca.mod, disp="sp")
points(cca.mod, disp="sites", col=prc.treat, type="l", select=T)
                                                 
prc.mod<-prc(prc.response, prc.treat,prc.var)
prc.mod
anova(prc.mod, strata = prc.var, first=TRUE)
x11();par(mar=c(3,3,1,3), mgp=c(1.5,.5,0))
plot(prc.mod, cex=.5, air=1.01, select=colSums(prc.response>0)>5, xlab="O2")



#taxonomic distance between paired samples
  keep<-rowSums(fpt)>20
t1.d<-sapply((1:nrow(fpt))[keep], function(n){
  vegdist(decostand(rbind(fpt[n,],fp1[n,]), "total"))
})
plot(chem$O2[keep],t1.d)

keep<-rowSums(fp2)>3
t2.d<-sapply((1:nrow(fp2))[keep], function(n){
  vegdist(decostand(rbind(fp2[n,],fp1[n,]), "total"))
})
plot(chem$O2[keep],t2.d)

