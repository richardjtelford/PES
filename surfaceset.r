library(RODBC)
library(vegan)
library(Hmisc)
library(gdata)
setwd(  "\\\\nturt.uib.no\\gbsrt\\DATA\\pes\\")

con<-odbcConnectAccess("Pes_base_15April2010.mdb")
macroF<-sqlQuery(con, "select * from pes_macrofauna_at_forams_stations")
names(macroF)<-make.names(names(macroF))
stations<-sqlQuery(con,"select * from pes_stations")
names(stations)
foramsF<-sqlQuery(con,"select * from PES_foraminifera_species_data where slice_numeric>0")
names(foramsF)<-make.names(names(foramsF))
dim(foramsF)

chemF<-sqlQuery(con,"SELECT STAS, Chemical_species, Avg(IIf([Value]<-98,0,[value])) AS Val FROM PES_chemistry_data_flat_8nov09 WHERE (((DATE)>20050000) AND ((Slice_numeric)<1 Or (Slice_numeric) Is Null)) GROUP BY STAS, Chemical_species;")



close(con)

###map
library(maps)
library(mapdata)
map("worldHires",xlim=c(7.5,11.5), ylim=c(58,60))
points(stations$"X_COORD_DECIMALDEGREE",stations$"Y_COORD_DECIMALDEGREE", col=2)
box()


#chemistry data
names(chemF)
levels(chemF$STAS)
levels(chemF$STAS)<-trim(levels(chemF$STAS))
levels(chemF$STAS)

table(chemF$Chem)

wantstation<-chemF$STAS%in%levels(macroF$Station)|chemF$STAS%in%levels(foramsF$Station)
sum(wantstation)
chem<-data.frame(unclass(xtabs(Val~STAS+Chemical_species, data=chemF, subset=wantstation&!chemF$Chem%in%c("Chl c3","H2S"), drop=T)))

chem1<-data.frame(unclass(xtabs(rep(1,nrow(chemF))~STAS+Chemical_species, data=chemF, subset=wantstation&!chemF$Chem%in%c("Chl c3","H2S"), drop=T)))
chem[chem1==0]<-NA
dim(chem)
apply(chem,2,function(x)sum(is.na(x)))
apply(chem,1,function(x)sum(is.na(x)))
lapply(chem[,colSums(is.na(chem))>0],function(x)rownames(chem)[is.na(x)])

chem$Chl.a.allomer<-NULL
chem$Pheophorbide<-NULL

chem<-chem[rownames(chem)!="HV16",]
chem<-chem[rownames(chem)!="KV01",]
chem<-chem[rownames(chem)!="KRG",]

sum(is.na(chem))  

#standardise pigments by TOC
names(chem)
pig<- names(chem)%in%c("allo.xanthin","aphanizophyll","beta.carotene","cantha.xanthin","chl.a","chl.a.total..a.allom.","Chl.b","Chl.c2","diadino.xanthin","diato.xanthin","fuco.xanthin","lutein","peridinin","pheo.phytin.a","Pheophytin.b","Pyropheophytin.b", "violaxanthin","zea.xanthin")

chem[,pig]<-chem[,pig]/chem$TOC



x11();par(mfrow=c(5,5), mar=c(3,3,1,1), mgp=c(1.5,.5,0))
mapply(hist,x=chem, main=colnames(chem), xlab="")

names(chem)
chem[,pig]<-sapply(chem[,pig], function(x)log(x+.5*min(x[x>0])))
x11();par(mfrow=c(5,5), mar=c(3,3,1,1), mgp=c(1.5,.5,0))
mapply(hist,x=chem, main=colnames(chem), xlab="")


chem.pca<-prcomp(chem,scale.=T)


x11()
plot(chem.pca)
screeplot(chem.pca, bstick=T)
biplot(chem.pca)

pairs(chem, col=ifelse(rownames(chem)=="GRO50",2,1), gap=0)
plot(chem$TOC, chem$X)
scatter.smooth(chem$TOC, chem$X)


#pigments
RDA<-rda(chem[,pig]~., data=chem[,!pig], scale=T)
RDA0<-rda(chem[,pig]~1, data=chem[,!pig], scale=T)
plot(RDA0)
ordisurf(RDA0,chem$O2, add=T)

plot(chem$O2,scores(RDA0, disp="sites")[,1])
scatter.smooth(chem$O2,scores(RDA0, disp="sites")[,1])
identify(chem$O2,scores(RDA0, disp="sites")[,1], rownames(chem))

RDA1<-step(RDA0, reformulate(names(chem[,!pig])), test="perm")
RDA1
plot(RDA1)
RDAo<-rda(chem[,pig]~O2, data=chem[,!pig], scale=T)
RDAo
plot(RDAo)
RDAoc<-rda(chem[,pig]~O2+TOC, data=chem[,!pig], scale=T)
RDAoc

#macros

macro3r<-xtabs(Number.1~paste(Station_code,Replicate, sep=":")+CodeMacrofauna, data=macroF[macroF$DAT<20050000,])
macro8r<-xtabs(Number.1~paste(Station_code,Replicate, sep=":")+CodeMacrofauna, data=macroF[macroF$DAT>20050000,])

macro3g<-xtabs(Number.1~Station_code+CodeMacrofauna, data=macroF[macroF$DAT<20050000,])
macro8g<-xtabs(Number.1~Station_code+CodeMacrofauna, data=macroF[macroF$DAT>20050000,])



m83<-macro8r[rownames(macro8r)%in%rownames(macro3r),]
m38<-macro3r[rownames(macro3r)%in%rownames(macro8r),]
x11();boxplot(cbind(rowSums(m38),rowSums(m83)))#did the volume sampled remain constant?  --YES, Knockout in 2008
boxplot(cbind(rowSums(m38>0),rowSums(m83>0)))#species
x11();plot(cbind(rowSums(m38),rowSums(m83)), xlab="2003 no ind", ylab="2008 no ind")#species
abline(0,1)
identify(cbind(rowSums(m38),rowSums(m83)), labels=rownames(m83))
x11();plot(cbind(rowSums(m38>0),rowSums(m83>0)), xlab="2003 no spp", ylab="2008 no spp")#species
abline(0,1)
identify(cbind(rowSums(m38>0),rowSums(m83>0)), labels=rownames(m83))


plot(sort(rowSums(macro8), decreasing=T), log="y")
x11();plot(log(colSums(macro8g>0)),log(apply(macro8g,2,function(b)mean(b[b>0]))))#range abundance
abline(lm(log(apply(macro8g,2,function(b)mean(b[b>0])))~log(colSums(macro8g>0))))#range abundance
cor.test(log(colSums(macro8g>0)),log(apply(macro8g,2,function(b)mean(b[b>0]))))
      

decorana(log(macro8[,colSums(macro8>0)>2]+1))
                                           
macro.cca.mod<-(cca(log(macro8[,colSums(macro8>0)>1]+1)))
screeplot(macro.cca.mod, bstick=T)

macro.cca.mod<-(cca(sqrt(decostand(macro8g[rowSums(macro8g)>5,colSums(macro8g>0)>2], "total"))))
screeplot(macro.cca.mod, bstick=T)
plot(macro.cca.mod)
plot(scores(macro.cca.mod, disp="sites")[,1],rowSums(macro8g[rowSums(macro8g)>5,colSums(macro8g>0)>2]))


foram.cca.mod<-(cca(sqrt(decostand(foram, "total"))))
screeplot(foram.cca.mod, bstick=T)

macro8g<-macro8g[rowSums(macro8g)>0,]
macro8gc<-macro8g[rownames(macro8g)%in%rownames(chem),]
chemM<-chem[rownames(chem)%in%rownames(macro8g),]
identical(rownames(macro8gc), rownames(chemM))

mod<-step(cca(log(macro8gc[,colSums(macro8gc>0)>1]+1)~1, data=chemM), reformulate(names(chemM)), test="perm")
mod
x11();plot(mod)

cor(diversity(macro8gc),chemM$O2)


##########
#forams
foram3r<-xtabs(Number.1~paste(Station_code,Replicate, sep=":")+SpeciesForam, data=foramsF[foramsF$DAT<20050000,])
foram8r<-xtabs(Number.1~paste(Station_code,Replicate, sep=":")+SpeciesForam, data=foramsF[foramsF$DAT>20050000,])

dim(foram8r)
dim(macro8r)

quantile(rowSums(foram8r))
quantile(rowSums(macro8r))
div.f<-tapply(diversity(foram8r),unlist(strsplit(rownames(foram8r),":"))[c(T,F)],sd)
div.m<-tapply(diversity(macro8r),unlist(strsplit(rownames(macro8r),":"))[c(T,F)],function(x)sd(x[1:min(3,length(x))]))
div.m<-div.m[names(div.m)%in%names(div.f)]
div.f<-div.f[names(div.f)%in%names(div.m)]
identical(names(div.f),names(div.m))


boxplot(list(foram=div.f,macro=div.m), notch=T)
wilcox.test(div.f,div.m, paired=T)             #try paired analysis?

div.f<-tapply(rarefy(foram8r,10),unlist(strsplit(rownames(foram8r),":"))[c(T,F)],sd)
div.m<-tapply(rarefy(macro8r,10),unlist(strsplit(rownames(macro8r),":"))[c(T,F)],function(x)sd(x[1:min(3,length(x))]))
div.m<-div.m[names(div.m)%in%names(div.f)]
div.f<-div.f[names(div.f)%in%names(div.m)]
identical(names(div.f),names(div.m))

boxplot(list(foram=div.f,macro=div.m))
wilcox.test(div.f,div.m, paired=T)


foram3<-xtabs(Number.1~Station_code+SpeciesForam, data=foramsF[foramsF$DAT<20050000,])
foram8<-xtabs(Number.1~Station_code+SpeciesForam, data=foramsF[foramsF$DAT>20050000,])


foram3r<-xtabs(Number.1~paste(Station_code,Replicate, sep=":")+SpeciesForam, data=foramsF[foramsF$DAT<20050000,])
foram8r<-xtabs(Number.1~paste(Station_code,Replicate, sep=":")+SpeciesForam, data=foramsF[foramsF$DAT>20050000,])

foram3r<-foram3r[rowSums(foram3r>0)>0,]
foram8r<-foram8r[rowSums(foram8r>0)>0,]
                        
foram38<-foram3r[rownames(foram3r)%in%rownames(foram8r),]
foram83<-foram8r[rownames(foram8r)%in%rownames(foram3r),]

plot(rowSums(foram38),rowSums( foram83), xlab="2003 no individuals",ylab="2008 no individuals")
abline(0,1)
text(rowSums(foram38),rowSums( foram83), labels=rownames(foram83))


x11();plot(seq(10,rowSums(foram8)[1]/8,10),sapply(seq(10,rowSums(foram8)[1]/8,10),function(n)rarefy(foram8[1,,drop=F]/8,n)))


mod<-cca(log(foram8[rowSums(foram8)>0,colSums(foram8>0)>1]+1))
screeplot(mod, bstick=T)
plot(mod)
chemf<-chem[rownames(chem)%in%rownames(foram8),]
foram8c<-foram8[rownames(foram8)%in%rownames(chem),]
identical(rownames(foram8c), rownames(chemf))

cor(diversity(foram8c),chemf$O2)

modo2<-cca(

modf<-step(cca(log(foram8c[,colSums(foram8c>0)>1]+1)~1, data=chemf), reformulate(names(chemf)), test="perm")
modf
plot(modf)


macro8gf<-macro8g[rownames(macro8g)%in%rownames(foram8),]
foram8m<-foram8[rownames(foram8)%in%rownames(macro8g),]


plot(vegdist(macro8gf[rowSums(macro8gf)>10,]),vegdist(foram8m[rowSums(macro8gf)>10,]))
mantel(vegdist(macro8gf),vegdist(foram8m))

plot(sort(colSums(foram8), dec=T), log="y")
x11();
plot(sort(colSums(macro8), dec=T), log="y")

 #co=correspondance
 library(cocorresp)

macro8gf<-decostand(macro8g[rownames(macro8g)%in%rownames(foram8),],"total")
foram8m<-decostand(foram8[rownames(foram8)%in%rownames(macro8g),],"total")



identical(rownames(macro8gf),rownames(foram8m))

plot(vegdist(macro8gf[rowSums(macro8gf)>10,]),vegdist(foram8m[rowSums(macro8gf)>10,]))
mantel(vegdist(macro8gf),vegdist(foram8m))

plot(stepacross(vegdist(macro8gf), toolong=.9),stepacross(vegdist(foram8m), toolong=.9))
mantel(stepacross(vegdist(macro8gf), toolong=.9),stepacross(vegdist(foram8m), toolong=.9))

## predictive CoCA using SIMPLS and formula interface
coco.pred <- coca(macro8gf ~ ., data = as.data.frame(foram8m) )
summary(coco.pred)
plot(coco.pred)
## should retain only the useful PLS components for a parsimonious model
## Not run: 
## Leave-one-out crossvalidation - this takes a while
crossval(macro8gf,  as.data.frame(foram8m) )
## so 2 axes are sufficient
## permutation test to assess significant PLS components - takes a while
bp.perm <- permutest.coca(coco.pred, permutations = 99)
bp.perm
summary(bp.perm)

coco.pred <- coca(macro8gf ~ ., data = as.data.frame(foram8m) , n.axes=2)
coco.pred 


plot(diversity(macro8gf),diversity(foram8m))
abline(0,1)

library(rioja)
