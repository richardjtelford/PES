source("load.vincent.data.r")

#remove recolinisation sites

forams<-forams[!forams$Station_code%in%c("RC5","RC8","RC9"),]
#correlation of all forams diversity with the environment d1
#correlation of all foram community with the environment  c1

#compare d1 with diversity of different live fractions
#correlation of all foram community with the live fractions
#compare fossilisable diveristy with dead 

#compare c1 with community of different live fractions
#compare fossilisable community with dead


fall<-xtabs(Number.1~Station_code+SpeciesForam, data=forams, drop=T)

#size
unique(forams$Size)
levels(forams$Size)[levels(forams$Size)== "063-125"]<-"63-125"
fbig<-forams$Size=="125-500"
fsmall<-forams$Size=="63-125" 

f63<-xtabs(Number.1~Station_code+SpeciesForam, data=forams[fsmall,], drop=T)
f125<-xtabs(Number.1~Station_code+SpeciesForam, data=forams[fbig,], drop=T)

#depth
shallow<-forams$Slice%in%c( "0-0.5","0-1","0.5-1"  )
deep<-forams$Slice%in%c( "1-1.5" ,    "1-2" ,  "1.5-2" ,"1.5-2.5")
alldep<-forams$Slice%in%c( "0-0.5","0-1","0.5-1","1-1.5","1-2","1.5-2","1.5-2.5")

f01<-   xtabs(Number.1~Station_code+SpeciesForam, data=forams[shallow,], drop=T)
f12<-   xtabs(Number.1~Station_code+SpeciesForam, data=forams[deep,], drop=T)
#f02<-  xtabs(Number.1~Station_code+SpeciesForam, data=forams[alldep,], drop=T)
                                               
#foss
fos<-forams$Species%in%fspp$Species_foram[fspp$Foss=="Fossilizable"]
nofos<-forams$Species%in% fspp$Species_foram[fspp$Foss=="Not fossilizable"|fspp$Foss=="Probably fossilizable"]

ffos<-  xtabs(Number.1~Station_code+SpeciesForam, data=forams[fos,], drop=T)
fnofos<-  xtabs(Number.1~Station_code+SpeciesForam, data=forams[nofos,], drop=T)

flargefos<-xtabs(Number.1~Station_code+SpeciesForam, data=forams[fos&fbig,], drop=T)

#dead
levels(dead$Station_code)[levels(dead$Station_code)%in%c("G40","G50","G60")]<-c("GRO40","GRO50","GRO60")

fdead<- xtabs(Number.1~Station_code+SpeciesForam, data=dead)
fdead<-fdead[order(rownames(fdead)),]
identical(rownames(fall),rownames(fdead))

chem<-chem[rownames(chem)%in%rownames(fall),]

decorana(fall/rowSums(fall))
decorana(sqrt(fall/rowSums(fall)))

#############
fff<-list(  all=fall,  small=f63,  large=f125,  shallow=f01,  deep=f12,  fos=ffos,  nofos=fnofos,  largefos=flargefos,  dead=fdead)
sapply(fff,function(f)sum(rowSums(f)<50))

fdiv<-lapply(fff,function(f)exp(div(f)))
cor.d1<-c(sapply(fdiv,cor,fdiv$all),deadvsfos=cor(fdiv$dead,fdiv$fos))
lapply(fdiv,cor.test,chem$mO2,use="complete")
cor.o2<-c(sapply(fdiv,cor,chem$mO2, use="complete"),  deadvsfos=NA)

ford<-lapply(fff,function(f)rda(sqrt(f/rowSums(f))))
m2<-c(sapply(ford,function(f)sqrt(1-procrustes(ford$all,f, sym=T)$ss)),deadvsfos=sqrt(1-procrustes(ford$fos,ford$dead, sym=T)$ss))

of<-lapply(fff,function(f)rda(sqrt(f/rowSums(f))~mO2,data=chem,  na.action = na.omit))
o2c<-c(sapply(of,function(f)f$CCA$tot.chi/f$tot.chi),deadvsfos=NA)

round(cbind(cor.diversity.O2=cor.o2,cor.diversity.full.diversity=cor.d1, community.var.explained.by.O2=o2c,procrustes.m2=m2),2)

cor.test(div(fall),rarefy(fall,100))
cor(div(fall),rarefy(fall,100))^2

############
############ replicates
fallrep<-   xtabs(Number.1~paste(Station_code, Replicate, sep=":")+SpeciesForam, data=forams[forams$DAT>20050000,], drop=T)
fallrep.meta<-as.data.frame(matrix(unlist(strsplit(rownames(fallrep),":")), ncol=2, byrow=T))
names(fallrep.meta)<-c("station", "replicate")


fallrep.meta$mO2<-chem$mO2[sapply(fallrep.meta$station,function(g)which(rownames(chem)==g))]
fallrep.shan<-exp(div(fallrep))

x11()
plot(fallrep.meta$mO2,fallrep.shan, type="n", xlab=expression(paste("Two Year Minimum ",O[2])), ylab=expression(paste("Foram  ", exp(H*minute[bc]))))
sapply(((1:(length(fallrep.shan)/3))-1)*3+1, function(n)points(fallrep.meta$mO2[n:(n+2)],fallrep.shan[n:(n+2)], type="o", pch=20))
points(chem$mO2,exp(div(fall)), col=2, pch=4)


plot(chem$mO2,tapply((div(fallrep)),fallrep.meta$station,sd))

cbind(chem$mO2,tapply(exp(div(fallrep)),fallrep.meta$station,sd))[chem$mO2>2,-1, drop=F]


#- living in the >63 versus living in the >125 (is there any difference in terms of diversity? so compare results for H' and ES100, can a correlation be a  good analysis?)


x11(height=3,3.3,11);par(mar=c(3,3,.5,.5),mgp=c(1.5,.5,0))
plot(rowSums(f63),rowSums(f125), xlab="No forams 63-125", ylab="No forams 125+")
#plot(colSums(f63)+.1,colSums(f125)+.1, log="xy")
abline(0,1)
min(rowSums(f63))
min(rowSums(f125))
min(rowSums(fall))
plot(rarefy(f63,100),rarefy(f125,100), xlab="ES100 63-125", ylab="ES100 125+")
abline(0,1)
plot(rarefy(fall,100),rarefy(f125,100))
abline(0,1)


plot(div(f63),div(f125), xlab="H' 63-125", ylab="H' 125+")
abline(0,1)
plot(div(fall),div(f125))
abline(0,1)


#- living (>63) in the 0-1cm versus living (>63) in the 0-2cm (is there any difference in terms of diversity? so compare results for H' and ES100, can a correlation be a  good analysis?)

rowSums(f01)
rowSums(f12)
rowSums(f02)

plot(rowSums(f01),rowSums(f12), log="xy", xlab="No forams 0-1 cm", ylab="No forams 1-2 cm")
abline(0,1)

plot(rarefy(f01,100),rarefy(f12,100), xlab="ES100 0-1 cm", ylab="ES100 1-2 cm")
abline(0,1)


plot(div(f01),div(f12), xlab="H' 0-1 cm", ylab="H' 1-2 cm")
abline(0,1)

x11(height=3.5,3.5,11);par(mar=c(3,3,.5,.5),mgp=c(1.5,.5,0))
plot(div(f01),rarefy(f01,100), xlab=expression(paste("Foram  ", H*minute[bc])), ylab="ES100 0-1 cm")
title(main=paste(" r","=",round(cor(rarefy(f01,100),div(f01), use="pair"),2)), line=-1, adj=0)

plot(exp(div(f01)),rarefy(f01,100), xlab=expression(paste("Foram  ", exp(H*minute[bc]))), ylab="ES100 0-1 cm")

#- total living assemblages (>63, 0-2 cm) versus living fossilizable assemblages (>63, 0-2 cm) (is there any difference in terms of diversity? so compare results for H' and ES100, can a correlation be a good analysis?)

plot(rowSums(ffos),rowSums(fnofos), xlab="No Fossilizable", ylab="No non-fossilizable")
abline(0,1)
plot(rarefy(ffos,100),rarefy(fnofos,100), xlab="ES100 Fossilizable", ylab="ES100 non-fossilizable")
abline(0,1)
plot(div(ffos),div(fnofos), xlab="No Fossilizable", ylab="No non-fossilizable")
abline(0,1)

#- dead assemblages (>63, 0-2 cm) versus total living assemblages (>63, 0-2 cm) (is there any difference in terms of diversity? so compare results for H' and ES100, can a correlation be a  good analysis?)

levels(dead$Station_code)[levels(dead$Station_code)%in%c("G40","G50","G60")]<-c("GRO40","GRO50","GRO60")

fdead<- xtabs(Number.1~Station_code+SpeciesForam, data=dead)
fdead<-fdead[order(rownames(fdead)),]
fall<-fall[rownames(fall)%in%rownames(fdead),]
plot(rowSums(fall),rowSums(fdead), xlab="No living", ylab="No dead")
abline(0,1)
plot(rarefy(fall,100),rarefy(fdead,100), xlab="ES100 living", ylab="ES100 dead")
cor(rarefy(fall,100),rarefy(fdead,100))
abline(0,1)
plot(div(fall),div(fdead), xlab="H' living", ylab="H' dead" )
abline(0,1)
cor(div(fall),div(fdead))


#- dead assemblages (>63, 0-2cm) versus living fossilizable assemblages (>63, 0-2 cm) (is there any difference in terms of diversity? so compare results for H' and ES100, can a correlation be a  good analysis?)
ffos<-ffos[rownames(ffos)%in%rownames(fdead),]
plot(rowSums(ffos),rowSums(fdead), xlab="No fossilizable", ylab="No dead")
abline(0,1)
plot(rarefy(ffos,100),rarefy(fdead,100), xlab="ES100 fossilizable", ylab="ES100 dead")
cor(rarefy(ffos,100),rarefy(fdead,100))
abline(0,1)
x11(3.5,3.5);par(mar=c(3,3,1,1),mgp=c(1.5,.5,0))
plot(exp(div(ffos)),exp(div(fdead)), xlab=expression(paste("Fossilisable foram  ", exp(H*minute[bc]))), ylab=expression(paste("Dead foram  ", exp(H*minute[bc]))))
abline(0,1)
title(main=paste(" r","=",round(cor(div(ffos),div(fdead)),2)), line=-1, adj=0)




#In the second part of the paper, I want to make sure foraminifera and macrofauna respond in a similar way to environemntal changes . The thing is that I don't want to be to specific. What I mean is that in an other paper we will deal with specific response, classification in ecological groups, so I don't want to show too much in the actual one.

#Some months ago you sent correlation between H' values for forams and macrofauna and O2 and TOC, do you think it is strong enough to show a similar response to environemtn (actaully could you do it as well with ES100)? or should we present a multivariate analysis?

#Could you also look at possible correlation between betacarotenoid and fauna,  diatoxanthin and fauna, lutein and fauna, zooxanthin and fauna, alloxanthin and fauna, chloa+pheoa and fauna.

#correlation between macro and foram H'/ES100
setdiff(rownames(foram8),rownames(fdead))
setdiff(rownames(foram8),rownames(macro8g))
setdiff(rownames(foram8),rownames(chem))

setdiff(rownames(macro8g),rownames(foram8))
sitesfmc<-intersect(rownames(foram8),intersect(rownames(macro8g),rownames(chem)))
fmc<-foram8[rownames(foram8)%in%sitesfmc,]


mfc<-macro8g[rownames(macro8g)%in%sitesfmc,]
cfm<-chem[rownames(chem)%in%sitesfmc,]
dc<-fdead[rownames(fdead)%in%sitesfmc,]

cfm$xO2<-cfm$O2
cfm$O2<-cfm$MinO2_2_years
cfm$O2[is.na(cfm$O2)]<-cfm$xO2[is.na(cfm$O2)]

x11()
plot(cfm$xO2,cfm$MinO2_2_years)
abline(0,1)
identify(cfm$xO2,cfm$MinO2_2_years, labels=rownames(cfm))

all(rownames(fmc)==rownames(mfc))
all(rownames(fmc)==rownames(cfm))

x11();plot(diversity(fmc),div(fmc))
x11();plot(diversity(mfc),div(mfc))

plot(rowSums(fmc),rowSums(mfc), xlab="No Forams", ylab="No macro")

plot(exp(div(fmc)),exp(div(mfc)), xlab=expression(paste("Foram  ", exp(H*minute[bc]))), ylab=expression(paste("Macrofauna  ", exp(H*minute[bc]))), pch=ifelse(rowSums(mfc)>=50,1,3))
abline(0,1)
title(main=paste(" r","=",round(cor(exp(div(fmc))[rowSums(mfc)>=50],exp(div(mfc))[rowSums(mfc)>=50]),2)), line=-1, adj=1)
cor.test(exp(div(fmc))[rowSums(mfc)>=50],exp(div(mfc))[rowSums(mfc)>=50])





#correlation between forams/macro H'/ES100 and O2/TOC
x11(height=6,6,11);par(mfrow=c(2,2),mar=c(3,3,.5,.5),mgp=c(1.5,.5,0), tcl=-.4)
plot(cfm$O2,div(fmc), xlab=expression(O[2]), ylab=expression(paste("Foram  ", H*minute[bc])))
title(main=paste(" r","=",round(cor(cfm$O2,div(fmc), use="pair"),2)), line=-1, adj=0)
#plot(cfm$O2,exp(diversity(mfc)), xlab="O2", ylab="H' macro")
#title(main=paste(" r","=",round(cor(cfm$O2,diversity(mfc), use="pair"),2)), line=-1, adj=0)
#cor.test(cfm$O2,diversity(mfc), use="pair")
#cor.test(cfm$O2,apply(mfc,1,entropy.ChaoShen), use="pair")
plot(cfm$O2,div(mfc), xlab=expression(O[2]), ylab=expression(paste("Macrofauna  ", exp(H*minute[bc]))), pch=ifelse(rowSums(mfc)>=50,1,3))
title(main=paste(" r","=",round(cor(cfm$O2[rowSums(mfc)>=50],div(mfc)[rowSums(mfc)>=50], use="pair"),2)), line=-1, adj=1)

plot(cfm$O2,exp(div(fmc)), xlab=expression(O[2]), ylab=expression(paste("Foram  ", exp(H*minute[bc]))), cex=cfm$TOC/12)
title(main=paste(" r","=",round(cor(cfm$O2,exp(div(fmc)), use="pair"),2)), line=-1, adj=0)
#plot(cfm$O2,exp(diversity(mfc)), xlab="O2", ylab="H' macro")
#title(main=paste(" r","=",round(cor(cfm$O2,diversity(mfc), use="pair"),2)), line=-1, adj=0)
#cor.test(cfm$O2,diversity(mfc), use="pair")
#cor.test(cfm$O2,apply(mfc,1,entropy.ChaoShen), use="pair")
plot(cfm$O2,exp(div(mfc)), xlab=expression(O[2]), ylab=expression(paste("Macrofauna  ", exp(H*minute[bc]))), pch=ifelse(rowSums(mfc)>=50,1,3))
title(main=paste(" r","=",round(cor(cfm$O2[rowSums(mfc)>=50],exp(div(mfc))[rowSums(mfc)>=50], use="pair"),2)), line=-1, adj=1)

x11()

plot(cfm$O2,exp(div(fmc)), xlab=expression(O[2]), ylab=expression(paste("Foram  ", exp(H*minute[bc]))), cex=cn/5)

plot(cfm$O2,rowSums(mfc), xlab=expression(O[2]), ylab="No. macro")
title(main=paste(" r","=",round(cor(cfm$O2,rowSums(mfc), use="pair"),2)), line=-1, adj=0)

plot(rowSums(mfc), div(mfc)/diversity(mfc), ylab="H'", xlab="No. macro")
                          
plot(cfm$TOC,div(fmc), xlab="TOC", ylab=expression(paste("Foram  ", exp(H*minute[bc]))))
title(main=paste(" r","=",round(cor(cfm$TOC,div(fmc), use="pair"),2)), line=-1, adj=1)
plot(cfm$TOC,div(mfc), xlab="TOC", ylab=expression(paste("Macrofauna  ", exp(H*minute[bc]))))
title(main=paste(" r","=",round(cor(cfm$TOC,div(mfc), use="pair"),2)), line=-1, adj=1)


x11(height=3,6,11);par(mfrow=c(1,2),mar=c(3,3,.5,.5),mgp=c(1.5,.5,0))
plot(cfm$chl.a.total..a.allom.,div(fmc), xlab="chl.a.total..a.allom.", ylab=expression(paste("Foram  ", exp(H*minute[bc]))))
title(main=paste(" r","=",round(cor(cfm$chl.a.total..a.allom.,div(fmc), use="pair"),2)), line=-1, adj=0)
plot(cfm$chl.a.total..a.allom.,exp(div(mfc)), xlab="chl.a.total..a.allom.", ylab=expression(paste("Macrofauna  ", exp(H*minute[bc]))))
title(main=paste(" r","=",round(cor(cfm$chl.a.total..a.allom.,exp(div(mfc)), use="pair"),2)), line=-1, adj=0)

round(sort(sapply(cfm,cor,exp(div(fmc)),use="pair")),2)
round(sort(sapply(cfm,cor,cfm$O2,use="pair")),2)


#RDA of Pigments and O2

plot(cfm, gap=0)
names(cfm)
barplot(cfm$TOC/cfm$TN, names=rownames(cfm), las=2)
plot(cfm$O2,cfm$TOC)

#size vs env
f63mc<-f63[rownames(f63)%in%sitesfmc,]
f125mc<-f125[rownames(f125)%in%sitesfmc,]

x11(height=3,6,11);par(mfrow=c(1,2),mar=c(3,3,.5,.5),mgp=c(1.5,.5,0))
plot(cfm$O2,exp(div(f63mc)), xlab=expression(O[2]), ylab=expression(paste("Foram  ",63-125, mu,"  ", exp(H*minute[bc]))))
title(main=paste(" r","=",round(cor(cfm$O2,exp(div(f63mc)), use="pair"),2)), line=-1, adj=0)
plot(cfm$O2,exp(div(f125mc)), xlab=expression(O[2]), ylab=expression(paste("Foram >125",mu,"  ", exp(H*minute[bc]))))
title(main=paste(" r","=",round(cor(cfm$O2,exp(div(f125mc)), use="pair"),2)), line=-1, adj=0)

x11(height=3,6,11);par(mfrow=c(1,2),mar=c(3,3,.5,.5),mgp=c(1.5,.5,0))
plot(cfm$O2,rowSums(f63mc), xlab=expression(O[2]), ylab="No foram")
title(main=paste(" r","=",round(cor(cfm$O2,rowSums(f63mc), use="pair"),2)), line=-1, adj=0)
plot(cfm$O2,rowSums(f125mc), xlab=expression(O[2]), ylab="No foram")
title(main=paste(" r","=",round(cor(cfm$O2,rowSums(f125mc), use="pair"),2)), line=-1, adj=0)


#depth vs env
f01mc<-f01[rownames(f01)%in%sitesfmc,]
f12mc<-f12[rownames(f12)%in%sitesfmc,]

x11(height=3,6,11);par(mfrow=c(1,2),mar=c(3,3,.5,.5),mgp=c(1.5,.5,0))
plot(cfm$O2,exp(div(f01mc)), xlab=expression(O[2]), ylab=expression(paste("Foram ",0-1,"cm  ", exp(H*minute[bc]))))
title(main=paste(" r","=",round(cor(cfm$O2,exp(div(f01mc)), use="pair"),2)), line=-1, adj=0)
plot(cfm$O2,exp(div(f12mc)), xlab=expression(O[2]), ylab=expression(paste("Foram ",1-2,"cm  ", exp(H*minute[bc]))))
title(main=paste(" r","=",round(cor(cfm$O2,exp(div(f12mc)), use="pair"),2)), line=-1, adj=0)

#fos/nofos vs env
ffosmc<-ffos[rownames(ffos)%in%sitesfmc,]
fnofosmc<-fnofos[rownames(fnofos)%in%sitesfmc,]

x11(height=3,6,11);par(mfrow=c(1,2),mar=c(3,3,.5,.5),mgp=c(1.5,.5,0))
plot(cfm$O2,exp(div(ffosmc)), xlab=expression(O[2]), ylab=expression(paste("Fossilisable foram  ", exp(H*minute[bc]))))
title(main=paste(" r","=",round(cor(cfm$O2,exp(div(ffosmc)), use="pair"),2)), line=-1, adj=0)
plot(cfm$O2,exp(div(fnofosmc)), xlab=expression(O[2]), ylab=expression(paste("None fossilisable foram  ", exp(H*minute[bc]))), pch=ifelse(rowSums(fnofosmc)>=50,1,3))
title(main=paste(" r","=",round(cor(cfm$O2[rowSums(fnofosmc)>=50],exp(div(fnofosmc))[rowSums(fnofosmc)>=50], use="pair"),2)), line=-1, adj=0)

plot(cfm$O2,rowSums(ffosmc)/(rowSums(ffosmc)+rowSums(fnofosmc)), xlab="O2", ylab="prop fos")
cor.test(cfm$O2,rowSums(ffosmc)/(rowSums(ffosmc)+rowSums(fnofosmc)), use="pair")
range(rowSums(ffosmc)/(rowSums(ffosmc)+rowSums(fnofosmc)))

##
#largefossilisable
ff125<-  xtabs(Number.1~Station_code+SpeciesForam, data=forams[forams$DAT>20050000&fbig&fos,], drop=T)
ff125mc<-ffos[rownames(ffos)%in%sitesfmc,]

x11(height=3,3,11);par(mfrow=c(1,1),mar=c(3,3,.5,.5),mgp=c(1.5,.5,0))
plot(cfm$O2,exp(div(ff125mc)), xlab=expression(O[2]), ylab=expression(paste("Fossilisable foram >125",mu,"  ", exp(H*minute[bc]))))
title(main=paste(" r","=",round(cor(cfm$O2,exp(div(ff125mc)), use="pair"),2)), line=-1, adj=0)


#dead
fdeadmc<-fdead[rownames(fdead)%in%sitesfmc,]


x11(height=3,3,11);par(mfrow=c(1,1),mar=c(3,3,.5,.5),mgp=c(1.5,.5,0))
plot(cfm$O2[!rownames(cfm)%in%c("RC5","RC9")],exp(div(fdeadmc)), xlab=expression(O[2]), ylab=expression(paste("Dead foram  ",exp(H[bc]))))
title(main=paste(" r","=",round(cor(cfm$O2[!rownames(cfm)%in%c("RC5","RC9")],exp(div(fdeadmc)), use="pair"),2)), line=-1, adj=0)



#linear model
mod<-glm(exp(div(fmc))~O2+TOC+O2:TOC+ beta.carotene+allo.xanthin+ chl.a+ diato.xanthin+lutein+zea.xanthin+chl.a.total..a.allom., data=cfm, family=Gamma)
library(MASS)
mod2<-stepAIC(mod, k=log(nrow(fmc)), direction="both")
summary(mod2)
glm(div(fmc)~O2, data=cfm, family=gaussian)
mod<-glm(exp(div(fall))~O2+TOC, data=chem, family=Gamma)
summary(mod)

mod1<-glm(exp(div(f01mc))~O2*TOC, data=cfm, family=Gamma)
mod2<-glm(exp(div(f12mc))~O2*TOC, data=cfm, family=Gamma)

mod3<-glm(exp(div(f63mc))~O2*TOC, data=cfm, family=Gamma)
mod4<-glm(exp(div(f125mc))~O2*TOC, data=cfm, family=Gamma)

mod5<-glm(exp(div(ffosmc))~O2*TOC, data=cfm, family=Gamma)
mod6<-glm(exp(div(fnofosmc))~O2*TOC, data=cfm, family=Gamma)

summary(mod)
summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)
summary(mod6)








#ordinations 
decorana(log(fdeadmc+1))
mod.dead<-cca(log(fdeadmc+1)~O2,data=cfm[!rownames(cfm)%in%c("RC5","RC9"),], subset=!is.na(cfm$O2)[!rownames(cfm)%in%c("RC5","RC9")])
mod.fos<-cca(log(ffosmc+1)[!rownames(cfm)%in%c("RC5","RC9"),]~O2,data=cfm[!rownames(cfm)%in%c("RC5","RC9"),], subset=!is.na(cfm$O2)[!rownames(cfm)%in%c("RC5","RC9")])
mod.mac<-cca(log(mfc+1)[!rownames(cfm)%in%c("RC5","RC9"),]~O2,data=cfm[!rownames(cfm)%in%c("RC5","RC9"),], subset=!is.na(cfm$O2)[!rownames(cfm)%in%c("RC5","RC9")])
x11();plot(mod.dead, disp=c("sites","cn"))
x11();plot(mod.fos, disp=c("sites","cn"))
x11();plot(mod.mac, disp=c("sites","cn"))

p<-protest(mod.fos,mod.dead)
x11();
plot(p)
text(p)

mod.fos<-cca(log(ffosmc+1)[!rownames(cfm)%in%c("RC5","RC9"),colSums(ffosmc>0)>2]~O2+TOC,data=cfm[!rownames(cfm)%in%c("RC5","RC9"),], subset=!is.na(cfm$O2)[!rownames(cfm)%in%c("RC5","RC9")])
x11();plot(mod.fos)            
anova(mod.fos, by="margin")

                                                                         
mod.fos<-cca(sqrt(fmc/rowSums(fmc))[!rownames(cfm)%in%c("RC5","RC9"),colSums(fmc>0)>2]~O2,data=cfm[!rownames(cfm)%in%c("RC5","RC9"),], subset=!is.na(cfm$O2)[!rownames(cfm)%in%c("RC5","RC9")])
mod.fos2<-cca(sqrt(fmc/rowSums(fmc))[!rownames(cfm)%in%c("RC5","RC9"),colSums(fmc>0)>2]~TOC,data=cfm[!rownames(cfm)%in%c("RC5","RC9"),], subset=!is.na(cfm$O2)[!rownames(cfm)%in%c("RC5","RC9")])
x11();plot(mod.fos)
x11();plot(mod.fos2)
anova(mod.fos, by="margin")
anova(mod.fos2, by="margin")

mod.fos<-cca(sqrt(fmc/rowSums(fmc))[,colSums(fmc>0)>2]~exp(div(fmc)))
x11();plot(mod.fos)
anova(mod.fos)
   colSums(fmc>0)
   
strat.plot(100*(ffosmc/rowSums(ffosmc))[order(cfm$O2),colSums(ffosmc>0)>10], cfm$O2[order(cfm$O2)], wa.order="topleft", scale.percent=T)   
strat.plot(100*(ffosmc/rowSums(ffosmc))[order(cfm$TOC),colSums(ffosmc>0)>10], cfm$TOC[order(cfm$TOC)], wa.order="topleft", scale.percent=T)   

x11()
strat.plot(100*(fall/rowSums(fall))[order(chem$O2),colSums(fall>0)>10], chem$O2[order(chem$O2)], wa.order="topleft", scale.percent=T, ylabel=expression(O[2]), plot.line=F, lwd.bar=2, col.bar=1, cex.xlabel=.7)   
x11()
strat.plot(100*(fmc/rowSums(fmc))[order(cfm$TOC),colSums(fmc>0)>10], cfm$TOC[order(cfm$TOC)], wa.order="topleft", scale.percent=T,, ylabel="TOC", plot.line=F, lwd.bar=2, col.bar=1, cex.xlabel=.7, y.rev=T)   
x11()
strat.plot(100*(fall/rowSums(fall))[order(div(fall)),colSums(fall>0)>10], exp(div(fall))[order(div(fall))], wa.order="topleft", scale.percent=T, ylabel="exp(H')", plot.line=F, lwd.bar=2, col.bar=1, cex.xlabel=.7)   

x11()
strat.plot(100*(macro8g/rowSums(macro8g))[order(div(macro8g)),colSums(macro8g>0)>10], exp(div(macro8g))[order(div(macro8g))], wa.order="topleft", scale.percent=T, ylabel="exp(H')", plot.line=F, lwd.bar=2, col.bar=1, cex.xlabel=.7)   

wa<-WA(sqrt(fall/rowSums(fall))[!is.na(chem$mO2),colSums(fall>0)>4], chem$mO2[!is.na(chem$mO2)], tolDW = T)
opt<-wa$coef
opt<-opt[order(opt[,1]),]
x11();par(mar=c(3,8.5,1,.5), mgp=c(1.5,.5,0), tcl=-.4)
plot(opt[,1], 1:nrow(opt), yaxt="n", xlim=range(c(opt[,1]-opt[,2],opt[,1]+opt[,2])), pch=16, ylab="", ylim=c(0,nrow(opt)+1), yaxs="i", xlab="Effective number of species, exp(H')")
arrows(opt[,1]-opt[,2], 1:nrow(opt), x1 = opt[,1]+opt[,2],  code = 0)
axis(2,at=1:nrow(opt),rownames(opt), las=1, cex.axis=.65)

wa<-WA(sqrt(macro8g/rowSums(macro8g))[rowSums(macro8g>0)>4,colSums(macro8g>0)>5], exp(div(macro8g))[rowSums(macro8g>0)>4], tolDW = T)
opt<-wa$coef
opt<-opt[order(opt[,1]),]
x11();par(mar=c(3,8.5,1,.5), mgp=c(1.5,.5,0), tcl=-.4)
plot(opt[,1], 1:nrow(opt), yaxt="n", xlim=range(c(opt[,1]-opt[,2],opt[,1]+opt[,2])), pch=16, ylab="", ylim=c(0,nrow(opt)+1), yaxs="i", xlab="Effective number of species, exp(H')")
arrows(opt[,1]-opt[,2], 1:nrow(opt), x1 = opt[,1]+opt[,2],  code = 0)
axis(2,at=1:nrow(opt),rownames(opt), las=1, cex.axis=.65)




#In the 3rd part of the paper, I will compare values between foraminifera and macrofauna for diversity indices (both H' and ES100):

#- total living forams versus macrofauna

#- fossilizable forams versus macrofauna

#- dead forams versus macrofauna

#Correlation?

#Do forams and macrofauna agree to evaluate Ecological Quality Status of stations?

#       high     good  moderate    poor       bad
#H'     >3,8  3,0-3,8   1,9-3,0   0,9-1,9    <0,9
#ES100  >25     17-25     10-17      5-10     <5

#Borders between EcoQ have been established with macrofauna, I can see some disagreements between macrofauna and foraminifera, so maybe it worth trying to calibrate the borders for the forams.

library(rioja)
sp<-fmc/rowSums(fmc)
decorana(sqrt(sp[!is.na(cfm$O2),colSums(sp>0)>3]))
performance(mod<-crossval(WA(sqrt(sp[!is.na(cfm$O2),colSums(sp>0)>3]),cfm$O2[!is.na(cfm$O2)])))

plot(mod)

################
#how are AMBI classifications made
#how does BQI correlate with H'/ES100
#Does ISI have any utility
#fossil diversity biased high in low 02 environments because of 1) time averaging, 2) taphonoic processes
#Can anything be used to correct this bias in fossil diversity



