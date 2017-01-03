library(RODBC)
library(vegan)
library(Hmisc)
library(gdata)
setwd(  "\\\\nturt.uib.no\\gbsrt\\DATA\\pes\\")
library(entropy)

div<-function(x){apply(x,1,entropy.ChaoShen)}

con<-odbcConnectAccess("Pes_DB_28jan2011.mdb")
macroF<-sqlQuery(con, "select * from pes_DB_macrofauna_at_forams_stations where station_code='RC5' OR station_code='RC9'")
foramF<-sqlQuery(con, "select * from pes_db_foraminifera_species_data where station_code='RC5' OR station_code='RC9'")
names(foramF)<-make.names(names(foramF))
chemF<-sqlQuery(con, "select * from pes_DB_chemistry_data_flat where stas='RC5' OR stas='RC9'")
names(macroF)<-make.names(names(macroF))

close(con)

chemSed<-xtabs(Value~paste(chemF$STAS,chemF$DATE, chemF$Slice_numeric, chemF$REPLICATE, sep="|")+Chemical_species, data=chemF, exclude="O2", drop.unused.levels = FALSE)
chemWater<-chemF[chemF$Chemical_species=="O2",]

#macro
macro<-xtabs(Number.1~paste(macroF$Station_code,macroF$DAT, sep="|")+SpeciesMacrofauna, data=macroF)
decorana(log(macro[,colSums(macro>0)>1]+1))
x11();plot(rda(log(macro[,colSums(macro>0)>1]+1)))
round(exp(div(macro)),1)


#foram
foram<-xtabs(Number.1~paste(foramF$Station_code,foramF$DAT, foramF$Slice_numeric, sep=":")+SpeciesForam, data=foramF)


decorana(log(foram+1))
mod<-rda(log(foram+1))
mod<-rda(sqrt(foram/rowSums(foram))[,colSums(foram>0)>2])

screeplot(mod, bstick=T)
plot(mod)

foram.meta<-as.data.frame(matrix(unlist(strsplit(rownames(foram),":")), ncol=3, byrow=T))
names(foram.meta)<-c("site","date","depth")
foram.meta$site<-as.factor(foram.meta$site)
foram.meta$date<- as.Date(as.character(foram.meta$date), format="%Y%m%d")

x11()
plot(mod, type="n")
text(mod, disp="sp", col=4, cex=.7)
sapply(1:nlevels(foram.meta$site),function(n){
 N<-levels(foram.meta$site)[n]
  sapply(format(unique(foram.meta$date[foram.meta$site==N])),function(d){
    sel<-foram.meta$site==N&foram.meta$date==d
    points(mod, disp="site", select=sel, col=n, type="l")
    text(mod, disp="site", select=which.max(sel), labels=rep(format(as.Date(d), "%Y-%m"), nrow(foram)),col=n, cex=.7)
  })
})


x11()
plot(mod, type="n")
text(mod, disp="sp", col=4, cex=.7)
sapply(1:nlevels(foram.meta$site),function(n){
 N<-levels(foram.meta$site)[n]
  sapply(unique(foram.meta$depth[foram.meta$site==N]),function(d){
    sel<-foram.meta$site==N&foram.meta$depth==d
    points(mod, disp="site", select=sel, col=n, type="l")
    text(mod, disp="site", select=which.max(sel), labels=rep(d, nrow(foram)),col=n, cex=.7)
  })
})


#diversity               
rowSums(foram)
min(rowSums(foram))
rare<-cbind(foram.meta,rare=rarefy(foram,10), n=rowSums(foram), n1=exp(div(foram)), sh=div(foram))


x11()
plot(rare[rare$site=="RC5"&rare$depth==.25,c(2,4)], type="o", ylim=range(rare[,4]), main="Depth=0.25cm")
points(rare[rare$site=="RC9"&rare$depth==.25,c(2,4)], col=2,type="o")


x11();
plot(rare[rare$site=="RC5"&rare$depth==.75,c(2,4)], type="o", ylim=range(rare[,4]), main="Depth=0.75cm")
points(rare[rare$site=="RC9"&rare$depth==.75,c(2,4)], col=2,type="o")


x11();
plot(rare[rare$site=="RC5"&rare$depth==2.5,c(2,4)], type="o", ylim=range(rare[,4]), main="Depth=2.5cm")
points(rare[rare$site=="RC9"&rare$depth==2.5,c(2,4)], col=2,type="o")



x11()
plot(rare[rare$site=="RC5"&rare$depth==.25,c(2,5)], type="o", ylim=range(rare[,5]), main="Depth=0.25cm")
points(rare[rare$site=="RC9"&rare$depth==.25,c(2,5)], col=2,type="o")


x11();
plot(rare[rare$site=="RC5"&rare$depth==.75,c(2,5)], type="o", ylim=range(rare[,5]), main="Depth=0.75cm")
points(rare[rare$site=="RC9"&rare$depth==.75,c(2,5)], col=2,type="o")


x11();
plot(rare[rare$site=="RC5"&rare$depth==2.5,c(2,5)], type="o", ylim=range(rare[,5]), main="Depth=2.5cm")
points(rare[rare$site=="RC9"&rare$depth==2.5,c(2,5)], col=2,type="o")



x11()
plot(rare[rare$site=="RC5"&rare$depth==.25,c(2,6)], type="o", ylim=range(rare[,6]), main="Depth=0.25cm")
points(rare[rare$site=="RC9"&rare$depth==.25,c(2,6)], col=2,type="o")


x11();
plot(rare[rare$site=="RC5"&rare$depth==.75,c(2,6)], type="o", ylim=range(rare[,6]), main="Depth=0.75cm")
points(rare[rare$site=="RC9"&rare$depth==.75,c(2,6)], col=2,type="o")


x11();
plot(rare[rare$site=="RC5"&rare$depth==2.5,c(2,6)], type="o", ylim=range(rare[,6]), main="Depth=2.5cm")
points(rare[rare$site=="RC9"&rare$depth==2.5,c(2,6)], col=2,type="o")

###############################
#foram
foram2<-xtabs(Number.split~paste(foramF$Station_code,foramF$DAT,foramF$Slice_numeric,foramF$Replicate, sep=":")+CodeForam, data=foramF)

decorana(log(foram2+1))
mod<-rda(log(foram2+1))

screeplot(mod, bstick=T)
plot(mod)

foram.meta2<-as.data.frame(matrix(unlist(strsplit(rownames(foram2),":")), ncol=4, byrow=T))
names(foram.meta2)<-c("site","date","depth", "replicate")
foram.meta2$site<-as.factor(foram.meta2$site)
foram.meta2$date<- as.Date(as.character(foram.meta2$date), format="%Y%m%d")

x11()
plot(mod, type="n")
text(mod, disp="sp", col=4, cex=.7)
sapply(1:nlevels(foram.meta2$site),function(n){
 N<-levels(foram.meta2$site)[n]
  sapply(format(unique(foram.meta$date[foram.meta2$site==N])),function(d){
    sel<-foram.meta2$site==N&foram.meta2$date==d
    points(mod, disp="site", select=sel, col=n, type="l")
    text(mod, disp="site", select=which.max(sel), labels=rep(format(as.Date(d), "%Y-%m"), nrow(foram)),col=n, cex=.7)
  })
})


x11()
plot(mod, type="n")
text(mod, disp="sp", col=4, cex=.7)
sapply(1:nlevels(foram.meta$site),function(n){
 N<-levels(foram.meta$site)[n]
  sapply(unique(foram.meta$depth[foram.meta$site==N]),function(d){
    sel<-foram.meta$site==N&foram.meta$depth==d
    points(mod, disp="site", select=sel, col=n, type="l")
    text(mod, disp="site", select=which.max(sel), labels=rep(d, nrow(foram)),col=n, cex=.7)
  })
})


#diversity               
rowSums(foram2)
min(rowSums(foram2))
rare2<-cbind(foram.meta2,rare=rarefy(foram2,4), n=rowSums(foram2), n1=exp(div(foram2)))


x11()
plot(rare2[rare2$site=="RC5"&rare2$depth==.25,c(2,5)], type="o", ylim=range(rare2[,5]), main="Depth=0.25cm")
points(rare2[rare2$site=="RC9"&rare2$depth==.25,c(2,5)], col=2,type="o")


x11();
plot(rare2[rare2$site=="RC5"&rare2$depth==.75,c(2,5)], type="o", ylim=range(rare2[,5]), main="Depth=0.75cm")
points(rare2[rare2$site=="RC9"&rare2$depth==.75,c(2,5)], col=2,type="o")


x11();
plot(rare[rare$site=="RC5"&rare$depth==2.5,c(2,4)], type="o", ylim=range(rare[,4]), main="Depth=2.5cm")
points(rare[rare$site=="RC9"&rare$depth==2.5,c(2,4)], col=2,type="o")



x11()
plot(rare[rare$site=="RC5"&rare$depth==.25,c(2,5)], type="o", ylim=range(rare[,5]), main="Depth=0.25cm")
points(rare[rare$site=="RC9"&rare$depth==.25,c(2,5)], col=2,type="o")


x11();
plot(rare[rare$site=="RC5"&rare$depth==.75,c(2,5)], type="o", ylim=range(rare[,5]), main="Depth=0.75cm")
points(rare[rare$site=="RC9"&rare$depth==.75,c(2,5)], col=2,type="o")


x11();
plot(rare[rare$site=="RC5"&rare$depth==2.5,c(2,5)], type="o", ylim=range(rare[,5]), main="Depth=2.5cm")
points(rare[rare$site=="RC9"&rare$depth==2.5,c(2,5)], col=2,type="o")



x11()
plot(rare[rare$site=="RC5"&rare$depth==.25,c(2,6)], type="o", ylim=range(rare[,6]), main="Depth=0.25cm")
points(rare[rare$site=="RC9"&rare$depth==.25,c(2,6)], col=2,type="o")


x11();
plot(rare[rare$site=="RC5"&rare$depth==.75,c(2,6)], type="o", ylim=range(rare[,6]), main="Depth=0.75cm")
points(rare[rare$site=="RC9"&rare$depth==.75,c(2,6)], col=2,type="o")


x11();
plot(rare[rare$site=="RC5"&rare$depth==2.5,c(2,6)], type="o", ylim=range(rare[,6]), main="Depth=2.5cm")
points(rare[rare$site=="RC9"&rare$depth==2.5,c(2,6)], col=2,type="o")

