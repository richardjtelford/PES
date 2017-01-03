library(RODBC)
library(vegan)
library(Hmisc)
library(gdata)
library(entropy)
library(maps)
library(mapdata)

setwd(  "\\\\nturt.uib.no\\gbsrt\\DATA\\pes\\")
div<-function(x){apply(x,1,entropy.ChaoShen)}

con<-odbcConnectAccess("Pes_db_8Feb2011.mdb")

stations<-sqlQuery(con,"select * from pes_db_stations")
names(stations)

macro<-sqlQuery(con, "select * from pes_db_macrofauna_at_forams_stations")
names(macro)<-make.names(names(macro))

forams<-sqlQuery(con,"select * from PES_DB_foraminifera_species_data where slice_numeric>0 and not Size  = '>500' and dat>20050000 and slice_numeric<3")
names(forams)<-make.names(names(forams))
dim(forams)

dead<-sqlQuery(con,"select * from PES_DB_foraminifera_species_dead_data")
names(dead)<-make.names(names(dead))
dim(dead)



chem0<-sqlQuery(con,"SELECT STAS, Chemical_species, Avg(IIf([Value]<-98,0,[value])) AS Val FROM PES_db_chemistry_data_flat WHERE (((DATE)>20050000) AND ((Slice_numeric)<1 Or (Slice_numeric) Is Null)) GROUP BY STAS, Chemical_species;")

fspp<-sqlQuery(con,"select * from PES_DB_foraminifera_species_list")
names(fspp)<-make.names(names(fspp))


close(con)

###map
map("worldHires",xlim=c(7.5,11.5), ylim=c(58,60))
points(stations$"X_COORD_DECIMALDEGREE",stations$"Y_COORD_DECIMALDEGREE", col=2)
box()


#chemistry data
names(chem0)
levels(chem0$STAS)
levels(chem0$STAS)<-trim(levels(chem0$STAS))
levels(chem0$STAS)

table(chem0$Chem)

wantstation<-chem0$STAS%in%levels(macro$Station)|chem0$STAS%in%levels(forams$Station_code)
unique(chem0$STAS[!wantstation])
sum(wantstation)
drop.chem<-c("Zn","Cu","Pb", "Cd","H2S", "Chl c3")
chem1<-data.frame(unclass(xtabs(rep(1,nrow(chem0))~STAS+Chemical_species, data=chem0, subset=wantstation&!chem0$Chem%in%drop.chem, drop=T)))
range(chem1)
chem<-data.frame(unclass(xtabs(Val~STAS+Chemical_species, data=chem0, subset=wantstation&!chem0$Chem%in%drop.chem, drop=T)))

chem[chem1==0]<-NA
dim(chem)
apply(chem,2,function(x)sum(is.na(x)))
apply(chem,1,function(x)sum(is.na(x)))
lapply(chem[,colSums(is.na(chem))>0],function(x)rownames(chem)[is.na(x)])

chem$Chl.a.allomer<-NULL
chem$Pheophorbide<-NULL

#chem<-chem[rownames(chem)!="KRG",]

sum(is.na(chem))  

#standardise pigments by TOC
names(chem)         
pig<- names(chem)%in%c("allo.xanthin","aphanizophyll","beta.carotene","cantha.xanthin","chl.a","chl.a.total..a.allom.","Chl.b","Chl.c2","diadino.xanthin","diato.xanthin","fuco.xanthin","lutein","peridinin","pheo.phytin.a","Pheophytin.b","Pyropheophytin.b", "violaxanthin","zea.xanthin")

chem[,pig]<-chem[,pig]/chem$TOC

plot(chem$MinO2,chem$O2)

chem$mO2<-chem$MinO2_2_years
chem$mO2[is.na(chem$mO2)]<-chem$O2[is.na(chem$mO2)]


#macros

macro3r<-xtabs(Number.1~paste(Station_code,Replicate, sep=":")+CodeMacrofauna, data=macro[macro$DAT<20050000,])
macro8r<-xtabs(Number.1~paste(Station_code,Replicate, sep=":")+CodeMacrofauna, data=macro[macro$DAT>20050000,])

macro3g<-xtabs(Number.1~Station_code+CodeMacrofauna, data=macro[macro$DAT<20050000,], drop=T)
macro8g<-xtabs(Number.1~Station_code+CodeMacrofauna, data=macro[macro$DAT>20050000,], drop=T)

rownames(macro3g)
rownames(macro8g)


rownames(macro3g)[!rownames(macro3g)%in%rownames(chem)]
rownames(macro8g)[!rownames(macro8g)%in%rownames(chem)]

rownames(chem)[!rownames(chem)%in%rownames(macro8g)]

##########
#forams
foram3r<-xtabs(Number.1~paste(Station_code,Replicate, sep=":")+SpeciesForam, data=forams[forams$DAT<20050000,])
foram8r<-xtabs(Number.1~paste(Station_code,Replicate, sep=":")+SpeciesForam, data=forams[forams$DAT>20050000,])

dim(foram8r)
dim(macro8r)

foram3<-xtabs(Number.1~Station_code+SpeciesForam, data=forams[forams$DAT<20050000,], drop=T)
foram8<-xtabs(Number.1~Station_code+SpeciesForam, data=forams[forams$DAT>20050000,], drop=T)


