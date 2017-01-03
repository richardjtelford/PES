#"chem"     "chemf"          "chemM"      "foram8"      "foram8c"       "macro8g"     "macro8gc"   "mod.data"    "newdepths"        "stations"   

dim(chemM)
dim(chemf)
dim(foram8c)
dim(macro8gc)

colnames(chemM)
rm(chemf)
foram8c<-foram8c[rownames(foram8c)%in%rownames(macro8gc),]

foram<-foram8c
macro<-macro8gc
chem<-chemM
rm(foram8c,macro8gc, chemM)
mean(rownames(foram)==rownames(macro)&rownames(macro)==rownames(chem))

#n.ind
x11();par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(1.5,.5,0))
plot(chem$O2, rowSums(macro), xlab="O2", ylab="No. macrofauna")
plot(chem$O2, rowSums(foram), xlab="O2", ylab="No. forams")
plot(chem$TOC, rowSums(macro), xlab="TOC", ylab="No. macrofauna")
plot(chem$TOC, rowSums(foram), xlab="TOC", ylab="No. forams")

cor.test(chem$O2, rowSums(macro))
cor.test(chem$O2, rowSums(foram))
cor.test(chem$TOC, rowSums(macro))
cor.test(chem$TOC, rowSums(foram))


#n.spp
x11();par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(1.5,.5,0))
plot(chem$O2, rarefy(macro,100), xlab="O2", ylab="macro ES100", col=ifelse(rowSums(macro)<100,2,1))
plot(chem$O2, rarefy(foram,100), xlab="O2", , ylab="Foram ES100", col=ifelse(rowSums(foram)<100,2,1))
plot(chem$TOC, rarefy(macro,100), xlab="TOC", ylab="Macro ES100", col=ifelse(rowSums(macro)<100,2,1))
plot(chem$TOC, rarefy(foram,100), xlab="TOC",ylab="Foram ES100", col=ifelse(rowSums(foram)<100,2,1))

cor.test(chem$O2, rarefy(macro,100))
cor.test(chem$O2, rarefy(foram,100))
cor.test(chem$TOC, rarefy(macro,100))
cor.test(chem$TOC, rarefy(foram,100))

round(data.frame(forams=c(O2=cor(chem$O2, rarefy(foram,100)),TOC=cor(chem$TOC, rarefy(foram,100))),macros=c(cor(chem$O2, rarefy(macro,100)),cor(chem$TOC, rarefy(macro,100)))),2)


#diversity
x11();par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(1.5,.5,0))
plot(chem$O2, diversity(macro), xlab="O2", ylab="Macro Shannon", col=ifelse(rowSums(macro)<100,2,1))
plot(chem$O2, diversity(foram), xlab="O2", ylab="Foram Shannon", col=ifelse(rowSums(foram)<100,2,1))
plot(chem$TOC, diversity(macro), xlab="TOC", ylab="Macro Shannon", col=ifelse(rowSums(macro)<100,2,1))
plot(chem$TOC, diversity(foram), xlab="TOC", ylab="Foram Shannon", col=ifelse(rowSums(foram)<100,2,1))

cor.test(chem$O2, diversity(macro))
cor.test(chem$O2, diversity(foram))
cor.test(chem$TOC, diversity(macro))
cor.test(chem$TOC, diversity(foram))

round(data.frame(forams=c(O2=cor(chem$O2, diversity(foram)),TOC=cor(chem$TOC, diversity(foram))),macros=c(cor(chem$O2, diversity(macro)),cor(chem$TOC, diversity(macro)))),2)

library(rioja)
decorana(log(foram+1))
decorana(log(macro+1))

f02<-cca(log(foram+1)~O2, data=chem)
ftoc<-cca(log(foram+1)~TOC, data=chem)
m02<-cca(log(macro+1)~O2, data=chem)
mtoc<-cca(log(macro+1)~TOC, data=chem)

round(data.frame(forams=c(O2=f02$CCA$tot.chi/f02$tot.chi*100,TOC=ftoc$CCA$tot.chi/ftoc$tot.chi*100),macros=c(m02$CCA$tot.chi/m02$tot.chi*100,mtoc$CCA$tot.chi/mtoc$tot.chi*100)),2)


waO2f<-WA(sqrt(decostand(foram[,colSums(foram>0)>2], "total")),chem$O2)
watocf<-WA(sqrt(decostand(foram[,colSums(foram>0)>2], "total")),chem$TOC)
fo<-round(cbind(foram.O2=coef(waO2f), foram.TOC=coef(watocf)),2)
apply(fo,2,range)

waO2m<-WA(sqrt(decostand(macro[,colSums(macro>0)>2], "total")),chem$O2)
watocm<-WA(sqrt(decostand(macro[,colSums(macro>0)>2], "total")),chem$TOC)

mo<-round(cbind(macro.O2=coef(waO2m), macro.TOC=coef(watocm)),2)
apply(mo,2,range)

tf.perf<-rbind(
performance(crossval(waO2f))$crossval[1,,drop=F],
performance(crossval(waO2m))$crossval[1,,drop=F],
performance(crossval(watocf))$crossval[1,,drop=F],
performance(crossval(watocm))$crossval[1,,drop=F]
)
rownames(tf.perf)<-c("foram.O2", "macro.O2", "foram.toc", "macro.toc")
round(tf.perf,2)
sd(chem$TOC)
sd(chem$O2)