##drillcutting
library(vegan)
dat<-read.table("Drill_cutting_samples_for_Richard.txt", sep="\t", header=T)
meta<-dat[,1:5]
spp<-dat[,-c(1:5)]
rownames(spp)<-paste(meta[,4],c(0.5,1.5))                                                            
meta2<-meta[rowSums(spp)>0,]  
spp<-spp[rowSums(spp)>0,]                                                              
spp<-log(spp+1)
spp<-spp[,colSums(spp>0)>2]
#spp<-spp/rowSums(spp)

decorana(spp)
mod0<-cca((spp)~1, data=meta2,subset=meta2$Treatment!="Control"&meta2$Depth.in.sediment..cm.=="0-1")
mod1<-cca((spp)~Treatment, data=meta2, subset=meta2$Treatment!="Control"&meta2$Depth.in.sediment..cm.=="0-1")
anova(mod1)
mod2<-cca((spp)~Thickness, data=meta2,subset=meta2$Treatment!="Control"&meta2$Depth.in.sediment..cm.=="0-1")
anova(mod2)

mod2<-cca((spp)~Thickness+Condition(Treatment), data=meta2, subset=meta2$Treatment!="Control"&meta2$Depth.in.sediment..cm.=="0-1")
anova(mod2)

mod2<-cca((spp)~Treatment+Condition(Thickness), data=meta2, subset=meta2$Treatment!="Control"&meta2$Depth.in.sediment..cm.=="0-1")
anova(mod2)

mod3<-cca((spp)~Treatment+Thickness, data=meta2, subset=meta2$Treatment!="Control"&meta2$Depth.in.sediment..cm.=="0-1")
anova(mod3, by="marg")
plot(mod3)

mod3<-cca((spp)~Thickness+Treatment, data=meta2, subset=meta2$Treatment!="Control"&meta2$Depth.in.sediment..cm.==" 1-2")
anova(mod3, by="marg")


vp=varpart(sqrt(spp), ~Treatment+Thickness, data=meta2, subset=meta2$Treatment!="Control"&meta2$Depth.in.sediment..cm.=="0-1")
plot(vp)

plot(mod3)
                     
plot(mod0)
ordihull(mod0, meta$Treat, label=T)


##########diversity
sppr<-dat[,-c(1:5)]  
x11();par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(meta$Thickness[meta$Depth=="0-1"],exp(div(sppr))[meta$Depth=="0-1"], col=meta$Treat[meta$Depth=="0-1"])
plot(meta$Thickness[meta$Depth==" 1-2"],exp(div(sppr))[meta$Depth==" 1-2"], col=meta$Treat[meta$Depth==" 1-2"])


x11();par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(meta$Thickness[meta$Depth=="0-1"],rowSums(sppr)[meta$Depth=="0-1"], col=meta$Treat[meta$Depth=="0-1"], log="y")
plot(meta$Thickness[meta$Depth==" 1-2"],rowSums(sppr)[meta$Depth==" 1-2"], col=meta$Treat[meta$Depth==" 1-2"], log="y")


x11();par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(meta$Thickness[meta$Depth=="0-1"],sppr$Textularia.earlandi[meta$Depth=="0-1"], col=meta$Treat[meta$Depth=="0-1"], log="y")
plot(meta$Thickness[meta$Depth==" 1-2"],sppr$Textularia.earlandi[meta$Depth==" 1-2"], col=meta$Treat[meta$Depth==" 1-2"], log="y")
                                                                 
                                                                 
x11();par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(meta$Thickness[meta$Depth=="0-1"],sppr$Adercotryma.wrighti.glomeratum[meta$Depth=="0-1"], col=meta$Treat[meta$Depth=="0-1"], log="y")
plot(meta$Thickness[meta$Depth==" 1-2"],sppr$Adercotryma.wrighti.glomeratum[meta$Depth==" 1-2"], col=meta$Treat[meta$Depth==" 1-2"], log="y")


x11();par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(meta$Thickness[meta$Depth=="0-1"],sppr$Eggerelloides.scaber[meta$Depth=="0-1"], col=meta$Treat[meta$Depth=="0-1"], log="y")
plot(meta$Thickness[meta$Depth==" 1-2"],sppr$Eggerelloides.scaber[meta$Depth==" 1-2"], col=meta$Treat[meta$Depth==" 1-2"], log="y")

x11();par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(meta$Thickness[meta$Depth=="0-1"],sppr$Recurvoides.trochamminiforme[meta$Depth=="0-1"], col=meta$Treat[meta$Depth=="0-1"], log="y")
plot(meta$Thickness[meta$Depth==" 1-2"],sppr$Recurvoides.trochamminiforme[meta$Depth==" 1-2"], col=meta$Treat[meta$Depth==" 1-2"], log="y")




x11();par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(meta$Thickness[meta$Depth=="0-1"],sppr$Ammodiscus.gullmarensis[meta$Depth=="0-1"], col=meta$Treat[meta$Depth=="0-1"], log="y")
plot(meta$Thickness[meta$Depth==" 1-2"],sppr$Ammodiscus.gullmarensis[meta$Depth==" 1-2"], col=meta$Treat[meta$Depth==" 1-2"], log="y")

                        
x11();par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(meta$Thickness[meta$Depth=="0-1"],sppr$Textularia.kattegatensis[meta$Depth=="0-1"], col=meta$Treat[meta$Depth=="0-1"], log="y")
plot(meta$Thickness[meta$Depth==" 1-2"],sppr$Textularia.kattegatensis[meta$Depth==" 1-2"], col=meta$Treat[meta$Depth==" 1-2"], log="y")

                        
x11();par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(meta$Thickness[meta$Depth=="0-1"],sppr$Reophax.fusiformis[meta$Depth=="0-1"], col=meta$Treat[meta$Depth=="0-1"], log="y")
plot(meta$Thickness[meta$Depth==" 1-2"],sppr$Reophax.fusiformis[meta$Depth==" 1-2"], col=meta$Treat[meta$Depth==" 1-2"], log="y")

mod1<-glm(Eggerelloides.scaber~Thickness*Treatment, data=cbind(meta,sppr), family=quasipoisson, subset=Depth.in.sediment..cm.=="0-1"&Thickness>3)
summary(mod1)

mod2<-glm(Eggerelloides.scaber~Thickness*Treatment, data=cbind(meta,sppr), family=quasipoisson, subset=Depth.in.sediment..cm.==" 1-2"&Thickness>3)
summary(mod2)


mod1<-glm(Adercotryma.wrighti.glomeratum~Thickness*Treatment, data=cbind(meta,sppr), family=quasipoisson, subset=Depth.in.sediment..cm.=="0-1"&Thickness>3)
summary(mod1)

mod2<-glm(Adercotryma.wrighti.glomeratum~Thickness*Treatment, data=cbind(meta,sppr), family=quasipoisson, subset=Depth.in.sediment..cm.==" 1-2"&Thickness>3)
summary(mod2)

mod1<-glm(Textularia.earlandi~Thickness*Treatment, data=cbind(meta,sppr), family=quasipoisson, subset=Depth.in.sediment..cm.=="0-1"&Thickness>3)
summary(mod1)
mod2<-glm(Textularia.earlandi~Thickness*Treatment, data=cbind(meta,sppr), family=quasipoisson, subset=Depth.in.sediment..cm.==" 1-2"&Thickness>3)
summary(mod2)


##########################################
dat2<-read.table("Borkaks_all_spp_31Jan_2011.txt", sep="\t", header=T)
spp2<-dat2[,-(1:5)]

sort(colSums(spp2, na.rm=T))
sort(colSums(spp2>0, na.rm=T))

x11();par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(meta$Thickness[meta$Depth=="0-1"],spp2$Stainforthia.fusiformis[meta$Depth=="0-1"], col=meta$Treat[meta$Depth=="0-1"], log="y")
plot(meta$Thickness[meta$Depth==" 1-2"],spp2$Stainforthia.fusiformis[meta$Depth==" 1-2"], col=meta$Treat[meta$Depth==" 1-2"], log="y")

mod1<-glm(Stainforthia.fusiformis~Thickness*Treatment, data=cbind(meta,spp2), family=quasipoisson, subset=Depth.in.sediment..cm.=="0-1"&Thickness>3)
summary(mod1)
mod2<-glm(Stainforthia.fusiformis~Thickness*Treatment, data=cbind(meta,spp2), family=quasipoisson, subset=Depth.in.sediment..cm.==" 1-2"&Thickness>3)
summary(mod2)


x11();par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(meta$Thickness[meta$Depth=="0-1"],spp2$Bulimina.marginata[meta$Depth=="0-1"], col=meta$Treat[meta$Depth=="0-1"], log="y")
plot(meta$Thickness[meta$Depth==" 1-2"],spp2$Bulimina.marginata[meta$Depth==" 1-2"], col=meta$Treat[meta$Depth==" 1-2"], log="y")

mod1<-glm(Bulimina.marginata~Thickness*Treatment, data=cbind(meta,spp2), family=quasipoisson, subset=Depth.in.sediment..cm.=="0-1"&Thickness>3)
summary(mod1)
mod2<-glm(Bulimina.marginata~Thickness*Treatment, data=cbind(meta,spp2), family=quasipoisson, subset=Depth.in.sediment..cm.==" 1-2"&Thickness>3)
summary(mod2)

x11();par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(meta$Thickness[meta$Depth=="0-1"],spp2$Nonionellina.labradorica[meta$Depth=="0-1"], col=meta$Treat[meta$Depth=="0-1"], log="y")
plot(meta$Thickness[meta$Depth==" 1-2"],spp2$Nonionellina.labradorica[meta$Depth==" 1-2"], col=meta$Treat[meta$Depth==" 1-2"], log="y")

mod1<-glm(Nonionellina.labradorica~Thickness*Treatment, data=cbind(meta,spp2), family=quasipoisson, subset=Depth.in.sediment..cm.=="0-1"&Thickness>3)
summary(mod1)
mod2<-glm(Nonionellina.labradorica~Thickness*Treatment, data=cbind(meta,spp2), family=quasipoisson, subset=Depth.in.sediment..cm.==" 1-2"&Thickness>3)
summary(mod2)



Nonionellina.labradorica