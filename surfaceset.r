#### Analyses of surface data

## load packages
library("vegan")
library("cocorresp")

## load data
source("load.vincent.data.r")




#histogram of chemistry
ggplot(chem0, aes(x = Val)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~ Chemical_species, scales = "free_x")

#PCA of chemistry
chem.pca <- prcomp(chem, scale. = TRUE)


x11()
plot(chem.pca)
screeplot(chem.pca, bstick=TRUE)
biplot(chem.pca)
autoplot(chem.pca)

pairs(chem, col = ifelse(rownames(chem) == "GRO50", 2, 1), gap = 0)
plot(chem$TOC, chem$X)
scatter.smooth(chem$TOC, chem$X)


#pigments
RDA <- rda(chem[, pig] ~ ., data = chem[, !pig], scale = TRUE)
RDA0 <- rda(chem[, pig] ~ 1, data = chem[, !pig], scale = TRUE)
plot(RDA0)
ordisurf(RDA0, chem$O2, add = TRUE)

plot(chem$O2, scores(RDA0, disp = "sites")[, 1])
scatter.smooth(chem$O2, scores(RDA0, disp = "sites")[, 1])
identify(chem$O2, scores(RDA0, disp = "sites")[, 1], rownames(chem))

RDA1 <- step(RDA0, reformulate(names(chem[, !pig])), test = "perm")
RDA1
plot(RDA1)
RDAo <- rda(chem[, pig] ~ O2, data = chem[, !pig], scale = TRUE)
RDAo
plot(RDAo)
RDAoc <- rda(chem[, pig] ~ O2 + TOC, data = chem[, !pig], scale = TRUE)
RDAoc

#macros

x11();
boxplot(cbind(rowSums(m38),rowSums(m83)))#did the volume sampled remain constant?  --YES, Knockout in 2008
boxplot(cbind(rowSums(m38>0),rowSums(m83>0)))#species
x11();
plot(cbind(rowSums(m38),rowSums(m83)), xlab="2003 no ind", ylab="2008 no ind")#species
abline(0,1)
identify(cbind(rowSums(m38),rowSums(m83)), labels=rownames(m83))
x11();
plot(cbind(rowSums(m38>0),rowSums(m83>0)), xlab="2003 no spp", ylab="2008 no spp")#species
abline(0,1)
identify(cbind(rowSums(m38>0),rowSums(m83>0)), labels=rownames(m83))


plot(sort(rowSums(macro8), decreasing=TRUE), log="y")
x11()
plot(log(colSums(macro8g>0)),log(apply(macro8g,2,function(b)mean(b[b>0]))))#range abundance
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
x11()
plot(mod)

cor(diversity(macro8gc),chemM$O2)


##########
#forams
quantile(rowSums(foram8r))
quantile(rowSums(macro8r))
div.f<-tapply(diversity(foram8r),unlist(strsplit(rownames(foram8r),":"))[c(T,F)],sd)
div.m<-tapply(diversity(macro8r),unlist(strsplit(rownames(macro8r),":"))[c(T,F)],function(x)sd(x[1:min(3,length(x))]))
div.m <- div.m[names(div.m) %in% names(div.f)]
div.f <- div.f[names(div.f) %in% names(div.m)]
identical(names(div.f), names(div.m))


boxplot(list(foram=div.f,macro=div.m), notch=TRUE)
wilcox.test(div.f,div.m, paired=TRUE)             #try paired analysis?

div.f<-tapply(rarefy(foram8r,10),unlist(strsplit(rownames(foram8r),":"))[c(T,F)],sd)
div.m<-tapply(rarefy(macro8r,10),unlist(strsplit(rownames(macro8r),":"))[c(T,F)],function(x)sd(x[1:min(3,length(x))]))
div.m <- div.m[names(div.m) %in% names(div.f)]
div.f <- div.f[names(div.f) %in% names(div.m)]
identical(names(div.f), names(div.m))

boxplot(list(foram = div.f, macro = div.m))
wilcox.test(div.f, div.m, paired = TRUE)



plot(rowSums(foram38),rowSums( foram83), xlab="2003 no individuals",ylab="2008 no individuals")
abline(0,1)
text(rowSums(foram38),rowSums( foram83), labels=rownames(foram83))


x11()
plot(seq(10,rowSums(foram8)[1]/8,10),sapply(seq(10,rowSums(foram8)[1]/8,10),function(n)rarefy(foram8[1,,drop=FALSE]/8,n)))


mod<-cca(log(foram8[rowSums(foram8)>0,colSums(foram8>0)>1]+1))
screeplot(mod, bstick=TRUE)
plot(mod)
chemf<-chem[rownames(chem)%in%rownames(foram8),]
foram8c<-foram8[rownames(foram8)%in%rownames(chem),]
identical(rownames(foram8c), rownames(chemf))

cor(diversity(foram8c),chemf$O2)

modf<-step(cca(log(foram8c[,colSums(foram8c>0)>1]+1)~1, data=chemf), reformulate(names(chemf)), test="perm")
modf
plot(modf)





plot(vegdist(macro8gf[rowSums(macro8gf)>10,]),vegdist(foram8m[rowSums(macro8gf)>10,]))
mantel(vegdist(macro8gf),vegdist(foram8m))

plot(sort(colSums(foram8), dec=TRUE), log="y")
x11();
plot(sort(colSums(macro8), dec=TRUE), log="y")

 #co=correspondance

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

coco.pred <- coca(macro8gf ~ ., data = as.data.frame(foram8m) , n.axes = 2)
coco.pred 


plot(diversity(macro8gf),diversity(foram8m))
abline(0,1)