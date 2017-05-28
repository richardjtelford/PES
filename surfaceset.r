#### Analyses of surface data

## load packages
library("vegan")
library("cocorresp")
library("ggfortify")
devtools::install_github("richardjtelford/ggbiplot", ref = "experimental")#check latest version of ggbiplot
library("ggbiplot")
library("GGally")
library("ggvegan")
library("ggrepel")

## load data
source("load.vincent.data.r")

#load functions
source("fortify.envfit.R")


#histogram of chemistry
ggplot(chem0, aes(x = Val)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~ Chemical_species, scales = "free_x")

#PCA of chemistry
chem.pca <- rda(select(chem_complete, -Station_code), scale = TRUE)
screeplot(chem.pca, bstick=TRUE)

f.pca <- fortify(chem.pca, scaling = "symmetric")
f.pca.species <- f.pca %>% filter(Score == "species") %>%
  mutate(Label = recode(Label, "O2" = "O[2]", "%<63" = "'%'~'<'~63~mu*m", "Depth Below Threshold" = "Depth~Below~Threshold"))
f.pca.sites <- f.pca %>% 
  filter(Score == "sites") %>% 
  mutate(Label = chem_complete$Station_code)

propVar <- eigenvals(chem.pca) / sum(eigenvals(chem.pca)) * 100
propVar <- format(propVar, digits = 1)

ggplot(f.pca.sites, aes(x = Dim1, y = Dim2, label = Label)) + 
  geom_point() +
  geom_text_repel(size = 3.5) +
  geom_axis(data = f.pca.species, parse = TRUE) +
  coord_equal() +
  scale_x_continuous(expand = c(0.3, 0.3), minor_breaks = .Machine$double.eps) + 
  scale_y_continuous(expand = c(0.3, 0.3), minor_breaks = .Machine$double.eps) +
  labs(x = paste0("PC1 (",propVar[1],"%)"), y = paste0("PC2 (",propVar[2],"%)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank())
  
ggsave("Chemistry_PCA.eps", width = 90, height = 90, units = "mm")#save plot
#autoplot(chem.pca, scaling = "symmetric", legend.position = "none")
 

#pairs plots of chemistry
setNames(chem, make.names(names(chem))) %>% 
  ggpairs(columns = 2:ncol(chem))

#pigments
pig <- names(chem) %in% pigments
RDA <- rda(chem[, pig] ~ ., data = select(chem, -which(pig), -Station_code), scale = TRUE)
RDA0 <- rda(chem[, pig] ~ 1, data = select(chem, -which(pig), -Station_code), scale = TRUE)
plot(RDA0)
ordisurf(RDA0, chem$O2, add = TRUE)

plot(chem$O2, scores(RDA0, disp = "sites")[, 1])
scatter.smooth(chem$O2, scores(RDA0, disp = "sites")[, 1])
identify(chem$O2, scores(RDA0, disp = "sites")[, 1], chem$Station_code)

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

#procrustes analysis
foram.cca <- cca(foram8m30 ~ 1)
macro.cca <- cca(macro8f30 ~ 1)

pt <- protest(foram.cca, macro.cca, permutations = 999)
pt
plot(pt)


 #co=correspondance

plot(vegdist(macro8gf[rowSums(macro8gf)>10,]),vegdist(foram8m[rowSums(macro8gf)>10,]))
mantel(vegdist(macro8gf),vegdist(foram8m))

plot(stepacross(vegdist(macro8gf), toolong=.9),stepacross(vegdist(foram8m), toolong=.9))
mantel(stepacross(vegdist(macro8gf), toolong=.9),stepacross(vegdist(foram8m), toolong=.9))

## predictive CoCA using SIMPLS and formula interface

coco.pred <- coca(macro8f30 ~ ., data = foram8m30)
coco.pred
plot(coco.pred)
## should retain only the useful PLS components for a parsimonious model

## Leave-one-out crossvalidation - this takes a while
cv <- crossval(macro8f30,  foram8m30)
data_frame(cv = cv$CVfit, x = 1:length(cv)) %>% 
  ggplot(aes(x = x, y = cv)) +
    geom_point() +
    labs(x = "Number of axes", y = "Cross-validatory %fit")
## so 2 axes are sufficient
## permutation test to assess significant PLS components - takes a while
bp.perm <- permutest(coco.pred, permutations = 999, n.axes = 5)
bp.perm
as_data_frame(unclass(bp.perm)[1:7]) %>% 
  mutate(x = 1:n()) %>%
  ggplot(aes(x = x, y = pcent.fit, colour = pval < 0.05 )) + 
  geom_point()

summary(bp.perm)

coco.pred <- coca(macro8f30 ~ ., data = foram8m30 , n.axes = 2)
coco.pred 

#envfit
sol.pred <- envfit(coco.pred, env = select(chem30, O2, TN, TOC), which = "predictor")
sol.resp <- envfit(coco.pred, env = select(chem30, O2, TN, TOC), which = "response")


par(mfrow = c(1, 2))
plot(coco.pred, which = "predictor")
plot(sol.pred)
plot(coco.pred, which = "response")
plot(sol.resp)

# sco <- scores(coco.pred)
# sco$sites <- bind_rows(
#   sco$sites$X %>% as_data_frame() %>% mutate(which = "predictor", site = Station_code30),
#   sco$sites$Y %>% as_data_frame() %>% mutate(which = "response", site = Station_code30)
#   ) 
# sco$species <- bind_rows(
#   sco$species$X %>% as_data_frame() %>% mutate(which = "predictor"),
#   sco$species$Y %>% as_data_frame() %>% mutate(which = "response")
# ) 

#plot function for coco-predict
autoplot.coco <- function(x, which = c("response", "predictor"), th = theme_bw(), envfit = NULL){
  which <- match.arg(which)
  WHICH <- ifelse(which == "response", "Y", "X")
  sco <- scores(x)
  sco <- lapply(sco, `[[`, WHICH)
  
  g <- ggplot(sco$sites, aes(x = `Comp 1`, y = `Comp 2`)) + 
    geom_point() +
    geom_point(data = sco$species, colour = "red", shape = 2)+
    coord_equal() +
    scale_x_continuous(minor_breaks = .Machine$double.eps * 10) +
    scale_y_continuous(minor_breaks = .Machine$double.eps * 10) +
    th +
    theme(panel.grid.major = element_blank()) +
    labs(x = "CoCA axis 1", y = "CoCA axis 2")

  if(!is.null(envfit)){
    envfit <- fortify(envfit)

    #arrow.mul
    fill = 0.75 
    u <- c(range(c(sco$sites[,'Comp 1'], sco$species[,'Comp 1'])), range(c(sco$sites[,'Comp 2'], sco$species[,'Comp 2'])))
    r <- c(range(envfit[, 1], na.rm = TRUE), range(envfit[, 2], na.rm = TRUE))
    u <- u/r
    u <- u[is.finite(u) & u > 0]
    arrow.mul <- fill * min(u)

    g <- g +  geom_axis(data = envfit, mapping = aes(x = Comp.1 * arrow.mul, y = Comp.2 * arrow.mul, label = labs), parse = TRUE)
  }
  
  g
}

a <- autoplot.coco(coco.pred, which = "predictor", envfit = sol.pred) + ggtitle("Predictor - forams")
b <- autoplot.coco(coco.pred, which = "response", envfit = sol.resp) + ggtitle("Response - macrofauna")
gridExtra::grid.arrange(a, b, nrow = 1)


plot(diversity(macro8gf),diversity(foram8m))
abline(0,1)