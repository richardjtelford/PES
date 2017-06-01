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
  mutate(Label = recode(Label, "O2" = "O[2]", "pc_lt_63" = "'%'~'<'~63~mu*m", "DBT" = "Depth~Below~Threshold"))
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
pig <- names(chem_complete) %in% pigments
RDA <- rda(chem_complete[, pig] ~ ., data = select(chem_complete, -which(pig), -Station_code), scale = TRUE)
RDA0 <- rda(chem_complete[, pig] ~ 1, data = select(chem_complete, -which(pig), -Station_code), scale = TRUE)
plot(RDA0)
ordisurf(RDA0, chem_complete$O2, add = TRUE)

plot(chem_complete$O2, scores(RDA0, disp = "sites")[, 1])
scatter.smooth(chem_complete$O2, scores(RDA0, disp = "sites")[, 1])
identify(chem_complete$O2, scores(RDA0, disp = "sites")[, 1], chem_complete$Station_code)

RDA1 <- ordistep(RDA0, reformulate(names(select(chem_complete, -which(pig), -Station_code))), permutations = how(nperm = 999))
RDA1
plot(RDA1)
anova(RDA1, by = "margin")
RDAo <- rda(chem_complete[, pig] ~ O2, data = chem_complete[, !pig], scale = TRUE)
RDAo
plot(RDAo)
RDAoc <- rda(chem_complete[, pig] ~ O2 + TN, data = chem_complete[, !pig], scale = TRUE)
RDAoc

#macros
fortify_CCA <- function(x){
  x <- fortify(x, scaling = "symmetric")
  species <- x %>% filter(Score == "species")
  sites <- x %>% 
    filter(Score == "sites") %>% 
    mutate(Label = Station_code30)
  list(sites = sites, species = species)
}

plot.CCA <- function(x, xlab = "CA1", ylab = "CA2", exp = 0.2){
  ggplot(x$sites, aes(x = Dim1, y = Dim2, label = Label)) + 
    geom_point() +
    geom_text_repel(size = 3.5) +
    geom_point(data = x$species, shape = 3, colour = "red") +
    coord_equal() +
    scale_x_continuous(expand = c(exp, 0), minor_breaks = .Machine$double.eps) + 
    scale_y_continuous(expand = c(exp, 0), minor_breaks = .Machine$double.eps) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    theme(panel.grid.major = element_blank())
}

ggpairs(select(chem30, -one_of(pigments)))


#unconstrained macro
decorana(macro8f30)#long gradient

macro.ca <- cca(macro8f30 ~ 1, data = select(chem30, -one_of(pigments)))
screeplot(macro.ca, bstick = TRUE)
plot(macro.ca)

fmacro <- fortify_CCA(macro.ca)
plot.CCA(fmacro)

#unconstrained forams
ftaxa <- c("Cassidulina laevigata", "Liebusella. goësi", "Micrometula. hyalostriata", "Phainogulmia. aurata", "Textularia. earlandi" and "Recurvoides. trochamminiforme", "Bulimina marginata", "Cribrostomoides. bertheloti", "Cylindrogulmia. alba", "Leptohalysis. scottii" and "Spiroplectamina. biformis", "Stainforthia fusiformis", "Fissurina. sp.", "Bolivina. pseudopunctata") 

decorana(foram8m30)#shortish gradient

foram.ca <- cca(foram8m30 ~ 1, data = select(chem30, -one_of(pigments)))
screeplot(foram.ca, bstick = TRUE)
plot(foram.ca)

fforam <- fortify_CCA(foram.ca)
plot.CCA(fforam)



macro.mod <- ordistep(macro.ca, reformulate(names(select(chem30, -one_of(pigments)))))
macro.mod

plot(macro.mod)


foram.mod <- ordistep(foram.ca, reformulate(names(select(chem30, -one_of(pigments)))))
foram.mod

plot(foram.mod)



##########
#forams
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


#procrustes analysis
foram.cca <- cca(foram8m30 ~ 1)
macro.cca <- cca(macro8f30 ~ 1)

pt <- protest(foram.cca, macro.cca, permutations = 999)
pt
plot(pt)


#Mantel tests
dist_macro <- vegdist(macro8f30)
dist_foram <- vegdist(foram8m30)
ggplot(data_frame(macro = as.vector(dist_macro), 
                  foram = as.vector(dist_foram)), 
       aes(x = macro, y = foram)) +
  geom_point() +
  geom_smooth()

mantel(dist_macro, dist_foram)


#co=correspondance
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