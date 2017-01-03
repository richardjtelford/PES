lapply(list(mO2=chem$mO2, TOC=chem$TOC,CN=chem$TOC/chem$TN, grain.size=chem$X..63),function(env){
  fdiv<-lapply(fff,function(f)exp(div(f)))
  cor.d1<-c(sapply(fdiv,cor,fdiv$all),deadvsfos=cor(fdiv$dead,fdiv$fos))
  lapply(fdiv,cor.test,env,use="complete")
  cor.o2<-c(sapply(fdiv,cor,env, use="complete"),  deadvsfos=NA)
  
  ford<-lapply(fff,function(f)rda(sqrt(f/rowSums(f))))
  m2<-c(sapply(ford,function(f)sqrt(1-procrustes(ford$all,f, sym=T)$ss)),deadvsfos=sqrt(1-procrustes(ford$fos,ford$dead, sym=T)$ss))
  
  of<-lapply(fff,function(f)rda(sqrt(f/rowSums(f))~env,  na.action = na.omit))
  o2c<-c(sapply(of,function(f)f$CCA$tot.chi/f$tot.chi),deadvsfos=NA)
  
  round(cbind(cor.diversity.env=cor.o2,cor.diversity.full.diversity=cor.d1, community.var.explained.by.env=o2c,procrustes.m2=m2),2)
})

