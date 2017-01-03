data.dir<-"d:/downloaded_data/pes/"
cores<-paste("core", c("21-91","4-89","9","30-91B","48-91","43-91","20-91B"))

core.data<-lapply(cores,function(core){
  print(core)
  tab<-read.table(paste(data.dir,core,".txt",sep=""), header=T, row.names=1, sep="\t", skip=1)
  tab<-tab[order(rownames(tab)),]
  colnames(tab)<-sub("X","",sub("-","_",paste(core,colnames(tab), sep="_"), ex=F))
  rownames(tab)<-sub("*","", rownames(tab), extended=F)
  tab<-t(tab)
  as.data.frame(tab)
  
})


sapply(core.data,ncol)

identical(sort(colnames(core.data[[1]])),sort(colnames(core.data[[2]])))

cbind(colnames(core.data[[1]]),colnames(core.data[[2]]))
colnames(core.data[[1]])[!colnames(core.data[[1]])%in%colnames(core.data[[2]])]
colnames(core.data[[2]])[!colnames(core.data[[2]])%in%colnames(core.data[[1]])]

cn<-unlist(lapply(core.data,colnames))
length(unique(cn))
table(cn)[table(cn)<length(core.data)]
length(.Last.value)

library(analogue)
all.cores<-join(core.data[[1]],core.data[[2]],core.data[[3]],core.data[[4]],core.data[[5]],core.data[[6]],core.data[[7]], split=F)
core.data<-join(core.data[[1]],core.data[[2]],core.data[[3]],core.data[[4]],core.data[[5]],core.data[[6]],core.data[[7]], split=T)
rownames(all.cores)
cr<-rep(cores,sapply(core.data,nrow))
plot(sort(rowSums(all.cores)))

rare<-lapply(core.data,function(cd){
  depth<-as.numeric(sapply(strsplit(rownames(cd),"_"),function(b)b[length(b)]))
  n=rowSums(cd)
  rare<-rarefy(cd,50)
  rare[n<50]<-NA
  cbind(depth,rare)
})
x11();par(mfrow=c(7,1), mar=c(3,4,1,1))
for(i in 1:length(core.data)){
  plot(rare[[i]], type="b", main=cores[i])
}

decorana(sqrt((all.cores/rowSums(all.cores))[rowSums(all.cores)>50,]))

ca.all<-cca(sqrt((all.cores/rowSums(all.cores))[rowSums(all.cores)>50,]))

plot(ca.all, type="n")
for(i in 1:length(cores)){
  points(ca.all, select=cr[rowSums(all.cores)>50]==cores[i], type="o", col=i, pch=20)
  points(ca.all, select=which.max(cr[rowSums(all.cores)>50]==cores[i]), type="p", col=i,cex=1.5, pch=16)
}
legend("bottomright", legend=cores, lty=1, pch=20, col=1:length(cores))
