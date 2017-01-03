#sensitivity values

#1 WA optima
#2 WA optima - two Wsd
#3 min 5 occurences - ISI
#4 5% of distribution - BQI
#5 5% of maximum


library(rioja)
data(SWAP)
x<-SWAP$pH
y<-SWAP$spec
y<-y[,colSums(y)>20]

sense<-list()

wamod<-WA(y,x, tolDW=T)

sense$waopt<-wamod$coef[,1]
sense$waopt2tol<-wamod$coef[,1]-2*wamod$coef[,2]

sense$ISI<-apply(y,2,function(sp)sort(x[sp>0])[5])

sense$BQI<-apply(y,2,function(sp){
  mod<-glm(sp/100~x+I(x^2),family=quasibinomial, weights=rep(100,length(x)))
  nd<-data.frame(x=seq(min(x),max(x),length=100))
  pred<-predict(mod,newdata=nd, type="response")
  spred<-sum(pred)
  cs<-cumsum(pred)
  if(cs[1]>0.05*spred)comm<-min(x)
  else{
    x2<-nd$x[cs<0.05*spred]
    comm<-x2[length(x2)]  
  }
  comm
})




sense$comm<-apply(y,2,function(sp){
  mod<-glm(sp/100~x+I(x^2),family=quasibinomial, weights=rep(100,length(x)))
  nd<-data.frame(x=seq(min(x),max(x),length=100))
  plot(sp/100~x)
  points(x,fitted(mod), col=2, pch=16)
  pred<-predict(mod,newdata=nd, type="response")
 # browser()
  if(pred[1]>0.05*max(pred))comm<-min(x)
  else{
    p2<-pred<0.05*max(pred)
    p3<-cumsum(!p2)
    x2<-nd$x[p2&p3==0]
    comm<-x2[length(x2)]
    
  }
  
  comm
})


pairs(sense)
                           