`MaRaTe.internal`<-function(X,Y=NULL,test)
{
 n<-dim(X)[1]
 p<-dim(X)[2]
 Xsigns<-sign(X)
 Xranks<-apply(abs(X),2,rank)
 if(is.null(Y))  # one sample
 {
  scores<-switch(test,"sign"= Xsigns, "normal"=Xsigns*qnorm(0.5+Xranks/(2*(n+1))), "rank"=Xsigns*Xranks)
  A<-t(scores)%*%scores/switch(test,"rank"=(n+1)^2,1)
  S<-colSums(scores)/switch(test,"rank"=(n+1),1)
  test.statistic<-as.numeric(t(S)%*%solve(A)%*%S)
  p.value<-1-pchisq(test.statistic,p)
  return(list(test.statistic=test.statistic,p.value=p.value))
 }
 # else two samples
 m<-dim(Y)[1]
 
 XYdata<-rbind(X,Y)
 XYranks<-apply(XYdata,2,rank)
 E<-switch(test, "sign"={medians<-apply(XYdata,2,median); 
                          as.matrix(ifelse(XYdata<=medians,1,0));},
                  "rank"=XYranks/(n+m+1),
                  "normal"=qnorm(XYranks/(n+m+1)))
 E.X<-E[1:n,]
 E.Y<-E[(n+1):(n+m),]
 T.X<-colMeans(as.matrix(E.X))
 T.Y<-colMeans(as.matrix(E.Y))
 E.bar<-colMeans(E)
 W.x <- t(E.X)%*%E.X - n*(E.bar%*%t(E.bar))
 W.y <- t(E.Y)%*%E.Y - m*(E.bar%*%t(E.bar))
 W<-(1/(n+m))*(W.x+W.y)
 W.inv<-solve(W)
 L<-as.numeric(n*t(T.X-E.bar)%*%W.inv%*%(T.X-E.bar)+m*t(T.Y-E.bar)%*%W.inv%*%(T.Y-E.bar))
                        
 p.value<-1-pchisq(L,p)

 list(test.statistic=L,p.value=p.value,Cov=W,df=p)
}
