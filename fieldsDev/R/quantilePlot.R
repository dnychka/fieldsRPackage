# testing 
# 
# y<- rnorm( 60);  quantilePlot(y)
# y<- rgamma( 60, shape=1);
# quantilePlot(y, quantileFunction="qgamma", quantileArgs=list(shape=1))
# abline(0,1)
quantilePlot<- function(y,
                        quantileFunction="qnorm", 
                         quantileArgs= NULL,alpha=.05, M=100,
                        smoothQuantiles=FALSE,dfSmooth=NA,
                        col = "thistle1", lineCol = "thistle3",
                        ...){
  n<- length(y)
  sortedY<- sort(y)
  
  P<- ((1:n)-.5)/n
  Q<- do.call(quantileFunction, c(list( p=P ), quantileArgs))
  # use the fact  quantile of uniform RV gives a draw from the distribution of
  # the quantile. 
  bigZ<- runif(n*M)
  bigZ<- do.call(quantileFunction, c(list( p=bigZ ), quantileArgs) )
  
  bigZ<- matrix( bigZ, n,M)
  bigZ<- apply( bigZ, 2, sort)
  
  # rely on CLT even for extreme quantiles!
  sdZ<- apply( bigZ, 1, sd)
  meanZ<- apply( bigZ, 1, mean)
  
  if( smoothQuantiles){
    sdZ<- c(sreg( Q, sdZ, df=dfSmooth)$fitted.values)
    meanZ<- c(sreg( Q,meanZ, df=dfSmooth)$fitted.values)
  }
  
  # simultaneous band this is exact up to Monte Carlo error
  U<- abs( (bigZ-meanZ)/sdZ)
  UMax<- apply( U, 2, max)
  # sup bound that contains 1-alpha of the simulated samples. 
  BSimultaneous<- quantile( UMax, (1-alpha))
# normal  family here assuming quantile distribution is normal 
  upper<- meanZ+ qnorm(1-alpha/2)*sdZ
  lower<-  meanZ -  qnorm(1-alpha/2)*sdZ
  
  upperS<- meanZ+ BSimultaneous*sdZ
  lowerS<-  meanZ -  BSimultaneous*sdZ
  plot(Q, sortedY,type="n", xlab="Theoretical quantiles",
       ylab="Data quantiles",... )
  envelopePlot(Q, lowerS,Q, upperS,col =col, lineCol = lineCol)
  
  matlines( Q, cbind( lower, upper), lty=2, lwd=1.5, col=lineCol)
  
  #matlines( Q,bigZ,lty=1, col="red1" )
  points( Q,sortedY, pch=16, col="grey30")
  
  
}