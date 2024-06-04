plotMatrix<- function( A, M=5,N=5,...){
  
  n<- nrow(A)
  m<- ncol( A)
  
  obj<-  ((A)[n:1,] )
  
  imagePlot( obj,  axes=F, xlab="columns", ylab="rows",
             ...)
   nG<- round(pretty( 1:n, n=N) )
  axis(side=1,  at= (nG-1)/(n-1), labels= format( nG) )
  mG<-  seq( 1,m, round(m/M) )
  mG<- round( pretty( 1:m, n=M))
  axis(side=2,  at= 1- (mG-1)/(m-1),
       labels= format( mG), gap.axis=0, las=2 )
  box()
}