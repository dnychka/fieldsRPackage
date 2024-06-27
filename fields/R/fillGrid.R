fillGrid<- function( gridList, M){
  deltaX<-gridList$x[2]- gridList$x[1]
  m1<- length( gridList$x)
  xSeq<- seq( min(gridList$x),max(gridList$x),
              deltaX/M)
  deltaY<-gridList$y[2]- gridList$y[1]
  m2<- length( gridList$y)
  ySeq<- seq( min(gridList$y), max(gridList$y),
              deltaY/M)
  return( 
    list( x= xSeq, y=ySeq)
  )
}