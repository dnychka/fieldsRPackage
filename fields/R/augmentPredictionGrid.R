augmentPredictionGrid <- function(gridList=NULL, s, nx=NULL, ny=NULL, NNSize,
                                  verbose=FALSE){
  #
  # NNSize used to add additional points so offGridWeights
  # has enough neighbors.
  #
  M<- ncol( s)
  if( M != 2){
    stop("augment grid only works for 2D")
  }
  
  
  # check if data is outside range of the grid  and if augmentation is  even needed !
  # 
  gridListOld<- gridList
  
  
  if( !is.null(gridList)){
    testforxyNames<- all( names(gridList )%in% c( "x","y"))
    if( !testforxyNames){
      print(names(gridList ) )
      stop("gridList has to have component names x and y ")
    }
    xr<- range( gridList$x)
    dx<- (gridList$x[2] - gridList$x[1] )
    yr<- range( gridList$y)
    dy<- (gridList$y[2] - gridList$y[1] )
# data outside grid ???
    if( (min(s[,1])< xr[1])|(max(s[,1])>xr[2])|
        (min(s[,2])< yr[1])|(max(s[,2])>yr[2])) {
      cat("range s[,1]", range(s[,1]),fill=TRUE )
      cat("range gridList$x", xr,fill=TRUE )
      cat("range s[,2]", range(s[,1]),fill=TRUE )
      cat("range gridList$y", yr,fill=TRUE )
      stop(" some observations outside grid")
    }
# gridList OK as is ?    
    indLx<- (min(s[,1]) - xr[1])/dx 
    indRx<- (xr[2]- max(s[,1]))/dx 
    indLy<- (min(s[,2]) -   yr[1])/dy 
    indRy<- (yr[2] - max(s[,2]))/dy
    allMargins<- c(indLx,indRx,indLy,indRy    )
    # take into account round off. 
    OKGrid<-  all( abs(allMargins -NNSize)<=1e-10)
    if(OKGrid)
    { 
      # grid is OK just use it!
      return(
        list(predictionGrid= gridList,
             
             indexSubset= list( x = 1:length(gridList$x),
                                y = 1:length(gridList$y)),
             expandGrid=FALSE)
      )
      }
  }
  
# creating a new grid or expanding the one passed 
  
  if(is.null(gridList) ){
    if(is.null(nx)|is.null(ny)){
      stop("need to specify nx and ny")
    }
    sDimension <- ncol(s)
    xr <- range(s[, 1])
    dx <- (xr[2] - xr[1]) / (nx-1) 
    yr <- range(s[, 2])
    dy <- (yr[2] - yr[1]) / (ny-1)
    gridListOld<- list( x= seq( xr[1], xr[2],dx),
                        y= seq( yr[1], yr[2],dy)
                       )
  }
  else{
    nx<- length(gridList$x)
    xr<- range( gridList$x)
    dx<- (gridList$x[2] - gridList$x[1] )
    ny<- length(gridList$y)
    yr<- range( gridList$y)
    dy<- (gridList$y[2] - gridList$y[1] )
  }     
  #
  # now create expanded grid. 
  # right side has an extra row/column of grid boxes to avoid 
  # how largest locations are referenced to a grid box 
  # there may be other routes to deal with this but this seems the most straight 
  # forward
  # in approximateCovariance2D here is how the lower left index is found:
  #
  # cbind( 
  #   trunc( (sObs[,1]- predictionGrid$x[1] )/dx) + 1 ,
  #   trunc( (sObs[,2]- predictionGrid$y[1] )/dy) + 1
  # ) 
  # 
  xg <- seq( xr[1] - (NNSize-1)*dx, xr[2] + (NNSize)*dx, dx)
  yg <- seq( yr[1] - (NNSize-1)*dy, yr[2] + (NNSize)*dy, dy)
  predictionGrid <- list(x=xg, y = yg)
  indexSubset<- list( x= (1:nx)+(NNSize-1), y= (1:ny)+(NNSize-1))
return(
       list(predictionGrid= predictionGrid,
            indexSubset=indexSubset,
            expandGrid=TRUE,
            gridList=gridListOld)
      )
}

