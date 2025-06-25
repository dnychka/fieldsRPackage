#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2024 Colorado School of Mines
# 1500 Illinois St., Golden, CO 80401
# Contact: Douglas Nychka,  douglasnychka@gmail.com,
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2
##END HEADER\
approximateCovariance2D<-function(s, gridList, NNSize=4,
                         mKrigObject=NULL, 
                         Covariance=NULL, covArgs=NULL,
                         aRange=NULL, sigma2=NULL, 
                         giveWarnings=TRUE,
                         debug=FALSE,
                         verbose=FALSE
                   ){
  np<- NNSize
  #
  # function works with  the grid as 
  # integer locations and 1:m by 1:n
  # thus the grid and off grid locations need to be transformed to that scale
  # 
  # also assumes the grid extends NNSize cells beyond any off grid (data) locations
  # e.g.  the s  coordinates should be between 
  # NNSize and m- NNSize -1  and NNSize and n- NNSize -1
  #
  # setup the grid info from which to interpolate
  m<- length( gridList$x)
  n<- length( gridList$y)
  
  dx<- gridList$x[2]- gridList$x[1]
  dy<- gridList$y[2]- gridList$y[1]
  
  # If mKrigObject (result of fitting model) is given 
  # extract all the covariance information from it. 
  # For the Matern family besides aRange and sigma2 is the 
  # smoothness
  if( !is.null( mKrigObject)){
    sigma2<- mKrigObject$summary["sigma2"]
    aRange<- mKrigObject$summary["aRange"]
    Covariance<- mKrigObject$args$Covariance
    if( is.null(Covariance)){
      Covariance<- "Matern"
      
    }
    
    covArgs<-mKrigObject$args 
    
  # some R arcania -- strip out all arguments used by say stationary.cov
  # but not used by the Covariance function 
  # Do not want to call the covariance function with these extra args. 
    if( !is.null( covArgs) ){
      argNames<- names( as.list( get(Covariance)))
      argNames<- argNames[ -length(argNames)]
      ind<- match(  names(covArgs), argNames)
      covArgs[is.na( ind)] <- NULL
    }
  }
  #
  # wipe out aRange to prevent it from being used twice
  # in do.call below assumed to be 1.0 with scaled distances
  #
  covArgs$aRange<- NULL
  
 
  M<- nrow( s)
  # lower left corner of grid box containing the observation locations
  s0<-  cbind( 
               trunc( (s[,1]- gridList$x[1] )/dx) + 1 ,
               trunc( (s[,2]- gridList$y[1] )/dy) + 1
               ) 
  
  # index  of locations when 2D array is unrolled
  # for the grid point in the lower left of each observation
  s0Index<- as.integer( s0[,1] + (s0[,2]-1)*m)
  # check for more than one obs in a grid box
    tableLoc<- table( s0Index)
    allSingle<- all( tableLoc ==1 ) 
  # np=2
  # specific to 2nd degree neighborhood
  #   (2*np)^2 = 16  points total 
  #xshift<- rep( c(0,1,2,3), 4)
  #yshift<- rep( c(0,1,2,3), each=4 )
  
  theShift<- (0:(2*np-1)) - (np-1)
  xshift<- rep( theShift, 2*np)
  yshift<- rep( theShift, each=2*np )
  
  nnXY<- cbind( xshift, yshift)
  nnXYCoords<- cbind( xshift*dx, yshift*dy)
  
  #
  #  sX and sY are M by (2*np)^2 matrices  where each row is
  #  the unrolled row and column indices 
  #  for  NNSize=2
  #  yields 16 nearest neighbors and so they form 
  #  nine  3X3 grid boxes with the observation location in the center one. 
  
  sX<- s0[,1] + matrix( rep( xshift,M),
                        nrow=M, ncol=(2*np)^2, byrow=TRUE)
  
  sY<- s0[,2] + matrix( rep( yshift,M),
                        nrow=M, ncol=(2*np)^2, byrow=TRUE)
  
  
  n<- length( gridList$y)
  
  
  dy<- gridList$y[2]- gridList$y[1]
  
  if( any( (sX < 1)| (sX>m)) ) {
    cat( min( sX), max(sX), fill=TRUE)
    stop("sX outside range of 1:m, 
         increase the range of the x grid")
  }
  if( any( (sY < 1)| (sY>n)) ) {
     cat( min( sY), max(sY), fill=TRUE)
    stop("sY outside range for grid 1:n, 
         increase the range of the y grid")
  }
  # indices of all nearest neighbors for unrolled vector.
  # this is an M by (2*np)^2 matrix where indices go from 1 to m*n
  # these work for the unrolled 2D array 
  # 
  sIndex<-  sX + (sY-1)*m
  # differences between nearest neighbor grid points 
  # and the off grid locations
  # for both coordinates
  # convert from integer grid to actual units. 
  differenceX<- (sX-1)*dx + gridList$x[1] - s[,1]
  differenceY<- (sY-1)*dy + gridList$y[1] - s[,2]
  # print(  cbind( t( (sX-1)*dx + gridList$x[1])[,1],  
  #         t( (sY-1)*dy + gridList$y[1])[,1]))
  # all pairwise distances between each off grid and 
  # (2*np)^2  ( np=2 has 16) nearest neighbors 
  dAll<- sqrt(differenceX^2 + differenceY^2)
  # cross covariances
  Sigma21Star<- sigma2* do.call(Covariance,
                                c(list(d = dAll/aRange), 
                                         covArgs)) 
  # pairwise distance among nearest neighbors. 
  dNN<- rdist(nnXYCoords, nnXYCoords )
  # covariance among nearest neighbors 
  #
  Sigma11 <-  sigma2* do.call(Covariance,
                              c(list(d = dNN/aRange), 
                                covArgs))
  eigenSigma<- eigen( Sigma11, symmetric=TRUE)
  U<- eigenSigma$vectors
  dV<- eigenSigma$values
  if( verbose){
  cat("conditional number of Sigma11", fill=TRUE)
  cat(max(dV)/ min(dV), fill=TRUE)
  }
  dVInv<- ifelse(dV/max(dV) >= 1e-10, 1/dV,0)  # note this will exclude negative values
  
  if( verbose){
  if( any(dVInv==0 )){
    cat(sum( dVInv==0), " eigenvalues set to zero out of ", 
        length( dVInv), " Sigma11", fill=TRUE)
  }
  }
  
  Sigma11Inv <- U%*% diag( dVInv)%*%t(U)
  
  # each row of B are the weights used to predict off grid point
  B <- Sigma21Star%*%Sigma11Inv
  # create spind sparse matrix
  # note need to use unrolled indices to refer to grid points
  ind<- cbind( rep(1:M, each= (2*np)^2 ), c( t( sIndex)))
  ra<-  c( t( B))
  da<- c( M, m*n )
  BigB<-  list(ind=ind, ra=ra, da=da )
  # now convert to the more efficient spam format
  BigB<- spind2spam( BigB)
  
  
  
 if( debug){ 
    return(
      list( B= BigB, SE= NA, 
            Sigma11Inv = Sigma11Inv,
            Sigma21Star= Sigma21Star,
            s0Index = s0Index,
            s0 = s0,
            gridX = t( (sX-1)*dx + gridList$x[1]),
            gridY = t((sY-1)*dy + gridList$y[1]),
            gridList = gridList
            )
          )
 }
  else{
    return(
      list( gridList=gridList, B = BigB )
    )
  }
  }
