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
offGridWeights1D<-function(s, gridList, NNSize=2,
                         mKrigObject=NULL, 
                         Covariance=NULL, covArgs=NULL,
                         aRange=NULL, sigma2=NULL, 
                         giveWarnings=TRUE,
                         debug=FALSE,
                         verbose=FALSE
                   ){
  #
  # function assumes the grid is 
  # integer locations and 1:m by 1:n
  # grid and off grid locations need to be transformed to that scale
  # 
  # also assumes the grid extends two cells beyond any off
  # e.g. s  coordinates should be between 
  # 2 and m-3 and 2 and n-3
  #
  # If mKrigObject (result of fitting model) is given 
  # extract all the covariance information from it. 
  # For the Matern family besides aRange and sigma2 is the 
  # smoothness
  np<- NNSize
  if( !is.null( mKrigObject)){
    sigma2<- mKrigObject$summary["sigma2"]
    aRange<- mKrigObject$summary["aRange"]
    Covariance<- mKrigObject$args$Covariance
    if( is.null(Covariance)){
      Covariance<- "Exponential"
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
  
  m<- length( gridList$x)
  
  dx<- gridList$x[2]- gridList$x[1]
  
  M<- nrow( s)
  # lower left corner of grid box containing the points
  s0<-  cbind( 
               trunc( (s[,1]- gridList$x[1] )/dx) + 1
               ) 
  
  # index  of locations when 2D array is unrolled
  s0Index<- as.integer( s0[,1])
  # check for more than one obs in a grid box
    tableLoc<- table( s0Index)
    allSingle<- all( tableLoc ==1 ) 
  
  theShift<- (0:(2*np-1)) - (np-1)
  xshift<- theShift
  
  
  nnX<- cbind( xshift)
  nnXCoords<- cbind( xshift*dx)
  
  #
  #  sX 
  
  
  sX<- s0[,1] + matrix( rep( xshift,M),
                        nrow=M, ncol=(2*np), byrow=TRUE)
  
  
  if( any( (sX < 1)| (sX>m)) ) {
    stop( "sX outside range for grid")
  }
  
  # indices of all nearest neighbors for unrolled vector.
  # this is an M by (2*np)^2 matrix where indices go from 1 to m*n
  # these work for the unrolled 2D array 
  # 
  sIndex<-  sX 
  # differences between nn and the off grid locations
  # for both coordinates
  # convert from integer grid to actual units. 
  differenceX<- (sX-1)*dx + gridList$x[1] - s[,1]
  
  
  dAll<- abs(differenceX)
  # pairwise distance among nearest neighbors. 
  dNN<- rdist(nnXCoords, nnXCoords)
  # cross covariances
  Sigma21Star<- sigma2* do.call(Covariance,
                                c(list(d = dAll/aRange), 
                                         covArgs)) 
  # covariance among nearest neighbors 
  Sigma11 <-  sigma2* do.call(Covariance,
                              c(list(d = dNN/aRange), 
                                covArgs))
  Sigma11Inv <- solve( Sigma11)
  # each row of B are the weights used to predict off grid point
  B <- Sigma21Star%*%Sigma11Inv
  # create spind sparse matrix
  # note need to use unrolled indices to refer to grid points
  ind<- cbind( rep(1:M, each= (2*np) ), c( t( sIndex)))
  ra<-  c( t( B))
  da<- c( M, m )
  spindBigB<-  list(ind=ind, ra=ra, da=da )
  # now convert to the more efficient spam format
  BigB<- spind2spam( spindBigB)
  #
  # prediction variances  
  # use cholesky for more stable numerics
  cholSigma11Inv<- chol(Sigma11Inv)
  # create spind sparse matrix of sqrt variances
  # or covariances to simulate prediction error. 
  w <- Sigma21Star%*%t(cholSigma11Inv)
  predictionVariance <-  sigma2 - rowSums(w^2)
  # easiest case of just one obs in each grid box  
  #  sigma2 - diag(Sigma21Star%*%Sigma11Inv%*%t(Sigma21Star) )
  spindObjSE<- list(ind=cbind( 1:M, 1:M),
                      ra=sqrt(predictionVariance),
                      da= c( M,M)
                    )
  BigSE<- spind2spam( spindObjSE)
  if(allSingle){
    duplicateIndex<-NA
  }
  if( !allSingle){
    indDuplicates<- (tableLoc > 1)
    if( giveWarnings){
    cat("Found", sum(indDuplicates), 
        "grid box(es) containing more than 1 obs location",
        fill=TRUE)
    }
    
    duplicateIndex<-names( tableLoc) [indDuplicates]
    duplicateIndex<-  as.numeric(duplicateIndex)
# duplicateIndex is the unrolled indices for all grid boxes with 
# 2 or more observations
# following code is written assuming there are not many of these. 
    nBox<- length( duplicateIndex) 
    indDupSE<-NULL
    raDupSE<- NULL
    for( k in 1:nBox){
      theBox<- duplicateIndex[k]
      # the obs that are in this box
      indBox<- which(s0Index == theBox)
      nDup<- length( indBox)
      dDup<- rdist( s[indBox,], s[indBox,])
      sigmaMarginal<- sigma2* do.call(Covariance,
                                      c(list(d = dDup/aRange), 
                                        covArgs))
      A<- w[indBox,]
      localSE2<-  sigmaMarginal - A%*%t(A)
      localSE<- t(chol( localSE2 ))
      # localSE %*% rnorm(nDup) will generate correct corrected 
      # prediction errors for obs in this grid box ("theBox")
      indTmp<- cbind(rep( indBox, nDup), rep( indBox, each=nDup) )
      raTmp<- c(localSE)
      indDupSE<- rbind( indDupSE,indTmp)
      raDupSE<-      c(  raDupSE, raTmp)
    }
    #print( dim(indDupSE ))
    #print( length(raDupSE))
  BigSE[indDupSE]<- raDupSE
  }

 if( debug){ 
    return(
      list( B= BigB, SE= BigSE, 
            predictionVariance = predictionVariance,
            Sigma11Inv = Sigma11Inv,
            Sigma21Star= Sigma21Star,
            s0Index = s0Index,
            s0 = s0,
            gridX = t( (sX-1)*dx + gridList$x[1]),
            gridList = gridList,
            duplicateIndex= duplicateIndex
            )
          )
 }
  else{
    return(
      list( B = BigB, 
            SE = BigSE,
            predictionVariance = predictionVariance )
    )
  }
  }
