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
offGridWeights2D<-function(s, gridList, NNSize=2,
                         mKrigObject=NULL, 
                         Covariance=NULL, covArgs=NULL,
                         aRange=NULL, sigma2=NULL, 
                         giveWarnings=TRUE,
                         debug=FALSE,
                         findCov=TRUE,
                         verbose=TRUE
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
  # setup the grid info from which to interpolate
  np<- NNSize
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
  # lower left corner of grid box containing the points
  s0<-  cbind( 
               trunc( (s[,1]- gridList$x[1] )/dx) + 1 ,
               trunc( (s[,2]- gridList$y[1] )/dy) + 1
               ) 
  
  # index  of locations when 2D array is unrolled
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
  #  the unrolled row and column indices for  np=2
  #  yield 16 nearest neighbors 
  
  sX<- s0[,1] + matrix( rep( xshift,M),
                        nrow=M, ncol=(2*np)^2, byrow=TRUE)
  
  sY<- s0[,2] + matrix( rep( yshift,M),
                        nrow=M, ncol=(2*np)^2, byrow=TRUE)
  
  if( (min(sX) < 1)| (max(sX)>m) | (min(sY) < 1)| (max(sY)>n)  ) {
    cat(" grid must extend at least", NNSize + 1,
        "grid points beyond the observation locations", fill=TRUE)
  }
  # indices of all nearest neighbors for unrolled vector.
  # this is an M by (2*np)^2 matrix where indices go from 1 to m*n
  # these work for the unrolled 2D array 
  # 
  sIndex<-  sX + (sY-1)*m
  # differences between nn and the off grid locations
  # for both coordinates
  # convert from integer grid to actual units. 
  differenceX<- (sX-1)*dx + gridList$x[1] - s[,1]
  differenceY<- (sY-1)*dy + gridList$y[1] - s[,2]
  # print(  cbind( t( (sX-1)*dx + gridList$x[1])[,1],  
  #         t( (sY-1)*dy + gridList$y[1])[,1]))
  # all pairwise distances between each off grid and 
  # (2*np)^2  ( np=2 has 16) nearest neighbors 
  dAll<- sqrt(differenceX^2 + differenceY^2)
  # pairwise distance among nearest neighbors. 
  dNN<- rdist(nnXYCoords, nnXYCoords )
  # cross covariances
  #print( sigma2)
  Sigma21Star<- sigma2* do.call(Covariance,
                                c(list(d = dAll/aRange), 
                                         covArgs)) 
  # covariance among nearest neighbors 
  Sigma11 <-  sigma2* do.call(Covariance,
                              c(list(d = dNN/aRange), 
                                covArgs))
  # OLD code
  #  Sigma11Inv <- solve( Sigma11)
  # replaced with reduced rank interpolation. 
  
  eigenSigma<- eigen( Sigma11, symmetric=TRUE)
  U<- eigenSigma$vectors
  dV<- eigenSigma$values
  
 
  dVInv<- ifelse(dV/max(dV) >= 1e-10, 1/dV,0)  # note this will exclude negative values
  
  if( verbose){
    cat("conditional number of Sigma11", fill=TRUE)
    cat(max(dV)/ min(dV), fill=TRUE)
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
  
  # print( w[1:5,])
  # print( rowSums(w^2)[1:5])
  # 
  predictionVariance <-  sigma2 - rowSums(w^2)
  # set to NA negative values -- due to colinearity/roundoff
  predictionVariance<- ifelse( predictionVariance <= 0, NA,predictionVariance )
  if( findCov|debug){
  # don't find covariances if findCov is FALSE
  # but debug=TRUE overrides that 
  # this saves computation for the cases with multiple 
  # obs in a single grid box.
  
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
  }
  else{
    BigSE=NA
    
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
            gridY = t((sY-1)*dy + gridList$y[1]),
            gridList = gridList,
            duplicateIndex= duplicateIndex
            )
          )
 }
  else{
    return(
      list( B = BigB, 
            SE = BigSE,
            predictionVariance = predictionVariance,
            gridList= gridList
            )
    )
  }
  }
