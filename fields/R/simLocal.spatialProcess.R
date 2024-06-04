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
##END HEADER
"simLocal.spatialProcess" <- function(mKrigObject,  
    predictionGridList = NULL,
    simulationGridList = NULL, 
        gridRefinement = 1, 
                    np = 2, 
                     M = 1,
                    nx = 80,
                    ny = 80, 
               verbose = FALSE,
                 delta = NULL, giveWarnings=TRUE,
                  fast = FALSE, 
                NNSize = 5, 
                          ...)
#
# NOTE throughout $x is first dimension of the grid in a gridList but also 
# $x in mKrig object is the _matrix_ of locations
# 
    {
  
   nObs<- nrow( mKrigObject$x)
   sDimension<- ncol(mKrigObject$x)
   if ( sDimension > 2) {
        stop("conditional simulation only implemented for 1 and 2 dimensions")
   }
   
   if( sDimension == 1 & fast ){
     stop("fast prediction not implemented in 1 D")
   }
   
    # create prediction set of points based on what is passed
    # and if the grid is not specified 
   if (is.null(predictionGridList)) {
     # these adjustments insure there are enough grid
     # points beyond the range of the locations.
     # Put xr[1] in the middle of the  npth grid box
     # and xr[2] in to the  nx - np
     predictionGridList<- makePredictionGridList(
       mKrigObject=mKrigObject, 
       nx=nx, 
       ny=ny, 
       np=np
     )
   }
   
   nx <- length(predictionGridList$x)
   ny <- ifelse(sDimension>=2, length(predictionGridList$y) ,1 )
   
# 
# check that predictionGrid is equally spaced
# this is needed because of the fast simulation algorithm
    checkPredictGrid( predictionGridList) 
#
#
   if (is.null(simulationGridList)) {
     simulationGridList<- makeSimulationGrid(
       predictionGridList,
       gridRefinement)
   }
#
# #
# # round off the grids so that they match to 8 digits
# # that way prediction grid is  precisely a subset of
# # simulation grid
         predictionGridList$x<- signif(predictionGridList$x, 8)
         simulationGridList$x<- signif(simulationGridList$x, 8)
         if( sDimension ==2){
         predictionGridList$y<- signif(predictionGridList$y, 8)
         simulationGridList$y<- signif(simulationGridList$y, 8)
         }
         indexSubset<-  list( x=match(predictionGridList$x,
                                      simulationGridList$x))
         # # shortcut to avoid if statement for predicted in for
         # # loop
         # indexSubset$y = rep(1, length( indexSubset$x) )
         if( sDimension ==2){
         indexSubset$y=match(predictionGridList$y,
                                      simulationGridList$y)
         }
         else{
           indexSubset$y=1
         }
         
#  
# core covariance parameters from spatial model   
    tau <-    mKrigObject$summary["tau"]
    sigma2 <- mKrigObject$summary["sigma2"]
    aRange<-  mKrigObject$summary["aRange"]
    Covariance <- mKrigObject$args$Covariance
# wipe out some extraneous components that are not used by the Covariance
# function.
    covArgs0 <- mKrigObject$args
    covArgs0$Covariance<- NULL 
    covArgs0$distMat <- NULL
    covArgs0$onlyUpper<- NULL
    covArgs0$aRange<- NULL
    
#
# set up various arrays for reuse during the simulation
    nObs <- nrow(mKrigObject$x)
#  
    
    timeCESetup<- system.time(
    # set up object for simulating on a grid using circulant embedding
    CEObject<- circulantEmbeddingSetup(simulationGridList,
                                   cov.function = mKrigObject$cov.function,
                                       cov.args = mKrigObject$args,
                                          delta = delta )
    )[3]
    
    #
    if (verbose) {
        cat("dim of full circulant matrix ", CEObject$M, 
            fill = TRUE)
    }
#
# weights crucial to fast off grid simulation
#
      timeOffGridSetup <- system.time(
        offGridObject <- offGridWeights(
          mKrigObject$x,
          simulationGridList,
          mKrigObject,
          np = np,
          giveWarnings = giveWarnings
        )
      )[3]
   
    #
    # find conditional mean field from initial fit
      hHat <- predictSurface(mKrigObject,
                           gridList = predictionGridList,
                           fast=fast, 
                           NNSize= NNSize,
                            ...)$z
      
      sdNugget<- tau* sqrt(1/mKrigObject$weights)
#      
# setup output array to hold ensemble
#  in 1D case ny=1 
#    
    out <- array(NA, c( nx, ny, M))
    t1<-t2<- t3<- rep( NA, M)
    
##########################################################################################
### begin the big loop
##########################################################################################
    
     for (k in 1:M) {
       if( k%%10 ==0 ){
        cat(k, " ")
       }
        # simulate full field
        t1[k]<- system.time(
        hTrue<- as.matrix(sqrt(sigma2) * circulantEmbedding(CEObject))
        )[3]
        
        # NOTE: fixed part of model (null space) does not need to be simulated
        # because the  estimator is unbiased for this part.
        # the variability is still captured because the fixed part
        # is still estimated as part of the predict step below
        #
        t2[k]<- system.time(
        hData <- offGridObject$B%*%c(hTrue) + 
                    (offGridObject$SE)%*%rnorm(nObs) 
              )[3]
        ySynthetic <- hData + sdNugget*rnorm(nObs)
        #
        # predict at grid using these data
        # and subtract from synthetic 'true' value
        #
        t3[k]<-system.time(
            spatialError <- predictSurface.mKrig(mKrigObject,
                                 gridList = predictionGridList,
                                 ynew = ySynthetic,
                                 fast=fast, 
                                 NNSize= NNSize,
                                 giveWarnings = FALSE,
                                 ...)$z
         
             )[3]
 
        
        # add the error to the actual estimate  (conditional mean)
        # subset  hTrue to the prediction grid
        # note for 1D $y is 1. 
       
        out[,, k] <- hHat + (spatialError -  
                               hTrue[indexSubset$x,indexSubset$y])
       
     }
    
    cat(" ", fill=TRUE)
    
    return(list(x = predictionGridList$x,
                y = predictionGridList$y,
                z = out, 
                hHat= hHat,
                timing=c( CESetup=timeCESetup,
                          OffSetup=timeOffGridSetup,
                                  CE = median(t1), 
                             OffGrid = median(t2),
                                mKrig = median(t3)
                          ),
                gridRefinement=gridRefinement,
                M= CEObject$M,
                #simulationGridList= simulationGridList,
                timingFull = cbind( t1, t2,t3),
             call = match.call())
           )
}

