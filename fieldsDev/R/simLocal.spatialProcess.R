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
    gridList = NULL,
                     M = 1,
                    nx = 80,
                    ny = 80, 
                 delta = NULL, 
                  fast = FALSE, 
                NNSize = 4, 
         NNSizePredict = 4, 
              #truncate = NULL, 
         giveWarnings = FALSE, 
               verbose = FALSE,
                          ...)
#
# NOTE throughout $x is first dimension of the grid in a gridList object  but also 
# $x in mKrig object is the _matrix_ of locations
# 
    {
  
   nObs<- nrow( mKrigObject$x)
   sDimension<- ncol(mKrigObject$x)
   
   if ( sDimension !=2) {
        stop("conditional simulation using fast, local algorithm 
             only implemented for  2 dimensions")
   }
   # great circle distance is not a stationary model for a regular lon/lat 
   # grid
   if( !is.null(mKrigObject$args$Distance)){
     if( (mKrigObject$args$Distance!="rdist") & fast){
       cat( mKrigObject$args$Distance, fill=TRUE)
       stop("can not use fast predict option with nonEuclidean distance.")
     }
   }
   # 
   # this is needed because of the fast simulation algorithm 
     gridInfo<- augmentPredictionGrid(
     gridList = gridList,
            s = mKrigObject$x, 
           nx = nx, 
           ny = ny,
       NNSize = NNSize,
     )
  # this will fill in gridList from nx, ny if not passed. 
     gridList<- gridInfo$gridList
     predictionGrid<- gridInfo$predictionGrid
     indexSubset<- gridInfo$indexSubset
     
#
# core covariance parameters from spatial model   
#
    tau <-    mKrigObject$summary["tau"]
    sigma2 <- mKrigObject$summary["sigma2"]
    aRange<-  mKrigObject$summary["aRange"]
    Covariance <- mKrigObject$args$Covariance
# just use mKrigObject for covariance info
    
   
#
# set up various arrays for reuse during the simulation
#  
####### NOTE there are two and possibly (if fast ==TRUE ) three SETUP steps  to 
# make conditional  simulation loops efficients by reusing some objects. 
#  
# SETUP weights for  fast, local off grid simulation
# grid is expanded to include nearest neighbors.
#
    timeOffGridSetup <- system.time(
      localObject <-  offGridWeights2D( s=mKrigObject$x,
                                        gridList=predictionGrid,
                                        NNSize=NNSize,
                                        mKrigObject=mKrigObject,
                                        Covariance=Covariance,
                                        covArgs=NULL,
                                        aRange=aRange,
                                        sigma2=sigma2,
                                        giveWarnings=giveWarnings,
                                        debug=FALSE,
                                        verbose=verbose)
        )[3]
#
# SETUP object when  fast==TRUE for rapid prediction
#
    timefastSetup<- NA
    fastObject<-NULL
    if(fast){
      if( verbose){
        cat("covariance function info", fill=TRUE)
        print(mKrigObject$cov.function.name)
        print(mKrigObject$args)
      }
      timefastSetup<- system.time(
        fastObject<- mKrigFastPredictSetup(
                                mKrigObject,
                          gridList = gridList,
                            NNSize = NNSizePredict,
                      giveWarnings = giveWarnings,
                           verbose = verbose) )[3]
    }
#
# SETUP object for simulating on a grid using circulant embedding (CE)
#
    timeCESetup<- system.time(
    CEObject<- circulantEmbeddingSetup(predictionGrid,
                                   cov.function = mKrigObject$cov.function,
                                       cov.args = mKrigObject$args,
                                       #truncate=truncate,
                                   delta = delta
                                       )
    )[3]
    if( verbose){
      cat("info on circulant embedding", fill=TRUE)
      cat( "ratioRMSE", CEObject$ratioRMSE, fill=TRUE)
      cat("dim of full circulant matrix ", CEObject$M, 
            fill = TRUE)
    }
#
# find conditional mean field from initial fit
#.    
      if( verbose){
        print(mKrigObject$args)
      }
      hHat <- predictSurface.mKrig(mKrigObject,
                           gridList = gridList,
                               fast = fast, 
                        setupObject = fastObject,
                        #setupObject = NULL,
                        #     NNSize = NNSizePredict,
                                  ...)$z
      sdNugget<- tau* sqrt(1/mKrigObject$weights)
#      
# setup output array to hold ensemble
#    
    out <- array(NA, c( length(indexSubset$x),
                        length(indexSubset$y), 
                        M))
    t1<-t2<- t3<- rep( NA, M)
    #outData<- matrix( NA,nObs,M )
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
        # NOTE: linear, fixed part of model (null space) does not need to be simulated
        # because the estimator is unbiased for this part.
        # the variability is still captured because the fixed part
        # is still estimated as part of the predict step below
        #
        # unroll hTrue from grid and do sparse matrix multiply
        # 
        localError<- rnorm(nObs)
        # sparse matrix multiplication of weights by process on grid
        # and add in independent error 
        t2[k]<- system.time(
        hData <- localObject$B%*%c(hTrue) + 
                    (localObject$SE)%*%localError
              )[3]
        ySynthetic <- hData + sdNugget*rnorm(nObs)
        #outData[,k] <-  ySynthetic
        #
        # predict at grid using these data
        # and subtract from synthetic 'true' value
        #
        t3[k]<-system.time(
            spatialPrediction <- predictSurface.mKrig(mKrigObject,
                                 gridList = gridList,
                                 ynew = ySynthetic,
                                 fast=fast, 
                                 NNSize= NNSizePredict,
                                 setupObject= fastObject,
                                 giveWarnings = TRUE,
                                 ...)$z
         
             )[3]
       
        # # add the error to the actual estimate  (conditional mean)
        # # subset  hTrue to the prediction grid
        # # note for 1D $y is 1.
        # NOTE subset hTrue because it may be generated on a larger grid to include nearest neighbors
        # for local simulation. 
        tmpError<- (spatialPrediction  - hTrue[indexSubset$x,indexSubset$y])
        
        # hHat is of course the conditional mean ...
        out[,, k] <- (hHat + tmpError)
     }
    
    cat(" ", fill=TRUE)
    timingTable<- c(  CESetup = timeCESetup,
                      SetupSim = timeOffGridSetup,
                      SetupPredict = timefastSetup,
                      medCE = median(t1), 
                      medOffGridSim = median(t2),
                      medPredict = median(t3)
    )
    return(list(x = predictionGrid$x[indexSubset$x],
                y = predictionGrid$y[indexSubset$y],
                z = out, 
                hHat= hHat,
                timing=timingTable,
                M= CEObject$M,
                timingFull = cbind( t1, t2,t3),
                #outData=outData,
                #outhData = outhData,
             call = match.call(),
             CEObject = CEObject
             )
           )
}

