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
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2
##END HEADER
spatialProcess <- function(x, y,  weights = rep(1, nrow(x)),   Z = NULL,
                     ZCommon = NULL,
                  mKrig.args = NULL,
                cov.function = NULL, 
                  	cov.args = NULL,
                      parGrid = NULL, 
                      reltol = 1e-4,
                       na.rm = TRUE,
                  	 verbose = FALSE,
                        REML = FALSE, 
            cov.params.start = NULL,
                       gridN = 5,
               profileLambda = FALSE,
               profileARange = FALSE,
               profileGridN  = 15, 
                  gridARange = NULL,
                  gridLambda = NULL,
                     CILevel = .95,
                       iseed = 303, 
         collapseFixedEffect = TRUE,
               simpleKriging = FALSE,
                           ...) {
 
# Through out this code:  obj is the output list and will be added to throughout the computation
# This device is typical in fields functions where a list is built up in the process of the 
# computation. 
  
#  THE RULES: 
# Being in the cov.params.start list means a parameter will be optimized
# by maximum likelhood. 
# If a  parameter is in the cov.args list it will be fixed at that value
#  
  
# NOTE all the  ... extra arguents are assumed to be for the cov.args list
  GCV<- FALSE # A top level placeholder to add GCV search capability once
# algorithm is stable 
# this default choice is used in the next level functions.
############################################################ 
#### Setting some defaults for the covariance function and parameters
############################################################ 
  

  
# It is also convenient to just add a few arugments of the covariance driectly in the 
# this list collects those.  E.g. spatialProcess( x, y, Covariance="Exponential" )
# will pass this choice along to the default covariance function, stationary.cov
#
  
  extraArgs<- list(...)
  
   if( REML&GCV){
     stop("Cannot optimize for both REML and GCV!")
   }
  
# set defaults based on the model passed. 
# the following function fills in some typical models and 
# just avoids some typing and makes for some clean examples
#   
# note this also creates the initial components for the output list. 
#
   obj<- spatialProcessSetDefaults(x, 
                            cov.function = cov.function, 
                                cov.args = cov.args,
                        cov.params.start = cov.params.start,
                              mKrig.args = mKrig.args,
                               extraArgs = extraArgs,
                                 parGrid = parGrid,
                                   gridN = gridN,
                     collapseFixedEffect =  collapseFixedEffect,
                           simpleKriging = simpleKriging,
                                 verbose = verbose)
   #
   # obj$CASE 
   # 0 evaluate on passed cov parameters but MLEs for sigma, tau found from
   #   lambda
   #
   # 1 optimize loglikelihood over any parameters specified in 
   #   cov.params.start but not in cov.args
   #
   # 2  grid search over parameters using parGrid and generating starting values for 
   #   for the MLEs 
   #
   # 3 profile over lambda and/or aRange 
   # this is more computationally demanding. 
   
  if( verbose){
    cat(" The CASE:", obj$CASE, fill=TRUE)
    cat("Complete list of components in cov.args: ", "\n",
        names(obj$cov.args),fill=TRUE )
    cat("Complete list of components in mKrig.args: ", "\n",
        names(obj$mKrig.args),fill=TRUE )
  }
   
  if( verbose){
     cat( "Names cov.args:","\n", names( obj$cov.args), fill=TRUE)
     cat( "Names cov.params.start:","\n", names(cov.params.start), fill=TRUE)
     cat( "Names argsFull:","\n", names( obj$cov.argsFull),fill=TRUE  )
   }
   
   
####################################################################  
# CASE 2 grid search for starting values 
####################################################################  
   
  if( (obj$CASE == 2 )  ){
    if( verbose){
      
      cat("*****************************************", fill=TRUE)
      cat("***** Grid search ", fill=TRUE)
      cat("***** Grid search starting values: ", fill=TRUE)
      print(cov.params.start )
      cat("***** parGrid names: " , fill=TRUE)
      print( names( parGrid))
      cat("***** the grid to search:", fill=TRUE)
      print( names( obj$parGrid))
    }
    
    InitialGridSearch<- mKrigMLEGrid(x, y,  
                             weights = weights,
                                   Z = Z, 
                             ZCommon = ZCommon,
                          mKrig.args = mKrig.args,
                        cov.function = obj$cov.function, 
                           cov.args  = obj$cov.args,
                            par.grid = obj$parGrid, 
                              reltol = reltol,
                               na.rm = na.rm,
                             verbose = verbose,
                                REML = REML,
                                 GCV = GCV,
                    cov.params.start = cov.params.start)
  # use grid search to set starting values
    if( all(is.na(InitialGridSearch$summary) )) {
      cat("spatialProcess: Problems with optim in grid search", fill=TRUE)
      cat("returned object includes the likelihood evaluations up to the 
         error", fill=TRUE)
      InitialGridSearch$optimSuccess<- FALSE
      InitialGridSearch$call<- match.call()
      return(InitialGridSearch)
    }
    
    if( verbose){
      cat("results of intitial grid search: ", fill=TRUE)
      print(InitialGridSearch$indMax )
    }
    parNames<- names( obj$parGrid)
    if( is.null( cov.params.start)){
      cov.params.start<- obj$parGrid[InitialGridSearch$indMax,]
      names(cov.params.start )<- parNames
    }
    else{
    cov.params.start[parNames] <- 
       obj$parGrid[InitialGridSearch$indMax, parNames]
    }
 }
   
####################################################################  
# CASES 1 , 2 , 3, 4
#################################################################### 
   if(obj$CASE !=0 ){
     # optimze over all parameters 
     # where starting values are given or if
     # values in cov.args are omitted. 
     obj$cov.params.start<- cov.params.start
     MLEInfo <-mKrigMLEJoint(x, y,  weights = weights, Z = Z, 
                             ZCommon = ZCommon,
                             mKrig.args = obj$mKrig.args,
                             cov.function = obj$cov.function, 
                             cov.args  = obj$cov.args,
                             na.rm = na.rm,
                             reltol=reltol,
                             cov.params.start = cov.params.start,
                             REML = REML,
                             GCV = GCV,
                             hessian = TRUE,
                             verbose = verbose,
                             iseed = iseed)
       
     
   
   if( is.na(MLEInfo$summary[1])){
     cat("spatialProcess: Problems with optim in mKrigMLEJoint ", 
         fill=TRUE)
     cat("returned object includes the likelihood evaluations up to the 
         error", fill=TRUE)
     MLEInfo$optimSuccess<- FALSE
     MLEInfo$call<- match.call()
     return(MLEInfo)
   }
    
  covarianceMLE<-  -1*solve(MLEInfo$optimResults$hessian)
  HessianResults<- diag( covarianceMLE)
   if( any( HessianResults < 0) ) {
       warning("Numerical hessian from optim indicates
                   MLE is not a maximum")
     }
      
   }
   
   
################################################################################
# final fit 
# now fit spatial model with MLE(s) 
# or the value(s) supplied in the call
# reestimate the other parameters for simplicity to get the
# complete mKrig object
####################################################################
   
# if all parameters are fixed -- don't mess with cov.args   
   if( obj$CASE == 0){
     obj$cov.argsFull <-  obj$cov.args
   }
   else{
      dupParameters<- match( names(MLEInfo$pars.MLE ), names(cov.args) )
     if( all( is.na(dupParameters)) ){
       obj$cov.argsFull<- c( obj$cov.args,
                             as.list(MLEInfo$pars.MLE) )
     }
   }
   
  
  
  mKrigObj <- do.call( "mKrig", 
	                c( list(x=x,
	                        y=y,
	                  weights=weights,
	                        Z=Z,
	                  ZCommon=ZCommon),
	                  obj$mKrig.args,
	             list( na.rm=na.rm),
	             list(cov.function = obj$cov.function),
	                  obj$cov.argsFull 
	            	)
             	)

####################################################################
# sort out output object based on the different cases
# also copy some information from the call and within function
# the output list, obj
####################################################################
  obj$GCV    <- GCV
  obj$REML   <- REML
  obj$CILevel<- CILevel
  
  if( obj$CASE==0){
    obj$MLEInfo <- NULL
    obj$MLESummary <- mKrigObj$summary
    obj$InitialGridSearch<- NULL
    obj$parameterCovariance<- NULL
  }
  
  if( obj$CASE==1){
    InitialGridSearch<- NULL
  }
  
####################################################################
#  Fill in all info related to finding MLE (Not CASE 0)
####################################################################
  if( obj$CASE!=0){
    # logical to distinguish from optim failure
    obj$optimSuccess<-TRUE
    #
    obj$InitialGridSearch<- InitialGridSearch
    obj$MLEInfo<- MLEInfo
    obj$MLESummary <- MLEInfo$summary
    # NOTE: covariance for lambda and ARange based on log(lambda) and log(ARange)
    obj$parameterCovariance<- solve(
            -1*MLEInfo$optimResults$hessian)
####################################################################
# Approximate large sample confidence intervals on the transformed scale
# followed by
#	then transforming back to original scale (see MLEInfo$par.transform)
# These are filled with NAs when numerical Hessian is not
# positive definite
####################################################################	
    obj$CITable<- confidenceIntervalMLE(obj, CILevel)
  }
  
  # combine everything into the output list, mKrig components first. 
  obj <- c( mKrigObj,obj)
  # replace call in mKrig  object with the top level one
  # from spatialProcess
  obj$call<- match.call()	
  
  class(obj) <- c( "spatialProcess","mKrig")
  
####################################################################
# Profiling depends on complete spatial process obj 
# which is why this is last
####################################################################	
  
  if (profileLambda) {
    obj$profileSummaryLambda <- profileMLE(obj, "lambda",
                                           parGrid = gridLambda,
                                           gridN = profileGridN
                                           )$summary
  }
  else{
    obj$profileSummaryLambda <- NULL
  }
  if (profileARange) {
    obj$profileSummaryARange <- profileMLE(obj, "aRange",
                                           parGrid = gridARange,
                                           gridN = profileGridN,
                                           )$summary
  }
  else{
    obj$profileSummaryARange <- NULL
  }
  

	return(obj)
}
