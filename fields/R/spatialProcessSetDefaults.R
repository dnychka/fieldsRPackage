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
spatialProcessSetDefaults<- function( x, cov.function,
                                      cov.args,
                                      cov.params.start,
                                      parGrid,
                                      mKrig.args,
                                  extraArgs = NULL,
                                      gridN = 5,
                        collapseFixedEffect = TRUE,
                                      verbose=FALSE
                        )
{
  
  ## convenient defaults for GP fitting.
  ## and also sort what starting parameter values are provided
  # this code runs on  by perhaps it is useful to see all the defualts and 
  # logic in one place.
  # Note the stranger device below where mKrig.args is created and amended
  #  mKrig is the basic computational function for evaluating the likleihood and 
  # setting up  the Kriging predictions. 
  #
  # aRange and lambda  are  handled specially because are almost always 
  # estimated and this will simplify the call in  this top level function 
  #
  ###########################################
  ## Set some convenient default choices for a 
  ## stationary covariance function 
  ###########################################
  if( is.null( cov.function)){
    cov.function <- 'stationary.cov'
    if( is.null(cov.args) ){
      cov.args<- list()
    }
    
    if( is.null(cov.args$Covariance )&is.null(extraArgs$Covariance )){
      cov.args$Covariance<- "Matern"
      if( is.null(cov.args$smoothness ) 
           & is.null(cov.params.start$smoothness ) 
           & is.null(parGrid$smoothness) ){
        cov.args$smoothness<- 1.0
      }
    }
  } 
  ###########################################
  ## Set some convenient default choices for a 
  ## thin plate spline  
  ###########################################
  if( cov.function=='Tps.cov'){
  # determine cardinal points if not included in
  # cov.args
    dimX<- ncol( x)
  if( is.null( mKrig.args$m)){
    # m should satisfy  2*m-dimX >0
    mMin<- max(c(2, ceiling(dimX/2 + 0.1))) 
   
    mKrig.args<- list( mKrig.args, m=mMin )
  }
    
#  
  if( is.null(cov.args)){
    cov.args<- list()
  }
#     
  if( is.null(cov.args$cardinalX)){
    nterms <- choose((mKrig.args$m + dimX - 1), dimX)
    cardinalX<- cover.design(x, nterms, num.nn = 50 )$design
    cov.args$cardinalX<- cardinalX
  }
    cov.args$aRange<- NA
    
  }
  
  ###########################################
  # overwrite the default choices if some are passed as ...
  #  (some R arcania!)
  ###########################################
  if( !is.null( extraArgs)){
    if(!is.null(cov.args)){
      ind<- match( names(cov.args), names(extraArgs) ) 
      cov.args <- c( cov.args[is.na(ind)], (extraArgs) )
    }
    else{
      cov.args <- list(extraArgs)
    }
  }
 ###########################################
 # check for duplicate arguments in starting values and fixed values
 ###########################################
  covArgsNames <- names(cov.args)
  covStartNames<-names(cov.params.start)
  covParGridNames<- names( parGrid)
  #print( covParGridNames)
  if( length( intersect( covArgsNames,covStartNames))>0){
    cat("A problem with duplicate  parameters:", fill=TRUE)
    cat("Names cov.args:", fill=TRUE)
    print(covArgsNames)
    cat("Names cov.params.start :", fill=TRUE)
    print(covStartNames)
    stop("parameters must either have starting values ( in cov.params.start list) 
     or be specified as a covariance function argument (in cov.args list) ")
  }
  
  
  if( verbose){
    cat("Updated and passed cov.args", fill=TRUE)
    print( cov.args)
  }
  
  
  ########################################### 
  # Some logic to figure out how do MLE search over lambda and aRange
  ###########################################
  noLambda<- is.null( cov.args$lambda) & is.null(cov.params.start$lambda)
  noARange<- is.null( cov.args$aRange) & is.null(cov.params.start$aRange)
  makeDefaultGrid<- (noLambda | noARange) & is.null(parGrid)
# easy default search grid if lambda and/or aRange ahave not been specified
  if( makeDefaultGrid ){
  if( noLambda){
    lGrid<- 10**seq( -4, .5, length.out= gridN)
  }
  if( noARange){
      minX<- apply( x, 2, min)
      maxX<- apply( x, 2, max)
      xCorners<- rbind( minX,
                        maxX)
      if( is.null( cov.args$Distance)){
        dMax<-rdist( rbind(xCorners[1,]), rbind(xCorners[2,]))
        
      }
      else{
        dMax<- do.call(cov.args$Distance, list(
                         x1= rbind(xCorners[1,]),
                         x2= rbind(xCorners[2,]))
                       )
      }
      dMax<- c( dMax)
        aGrid<- seq( .1*dMax, .7*dMax, length.out= gridN)
  }
 # now create parGrid   
      if( noLambda & !noARange){
        parGrid<- data.frame( lambda= lGrid)
      }
      if( noLambda & noARange){
        parGrid<- expand.grid( lambda= lGrid, aRange = aGrid)
      }
      if( !noLambda & noARange){
        parGrid<- data.frame( aRange= aGrid)
      }
  }
  ########################################### 
  # Identify the Cases 0 - 4 to set defaults
  ########################################### 
  
  # CASE 0 is to evaluate at fixed lambda and aRange
  # and there are no other parameters to optimize over.
  
  if( !is.null( cov.args$lambda) & 
      !is.null( cov.args$aRange) &
       is.null( cov.params.start) 
       ){
    CASE<- 0
  }
  
  #CASE 1 is to find MLEs using starting values provided a grid has not been 
  # supplied for an initial grid search.
  
  if( !is.null(cov.params.start) & is.null(parGrid) ){
    CASE<- 1
  }
 
  if( !is.null(parGrid) ){
    
    CASE<- 2
  }

########################################### 
# Messing with mKrig
########################################### 
# Determine linear fixed model if not specified and add in how to find fixed part.
# collapseFixedEffect is important enough where it is handled at this level.
#
  
  
  if( is.null(mKrig.args)){
    
    mKrig.args<- list( m=2, collapseFixedEffect=collapseFixedEffect)
   
  }
  else{
    if( all(names( mKrig.args)!= "collapseFixedEffect")){
      
      mKrig.args<- c( mKrig.args, 
                      list(collapseFixedEffect= collapseFixedEffect))
    }
  }
  
# don't find effective df for optimization -- this would add extra computation that is not 
# needed 
  if( is.null(mKrig.args$find.trA) ){
    if( (CASE >=3)){
    mKrig.args<- c( mKrig.args, list(find.trA = FALSE))
  }
  else{
    mKrig.args<- c( mKrig.args, list(find.trA = TRUE))
  }
  }
  
  out<- 
    list(  
        cov.function = cov.function,
            cov.args = cov.args,
          mKrig.args = mKrig.args, 
                CASE = CASE,
             parGrid = parGrid
        )
  
 
 
  return(
         out
          )
}
