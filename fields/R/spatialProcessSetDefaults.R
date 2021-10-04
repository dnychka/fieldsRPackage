#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2021 Colorado School of Mines
# 1500 Illinois St., Golden, CO 80401
# Contact: Douglas Nychka,  douglasnychka@gmail.edu,
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
                                      extraArgs=NULL,
                                      gridN=5,
                                      verbose=FALSE)
{
  
  ## convenient defaults for GP fitting.
  ## and also sort what starting parameter values are provided
  #
  # aRange and lambda  are  handled specially because are almost always 
  # estimated and this will simplify the call in  this top level function 
  #
  if( is.null( cov.function)){
    cov.function <- 'stationary.cov'
    if( is.null(cov.args) ){
      cov.args<- list()
    }
    if( is.null(cov.args$Covariance )){
      cov.args$Covariance<- "Matern"
      if( is.null(cov.args$smoothness ) 
           & is.null(cov.params.start$smoothness ) 
           & is.null(parGrid$smoothness) ){
        cov.args$smoothness<- 1.0
      }
    }
  } 
  
  # overwrite the default choices if some are passed as ...
  #  (some R arcania!)
  
  if( !is.null( extraArgs)){
    if(!is.null(cov.args)){
      ind<- match( names(cov.args), names(extraArgs) ) 
      cov.args <- c( cov.args[is.na(ind)], (extraArgs) )
    }
    else{
      cov.args <- list(extraArgs)
    }
  }
  
 # check for duplicate arguments in starting values and fixed values
  
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
  
  
  
# linear fixed model if not specified. 
  if( is.null(mKrig.args)){
    mKrig.args<- list( m=2)
  }
  
# don't find eff df for optimization  
  if( is.null(mKrig.args$find.trA) ){
    if( (CASE >=3)){
    mKrig.args<- c( mKrig.args, list(find.trA = FALSE))
  }
  else{
    mKrig.args<- c( mKrig.args, list(find.trA = TRUE))
  }
  }
  
  #
  # tuck in starting value for lambda if missing
  # 
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
