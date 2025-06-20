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

mKrigFastPredict <- function(object, gridList, ynew = NULL,
                          derivative = 0, Z = NULL, drop.Z = FALSE,
                          NNSize=4, setupObject= NULL,
                          giveWarnings=TRUE,
                          verbose=FALSE) 
                           {
  
  #NOTE: covariance model is specified by the arguments in object$args
  # cov.args <- c( object$args, list(...) )
  # For convenience the Z covariates are already assumed to be 
  # in the unrolled form. But this may be awkward if this 
  # function is called directly 
  # See the code in predictSurface.mKrig for details. E.g. unrollZGrid
  
  if (derivative != 0) {
    stop("Derivatives not supported with fast prediction method")
  }
                               
 if( ncol(object$c.coef)>1 ){
     stop("Replicated fields currently not supported for fast predict.")
 }
                               
  names( gridList)<- c("x","y")
  
  NNSize                            
  xObs<- object$x
  
  nx<- length(gridList$x )
  ny<- length(gridList$y )
  
  if (!is.null(ynew)) {
    coef.hold <- mKrig.coef(object, ynew,
                            collapseFixedEffect=TRUE)
    c.coef <- coef.hold$c.coef
    beta <- coef.hold$beta
  }
  else {
    c.coef <- object$c.coef
    beta <- object$beta
  }
  # fixed part of the model this a polynomial of degree m-1
  # Tmatrix <- fields.mkpoly(xnew, m=object$m)
  # only do this if nt>0, i.e. there is a fixed part.
  #
  if (!drop.Z & (object$nZ > 0) & (derivative >0) ) {
    stop("derivative not supported with Z covariate included
         use drop.Z = FALSE to omit Z ")
  }
  if( object$nt>0){
    xnew<- make.surface.grid( gridList)
    if (derivative == 0) {
      if (drop.Z | object$nZ == 0) {
        # just evaluate polynomial and not the Z covariate
        temp1 <- fields.mkpoly(xnew, m = object$m) %*% 
          beta[object$ind.drift, ]
      }
      else {
        if( nrow( xnew) != nrow(as.matrix(Z)) ){
          stop(paste("number of rows of covariate Z",
                     nrow(as.matrix(Z)), 
                     " is not the same as the number of locations",
                     nrow( xnew) )
          )
        }
        temp0 <-  cbind(fields.mkpoly(xnew, m = object$m),as.matrix(Z)) 
        temp1 <- temp0 %*% beta
      }
    }
    else {
      temp1 <- fields.derivative.poly(xnew, m = object$m,
                                      beta[object$ind.drift,])
    }
   
  }  
  # add nonparametric part. Covariance basis functions
  # times coefficients.

  # enlarge the  evaulation grid if needed so that obs have NNSize grid point points on all margins. 
  
  if( (min(xObs[,1]) < gridList$x[1]) | (max(xObs[,1]) > gridList$x[nx] ) ) {
    stop( "x obs locations can not be outside the grid ")
  }
  if( (min(xObs[,2]) < gridList$y[1]) | (max(xObs[,2]) > gridList$y[ny])  ){
    stop( "y obs locations can not be outside the grid ")
  }
 
  # adjust grid if needed to include a margin of NNSize+1 grid points beyond xObs
  # these are the slightly larger grids by adding margins.
  # also create sparse matrices. 
     
  if( is.null(setupObject) ){
  setupObject<- mKrigFastPredictSetup(object, 
                                      gridList = gridList, 
                                        NNSize = NNSize,
                                  giveWarnings = giveWarnings,
                                       verbose = verbose)
                                                 
  }
  # the expanded grid to include extra neighbors 
  gridListNew<-  setupObject$marginInfo$gridListNew
  nxNew<- length(gridListNew$x )
  nyNew<- length(gridListNew$y )
  # indX and indY are the subset of indices that match gridList
  # ( gridListNew contains the grids in gridList)
  indX<- setupObject$marginInfo$indX
  indY<- setupObject$marginInfo$indY
  # 
  if( ((indX[2]- indX[1] + 1)!= nx)| ((indY[2]- indY[1] + 1)!= ny)) {
    cat(" indX, nx") 
    print( c(indX, nx) )
    cat(" indY, ny")
    print( c(indY, ny) )
    stop("mismatch between  subset of larger grid and gridList passed")
  }
  # funky  coefficients on the grid that stand in for the actual ones.
    c.coefStar<-  colSums( diag.spam( c(object$c.coef) ) %*% 
                              setupObject$offGridObject$B )
    # reshape as an image. 
    c.coefStar<- matrix( c.coefWghts, nxNew, nyNew )
    #
    # fast multiplication of covariances on the grid with
    # the coefficients, cStar, on the grid (via FFT)
    # note that there can be lots of zeroes in coefWghts since the 
    # approximation is local. 
     
    temp2<-  stationary.image.cov( Y=c.coefStar, cov.obj=setupObject$cov.obj)
    
    # cut down the size of temp2 trimming off margins using to approximate
    # exact covariance kernel. 
    temp2<- temp2[ indX[1]:indX[2], indY[1]:indY[2] ]
  
  # add fixed part and spatial  parts together and coerce to matrix
  if( object$nt>0){
    return( (matrix(temp1,nx,ny) + temp2) )
  }
  else{
    # return only the spatial part  because the fixed part is absent. 
    return(  temp2)
  }
}
