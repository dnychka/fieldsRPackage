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

mKrigFastPredict <- function(object, grid.list, ynew = NULL,
                          derivative = 0, Z = NULL, drop.Z = FALSE,
                          collapseFixedEffect = object$collapseFixedEffect, 
                          np=4,
                          ...) {
  # the main reason to pass new args to the covariance is to increase
  # the temp space size for sparse multiplications
  # other optional arguments that typically describe the covariance function 
  # from mKrig are passed along in the list object$args
  cov.args <- list(...)
  xObs<- object$x
  
  nx<- length(grid.list$x )
  ny<- length(grid.list$y )
  
  if (!is.null(ynew)) {
    coef.hold <- mKrig.coef(object, ynew,
                            collapseFixedEffect=collapseFixedEffect)
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
    xnew<- make.surface.grid( grid.list)
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
  
  
  # syntax is the name of the function and then a list with
  # all the arguments. This allows for different covariance functions
  # that have been passed as their name.
  
  # enlarge the grid if needed so that obs have np grid point points on all margins. 
         
  if( (min(xObs[,1]) < grid.list$x[1]) | (max(xObs[,1]) > grid.list$x[nx] ) ) {
    stop( "x obs locations can not be outside the grid ")
  }
  if( (min(xObs[,2]) < grid.list$y[1]) | (max(xObs[,2]) > grid.list$y[ny])  ){
    stop( "y obs locations can not be outside the grid ")
  }
  
  # adjust grid if needed to include a margin of np+1 grid points beyond xObs
   
  marginInfo<- addMarginsGridList(xObs, grid.list, np)

  # these are the slightly larger grids by adding margins.
  gridListNew<-  marginInfo$gridListNew
  nxNew<- length(gridListNew$x )
  nyNew<- length(gridListNew$y )
  
  if (derivative == 0) {
    offGridObject<- offGridWeights( xObs,
                                    gridListNew,
                                    mKrigObject = object,
                                    np=np,
                                    giveWarnings = TRUE
    )
    cov.obj<- stationary.image.cov( setup=TRUE, 
                                    grid=gridListNew,
                                    cov.args= object$args)
    c.coefWghts<-  colSums( diag.spam( c(object$c.coef) ) %*% 
                       offGridObject$B )
    c.coefWghts<- matrix( c.coefWghts, nxNew, nyNew )
    temp2<-  stationary.image.cov(Y=c.coefWghts, cov.obj=cov.obj)
    # cut down the size of temp2 trimming off margins. 
    temp2<- temp2[marginInfo$indX, marginInfo$indY]
  }
  else {
    stop("Derivatives not supported with fast prediction method")
  }
  # add two parts together and coerce to vector
  if( object$nt>0){
    return((matrix(temp1,nx,ny) + temp2))
  }
  else{
    # only the spatial part fixed part is absent. 
    return(  temp2)
  }
}
