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

"predictSurface.mKrig" <- function(object, gridList=NULL, grid.list=NULL, 
        ynew = NULL,                          
       extrap = FALSE, chull.mask = NA,
       nx = 80, ny = 80,
       xy = c(1,2),  verbose = FALSE,
       ZGrid=NULL, drop.Z= FALSE, just.fixed=FALSE, 
       fast=FALSE, NNSize=4, giveWarnings=FALSE,
       derivative=0, ...) {
#
      if( is.null(ZGrid) & !drop.Z & (!is.null(object$Z)) ) {
      stop("Need to specify covariate (Z) values or set drop.Z==TRUE")
      }
  
# grid.list is old syntax for fields gridList preferred 
  if (!is.null(grid.list)){ gridList<-  grid.list}
  
# create a default grid if it is not passed    
    if (is.null(gridList)) {
    # in more than 2-D 
    # default is 80X80 grid on first two variables
    # rest are set to median value of the x's
        gridList <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
            xy = xy)
       
    }
  
  #print( gridList)
# do some checks on Zgrid and also reshape as a matrix
# rows index grid locations and columns  are the covariates
# (as Z in predict).
# if ZGrid is NULL just returns NULL back  ...
    Z<- unrollZGrid( gridList, ZGrid) 
    xg <- make.surface.grid(gridList)
# NOTE: the predict function called will need to do some internal checks
# whether the evaluation of a large number of grid points (xg)  makes sense.
  if( verbose){
    print( dim( xg))
    print( nrow( xg))
    print( drop.Z)
    print( dim( Z))
  }
    if( nrow(xg) > 5e5){
      warning("number of grid points is large for exact prediction
              consider approximate prediction using fast==TRUE")
    }
    out<- rep( NA, nrow(xg))  
# if extrapolate is FALSE set all values outside convex hull to NA
    if (!extrap) {
      if( is.null( object$x)){
        stop("need an x matrix in object")
      }
      if (is.na(chull.mask)) {
        chull.mask <- unique.matrix(object$x[, xy])
      }
      indexGood<- in.poly(xg[, xy], xp = chull.mask, convex.hull = TRUE)
      if( verbose){
        print( sum( indexGood) )
      }
    }
    else{
      indexGood<- rep( TRUE, nrow( xg))
    }
      
    if(!fast){
# here is the heavy lifting 
    out[indexGood] <-  predict.mKrig(object, xnew=xg[indexGood,], ynew=ynew,
                               Z=Z[indexGood,], drop.Z= drop.Z, 
                               collapseFixedEffect = object$collapseFixedEffect,
                               just.fixed=just.fixed, ...)
  }
  else{
  # fast approximate method 
  # will always predict  to full grid.
    
  out<- mKrigFastPredict( object,
                          gridList= gridList, 
                          ynew = ynew,
                          derivative = derivative,
                          Z = Z, drop.Z = drop.Z,
                          NNSize=NNSize, 
                          giveWarnings=giveWarnings,
                          
                          ...
                           )
  # wipe out predictions outside convex hull of observations. 
  if(!extrap){
    out[ !indexGood]<- NA
  }
    
  }
# reshape as list with x, y and z components    
    out <-  as.surface( xg, out )
    #
    #
    return(out)
}
