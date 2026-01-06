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

"predictSurface.Krig" <- function(object, grid.list = NULL, 
       extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
       xy = c(1,2),  verbose = FALSE,
       XMatGrid=NULL, drop.XMat= FALSE, just.fixed=FALSE,  ...) {
  
      if( is.null(XMatGrid) & !drop.XMat & (!is.null(object$XMat)) ) {
      stop("Need to specify covariate (XMat) values or set drop.XMat==TRUE")
    }
# create a default grid if it is not passed    
    if (is.null(grid.list)) {
    # in more than 2-D 
    # default is 80X80 grid on first two variables
    # rest are set to median value of the x's
        grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
            xy = xy)
    }
# do some checks on XMatgrid and also reshape as a matrix
# rows index grid locations and columns  are the covariates
# (as XMat in predict).
# if XMatGrid is NULL just returns that back 
    XMat<- unrollXMatGrid( grid.list, XMatGrid) 
    xg <- make.surface.grid(grid.list)
# NOTE: the predict function called will need to do some internal  the checks
# whether the evaluation of a large number of grid points (xg)  makes sense.
  if( verbose){
    print( dim( xg))
    print( nrow( xg))
    print( drop.XMat)
    print( dim( XMat))
  }
# if extrapolate is FALSE set all values outside convex hull to NA
  if (!extrap) {
    if( is.null( object$x)){
        stop("need and x matrix in object")
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
    
  out<- rep( NA, nrow(xg))
# here is the heavy lifting    
    out[indexGood] <-  predict(object, x=xg[indexGood,], XMat=XMat[indexGood,], drop.XMat= drop.XMat,   
                     just.fixed=just.fixed, ...)
# reshape as list with x, y and XMat components    
    out <-  as.surface( xg, out )
    #
    #
    return(out)
}
