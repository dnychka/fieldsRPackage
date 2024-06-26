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
"predict.surface" <- function(object, ...) {
    UseMethod("predict.surface")
}

predict.surface.default<- function(object,...){
   cat("predict.surface is now the function predictSurface")
 }

"predictSurface"<- function( object,...){
  UseMethod("predictSurface")
}

"predictSurface.default" <- function(object, grid.list = NULL, 
       extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
       xy = c(1,2),  verbose = FALSE, ...) {
    # NOTE: 
    # without grid.list
    # default is 80X80 grid on first two variables
    # rest are set to median value of x.
    if (is.null(grid.list)) {
        grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
            xy = xy)
    } 
    # here is the heavy lifting
    xg <- make.surface.grid(grid.list)
# NOTE: the specific predict function called will need to do the checks
# whether the evaluation of a large number of grid points makes sense. 
    out <-  as.surface( xg, predict(object, xg,...) )
    #
    # if extrapolate is FALSE set all values outside convex hull to NA
    if (!extrap) {
        if( is.null( object$x)){
          stop("need and x matrix in object")
        }
        if (is.na(chull.mask)) {
            chull.mask <- unique.matrix(object$x[, xy])
        }
        out$z[!in.poly(xg[, xy], xp = chull.mask, convex.hull = TRUE)] <- NA
    }
    #
    return(out)
}


"predictSurface.fastTps" <- function(object, gridList=NULL, 
       extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
       xy = c(1,2),  verbose = FALSE, ...) {
# NOTE:  See  predictSurface.default for comments
    if (is.null(gridList)) {
        gridList <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
            xy = xy)
    } 
#  in the case of fastTps and not great circle distance
#  pass the grid list instead of the locations of grid points
#  (see xg in predictSurface.default)
    if( object$args$Dist.args$method =="greatcircle"){
      out <-  predict(object, 
                      xnew= make.surface.grid(gridList), 
                      xy=xy, ...)
    }
    else{
    out <-  predict(object, grid.list=gridList, xy=xy, ...)
    }
# coerce to image format   
    out <-  as.surface(gridList, out )
    #
    # if extrapolate is FALSE set all values outside convex hull to NA
    if (!extrap) {
        if (is.na(chull.mask)) {
            chull.mask <- unique.matrix(object$x[, xy])
        }
        xg<- make.surface.grid( gridList)
        out$z[!in.poly(xg[, xy], xp = chull.mask, convex.hull = TRUE)] <- NA
    }
    #
    return(out)
}
