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
"fastTps" <- function(x, Y, m = NULL, p = NULL, aRange, 
                      lon.lat = FALSE, find.trA=FALSE,
                      REML=FALSE, theta=NULL,
                      ...) {
  # theta argument has been depreciated.
  if( !is.null( theta)){
    aRange<- theta
  }
  x <- as.matrix(x)
  d <- ncol(x)
  if (is.null(p)) {
    if (is.null(m)) {
      m <- max(c(2, ceiling(d/2 + 0.1)))
    }
    p <- (2 * m - d)
    if (p <= 0) {
      warning(" m is too small to satisfy thin plate spline crierion \n
                    you must have 2*m - dimension >0 \n
                    smoothness of Wendland set at k =2")
    }
  }
  # special arguments to send to the wendland covariance/taper function.
  # see nearest.dist for some explanation of 'method'
  method <- ifelse(!lon.lat, 
                   "euclidean", "greatcircle")
  cov.args <- list(        k = max(c(p,2)),
                           Dist.args = list(method=method),
                           aRange = aRange
  )
  
  object<-spatialProcess(x, Y, 
                         cov.function = "wendland.cov", 
                         mKrig.args = list( m = m), 
                         cov.args = cov.args,  REML=REML,
                         ... )
  
  object$call<- match.call()
  class(object) <- c( "fastTps", class(object))
  return( object)
}
