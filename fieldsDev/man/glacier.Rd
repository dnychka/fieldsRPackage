%#
%# fields  is a package for analysis of spatial data written for
%# the R software environment.
%# Copyright (C) 2024 Colorado School of Mines
%# 1500 Illinois St., Golden, CO 80401
%# Contact: Douglas Nychka,  douglasnychka@gmail.edu,
%#
%# This program is free software; you can redistribute it and/or modify
%# it under the terms of the GNU General Public License as published by
%# the Free Software Foundation; either version 2 of the License, or
%# (at your option) any later version.
%# This program is distributed in the hope that it will be useful,
%# but WITHOUT ANY WARRANTY; without even the implied warranty of
%# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%# GNU General Public License for more details.
%#
%# You should have received a copy of the GNU General Public License
%# along with the R software environment if not, write to the Free Software
%# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%# or see http://www.r-project.org/Licenses/GPL-2
%##END HEADER
%##END HEADER

\name{glacier}
\alias{glacier}
\docType{data}
\title{Franke's Glacier Elevation Data}
\description{
A moderate size (about 8400 locations) spatial dataset that is well-known in the applied mathematics approximation literature for testing interpolation methods.
}
\usage{data(glacier)}
\format{
  The format of \code{glacier} is a list with two components:
\describe{
\item{loc:}{8338x2 matrix of the  locations (meters??).} 
\item{y:}{A vector of elevations (meters ??).} 
}
}
\details{ 
This  data set  appears in papers that develop interpolation methods for
scattered  data and serves as an interesting bridge to the examples 
in applied math that develop radial basis function surface fitting.
The data was originally used by R. Franke. 

Unfortunately at this time we can not find any background on where
these data were collected or indeed even the location of this glacier.
However, it is an interesting data set in that it appears that
the elevations are reported along lnes of equal elevation, i.e. 
contours, perhaps from a digitization of a
topographic map or survey. It is important to estimate the
surface in a way that the artifacts from discretization are not
present. In the example below the compactly supported kernel
interpolation still has some artifacts. 

The glacier data set is available at this website
\url{https://oleg-davydov.de/scat_data.html}

The examples below are useful for comparing different
approximations
to a Gaussian spatial process estimate for the elevation surface.
Of
course in using a stationary covariance ( e.g. the Matern or
Wendland)
these are also radial basis smoothing or interpolation of the
data. 
}
\examples{
data( glacier )
# EDA for raw obs:

bubblePlot( glacier$loc, glacier$y, highlight=FALSE, size=.5)

# identifying contour levels. Note this is reported at regular levels
# (Every 25m ???)

table( glacier$y)



# find sigma and rho by maximum likelihood 
# for a fixed range
#  the default is the Wendland covariance with k=2
# See help(Wendland)

# this takes about 5 minutes
# macbook pro Quad-Core Intel Core i5 8 GB

#options(spam.nearestdistnnz=c(5e7,1e3))
#system.time( 
# obj0<- fastTps(glacier$loc, glacier$y, 
#                       theta=2,
#                      profileLambda=TRUE) 
#)
 

# set.panel(2,2)
# plot( obj0)
# set.panel()

# just evaluate at MLE
# reset default matrix size that the spam pacakge will use.

\dontrun{

options(spam.nearestdistnnz=c(5e7,1e3))
system.time( obj1<- 
               fastTps(glacier$loc, glacier$y, 
                       theta=2,
                       lambda= 7.58e-5
                        ) 
)

system.time(
look1<- predictSurface( obj1, nx=150, ny=150)
)

imagePlot( look1)


system.time(
out<- simLocal.spatialProcess(obj1, M=3, nx=150, ny=150)
)
set.panel( 2,2)
imagePlot( look1)
zlim<- range( out$z, na.rm=TRUE)
for( k in 1:3){
imagePlot(out$x, out$y, out$z[,,k], zlim=zlim)
}

# near interpolation surface using Matern smoothness .5 
 system.time( 
 obj2<- spatialProcess(glacier$loc, glacier$y,
                          aRange = 1.5, 
                          lambda = 1e-5,
                          smoothness = .5)
 )
 
system.time(
out<- simLocal.spatialProcess(obj2, M=3, nx=150, ny=150,
fast=TRUE)
)

set.panel( 2,2)
imagePlot( look1)
zlim<- range( out$z, na.rm=TRUE)
for( k in 1:3){
imagePlot(out$x, out$y, out$z[,,k], zlim=zlim)
}

% test out fast predict algorithm verses exact 
% note speedup of about 15 times 
system.time(
look2<- predictSurface.mKrig( obj2, nx=150, ny=150,
                fast=TRUE, NNSize=5)
)

system.time(
look2B<- predictSurface( obj2, nx=150, ny=150,
                fast=FALSE)
)

err<- c((look2$z - look2B$z)/look2B$z)
stats( log10( abs(err) ) )

# some error plots ( percent relative error)
imagePlot(look2$x, look2$y, 100*(look2$z - look2B$z)/look2B$z  )

imagePlot(look2$x, look2$y, 100*(look1$z - look2B$z)/look2B$z  )

} % end do not run
} % end examples
\keyword{datasets}
