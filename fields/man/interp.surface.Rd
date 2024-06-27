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

\name{interp.surface}
\alias{interp.surface}
\alias{interp.surface.grid}
\alias{interp.surface.FFT}
\alias{fillGrid}
\title{
  Fast bilinear interpolator from a grid.
}
\description{
Uses bilinear weights to interpolate values on a rectangular
grid to arbitrary locations or to another grid.
}
\usage{
interp.surface(obj, loc)
interp.surface.grid(obj, grid.list)
interp.surface.FFT(obj, M)
fillGrid( gridList, M)
}
\arguments{
\item{obj}{
A list with components x,y, and z in the same style as used by contour,
persp, image etc. x and y are the X and Y grid values and z is a matrix
with the corresponding values of the surface. For the FFT method there must be an 
odd number of grid points. 
}
\item{loc}{
A matrix of (irregular) locations to interpolate. First column of loc
isthe X coordinates and second is the Y's. 
}
\item{grid.list}{ A list with components x and y 
describing the grid to interpolate. The grids do not need to be equally spaced.}

\item{gridList}{Just a newer name for a grid list object.}

\item{M}{ A muliple of the source grid to interpolate to create the target grid. E.g. if M=5 then a 10X20 source grid in \code{obj}
will give an interpolated grid of (10*5)X(20*5).
See the \code{fillGrid} function that creates the larger grid.}

}
\value{
 \strong{ interp.surface} An vector of interpolated values.  NA are returned for
regions of the  z matrix  that are NA and also for locations outside of the
range of the parent grid. 

\strong{ interp.surface.grid} Interpolated values using bilinear interpolation in the  list/image format  with comonents: \code{x}, \code{y}, \code{z}.

\strong{ interp.surface.FFT} Interpolated values using the FFT method in an image format with the grid refined by the factor \code{M}.

}
\details{

\strong{ interp.surface} Here is a brief explanation of
the interpolation:  Suppose that the location, (locx, locy)  lies in
between the first two grid points in both x an y. That is  locx is between
x1 and x2 and
locy is between y1 and y2.  Let ex= (l1-x1)/(x2-x1) ey= (l2-y1)/(y2-y1).
The
interpolant is

( 1-ex)(1-ey)*z11 + (1- ex)(ey)*z12 + ( ex)(1-ey)*z21 + ( ex)(ey)*z22  

Where the z's are the corresponding elements of the Z matrix. 

Note that bilinear interpolation can produce some artifacts related to
the grid and not reproduce higher behavior in the surface. For, example
the extrema of the interpolated surface will always be at the parent
grid locations.  There is nothing special about about interpolating to
another grid, this function just includes a \code{for} loop over one
dimension and a call to the function for irregular locations. It was
included in fields for convenience. since the grid format is so common. 

See also the akima package for fast interpolation from irrgeular locations. 
Many thanks to  Jean-Olivier Irisson for making this code more efficient and 
concise. 

\strong{ interp.surface.FFT} This version does the interpolation via the usual sin
/cosine basis and the FFT. This method only makes sense for interpolating to a
target grid that is a refinement of the source grid. E.g. if M=5 then a 10X20
source grid in \code{obj}
will give an interpolated grid of (10*5)X(20*5) = 50X200. endpoints of this grid match the endpoints of the source.

Note the FFT interpolation is
"C infinity " accurate which means in practice that this method will do well for
smooth fields. Sharp changes will induce the usual Gibbs oscillations around large
the changes. The interpolation is also peridic in both dimensions ( a torus) and so might give strange results for non-periodic fields. Pad the edges to mitigate this artifact. 

The algorithm in brief is 
1) FFT of the source image/matrix. 
2) Stuff this into corners of a larger matrix of zeroes and of the size of the refined grid.
3) Inverse FFT of stuffed matrix.

}
\seealso{
image.smooth, as.surface, as.image, imagePlot, image.plot fastTps
}
\examples{
#
# evaluate an image object  at a finer grid
# 

data( lennon)
# create an example in the right list format like image or contour or persp.
obj<- list( x= 1:21, y=1:21, z= lennon[ 201:221, 201:221])

set.seed( 123)
# lots of random points
N<- 500
loc<- cbind( runif(N)*20, runif(N)*20)
z.new<- interp.surface( obj, loc)
# compare the image with bilinear interpolation at scattered points
set.panel(2,2)
image.plot( obj)
quilt.plot( loc, z.new) 

# sample at 100X100 equally spaced points on a grid

grid.list<- list( x= seq( 1,20,,100), y=  seq( 1,20,,100))

interp.surface.grid( obj, grid.list)-> look

# this will give an error in the FFT version because 
# there are an even number of grid points for x.
# objTest<- list( x= 1:20, y=1:21, z= lennon[ 201:220, 201:221])
# look2<- interp.surface.FFT( objTest, M=10)

look2<- interp.surface.FFT( obj, M=20)-> look2
# take a look
set.panel(2,2)
image.plot( obj)
image.plot( look)
image.plot( look2)

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
