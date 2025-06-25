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

\name{predictSurface}
\alias{predictSurface}
\alias{predictSurface.default}
\alias{predictSurface.mKrig}
\alias{predictSurface.Krig}
\alias{predictSurface.fastTps}
\alias{predictSurfaceSE}
\alias{predictSurfaceSE.default}
\alias{mKrigFastPredict}

\title{
  Evaluates a fitted function or the prediction error as a surface that is suitable for plotting with 
the image, persp, or contour functions.
}
\description{
Evaluates a a fitted model or the prediction error on a 2-D grid keeping any other variables constant.
The resulting object is suitable for use with functions for viewing 3-d
surfaces such as image, imagePlot and contour. 
}
\usage{
\method{predictSurface}{default}(object, grid.list = NULL, 
                     extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
                        xy = c(1,2),  verbose = FALSE, ...)

\method{predictSurface}{fastTps}(object, gridList = NULL, 
       extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
       xy = c(1,2),  verbose = FALSE, ...)
       
\method{predictSurface}{Krig}(object, grid.list = NULL, extrap = FALSE, chull.mask = NA, 
nx = 80, ny = 80, xy = c(1, 2), verbose = FALSE, ZGrid = NULL, 
    drop.Z = FALSE, just.fixed=FALSE,  ...)
    
\method{predictSurface}{mKrig}(
object, gridList = NULL, grid.list = NULL, ynew =
                 NULL, extrap = FALSE, chull.mask = NA, nx = 80, ny =
                 80, xy = c(1, 2), verbose = FALSE, ZGrid = NULL,
                 drop.Z = FALSE, just.fixed = FALSE, fast = FALSE,
                 NNSize = 4, setupObject = NULL, giveWarnings = FALSE,
                 derivative = 0, ...)
		 
mKrigFastPredict(object, gridList, ynew = NULL, derivative = 0, Z =
                 NULL, drop.Z = FALSE, NNSize = 4, setupObject = NULL,
                 giveWarnings = TRUE, verbose = FALSE)
			  
\method{predictSurfaceSE}{default}(object, grid.list = NULL, extrap = FALSE, chull.mask =
                 NA, nx = 80, ny = 80, xy = c(1, 2), verbose = FALSE,
                 ZGrid = NULL, just.fixed = FALSE, ...)
}

\arguments{
\item{object}{
An object from fitting a function to data. In fields this is usually a
Krig, mKrig, or fastTps object. 
}
\item{gridList}{
A list with as many components as variables describing the surface. 
All components should have a single value except the two that give the 
grid points for evaluation. If the matrix or data frame has column names,  
these must appear in the grid list. See the grid.list help file for more
details. If this is omitted and the fit just depends on two variables the
grid will be made from the ranges of the observed variables. 
(See the function \code{fields.x.to.grid}.)
 
}
\item{grid.list}{Alternative to the  \code{gridList} argument. }

\item{giveWarnings}{If TRUE will warn when more than one observation is in a grid box. }

\item{extrap}{
 Extrapolation beyond the range of the data. If \code{FALSE} (the 
default) the predictions will be restricted to the convex hull of the observed 
data or the convex hull defined from the points from the argument chull.mask. 
This function may be slightly faster if this logical is set to 
\code{TRUE} to avoid checking the grid points for membership in the 
convex hull. For more complicated masking a low level creation of a bounding 
polygon and testing for membership with \code{in.poly} may be useful. 

}

\item{chull.mask}{
Whether to restrict the fitted surface to be on a convex hull, NA's
are assigned to values outside the
convex hull. chull.mask should be a sequence of points defining a convex
hull. Default is to form the convex hull from the observations if this
argument is missing (and extrap is false).  
}
  

\item{nx}{Number of grid points in X axis. }

\item{ny}{Number of grid points in Y axis. }

\item{xy}{ A two element vector giving the positions for the "X" and "Y"
variables for the surface. The positions refer to the columns of the x 
matrix used to define the multidimensional surface. This argument is 
provided in lieu of generating the grid list. If a 4 dimensional surface
is fit to data then \code{ xy= c(2,4)} will evaluate a surface using the 
second and fourth variables with  variables 1 and 3 fixed at their median 
values. NOTE: this argument is ignored if a grid.list argument is 
passed. }

\item{drop.Z}{If TRUE the fixed part of model depending on covariates is omitted.}

\item{just.fixed}{If TRUE the nonparametric surface is omitted.}

\item{fast}{If TRUE approximate predictions for stationary models are made using the FFT. For large grids( e.g. nx, ny > 200) this can be substantially faster and still accurate to several decimal places. }

\item{NNSize}{Order of nearest neighborhood used for fast prediction.  The default,\code{NSize = 5}, means an 11X11=121 set of grid points/covariance kernels are used to approximate the off-grid covariance kernel. }

\item{setupObject}{The object created explicitly using
\code{\link{mKrigFastPredictSetup}}. Useful for predicting multple surfaces with the same observation locations.
If NULL the setup is done each time the function is called. (See \code{\link{approximateCovariance2D}}.) }

\item{derivative}{Predict the  estimated derivatives of order \code{derivative}.}

 \item{ynew}{
 New data to use to refit the spatial model. Locations must be the same but if so this is efficient because the matrix decompositions are reused. 
 }
  
 \item{\dots}{
Any other arguments to pass to the predict function associated with the fit
object. 
Some of the usual arguments for several of the fields fitted objects include:
\describe{

\item{ynew}{ New values of y used to reestimate the surface.}

\item{Z}{A matrix of covariates for the fixed part of model.}
}
}

\item{ZGrid}{An array  or list form of covariates to use for
	 prediction. This must match the same dimensions from the 
\code{grid.list} / \code{gridList} argument.  

If ZGrid is an array then the first two indices are the x and y
 locations in the 
grid. The third index, if present, indexes the covariates. e.g. For
 evaluation on 
a 10X15 grid and with 2 covariates. \code{ dim( ZGrid) == c(10,15, 2)}.
If ZGrid is a list then the components x and y shold match those of grid list and
the z component follows the shape described above for the no list 
case. 
}

\item{Z}{The covariates for the grid unrolled as a matrix. Columns index
the variables and rows index the grid locations. E.g. For
 evaluation on  a 10X15 grid and with 2 covariates. \code{ dim( ZGrid) == c(10,15, 2)}.  and so \code{  dim( Z) = c(150, 2)} and  
 \code{ Z[,1] <- c( ZGrid[,,1])}
 }

 \item{verbose}{If TRUE prints out some imtermediate results for debugging.}

}

\value{
The usual list components for making image, contour, and perspective plots
(x,y,z) along with labels for the x and y variables. For
\code{predictSurface.derivative} the component \code{z} is a three
dimensional array with values( \code{nx}, \code{ny}, 2 ) 
 } 
\details{ These function evaluate the spatial process or thin plate spline estimates on a regualr grid of points 
The  grid can be specified using the grid.list/ gridList  information or just the sizes.

For the standard Krig and mKrig versions the steps are to create a matrix of locations the represent the grid, 
 call the predict function for the object with these
points and also adding any extra arguments passed in the ... section,
and then reform the results as a surface object (as.surface). To
determine the what parts of the prediction grid are in the convex hull
of the data the function \code{in.poly} is used. The argument
inflation in this function is used to include a small margin around
the outside of the polygon so that point on convex hull are
included. This potentially confusing modification is to prevent
excluding grid points that fall exactly on the ranges of the
data. Also note that as written there is no computational savings for
evaluting only the convex subset compared to the full grid.

For the "fast" option a stationary covariance function and  resulting surface
estimate is approximated by the covariance kernel restricted to the grid
locations. In this way the approximate problem becomes a 2-d convolution.
The evaluation of the approximate prediction surface uses a fast Fourier
transform to compute the predicted values at the grid locations.


The nearest
neighbor argument \code{NNSize} controls the number of covariance kernels
only evalauted at grid location used
to approximate a  covariance function at an off-grid location. We have
found good results  with \code{NNSize=5}. 


\code{predictSurface.fastTps} is a specific version ( m=2, and k=2) of
Kriging with a compact covariance kernel (Wendland). 
that can be much more efficient because it takes advantage of a low
level FORTRAN call to evaluate the  covariance function. Use
\code{predictSurface} or \code{predict} for other choices of m and k.


\code{predictSurface.Krig} is designed to also include covariates for the fixed in terms of grids.

\code{predictSurface.mKrig} Similar in function to the Krig prediction function but it more efficient using the \code{mKrig} fit object.

\code{mKrigFastpredict} Although this function might be called at the top is it easier to use through the wrapper,  \code{predictSurface.mKrig} and \code{fast=TRUE}. 
 
NOTE: \code{predict.surface} has been depreciated and just prints out
a warning when called. 

 }
\seealso{
Tps, Krig, predict, grid.list, make.surface.grid, as.surface, surface, 
in.poly
}
\examples{


data( ozone2)

  x<- ozone2$lon.lat 
  y<- ozone2$y[16,]
  
  obj<- Tps( x,y)
  
  # or try the alternative model: 
  #  obj<- spatialProcess(x,y)
  
  fit<- predictSurface( obj, nx=40, ny=40)
  imagePlot( fit)

# predicting a 2d surface holding other variables fixed.

fit<- Tps( BD[,1:4], BD$lnya)  # fit surface to data 

# evaluate fitted surface for  first two 
# variables holding other two fixed at median values

out.p<- predictSurface(fit)
surface(out.p, type="C") 

#
# plot surface for second and fourth variables 
# on specific grid. 

glist<- list( KCL=29.77, MgCl2= seq(3,7,,25), KPO4=32.13, 
                     dNTP=seq( 250,1500,,25))

out.p<- predictSurface(fit, glist)
surface(out.p, type="C")

out.p<- predictSurfaceSE(fit, glist)
surface(out.p, type="C")

## a test of the fast prediction algorithm for use with 
# mKrig/spatialProcess objects. 
# for explanation of the fit see the MLEs from 
# fitObject<- spatialProcess( x,y)
#  fitObject$summary

\dontrun{
data(NorthAmericanRainfall)

x<- cbind(NorthAmericanRainfall$longitude,     
          NorthAmericanRainfall$latitude)
y<- log10(NorthAmericanRainfall$precip)

mKrigObject<- mKrig( x,y,
                    lambda=.024,
                    cov.args= list(     aRange= 5.87,
                                    Covariance="Matern",
                                    smoothness=1.0),
                                        sigma2=.157
                                 )
gridList<- list( x = seq(-134, -51, length.out = 100),
                 y = seq( 23,   57, length.out = 100))                          
                 
# exact prediction 
system.time( 
gHat<- predictSurface( mKrigObject, gridList)
)

# aproximate 
system.time( 
gHat1<- predictSurface( mKrigObject, gridList,
                                fast = TRUE)
)
                                
# don't worry about the warning ...
# just indicates some observation locations are located 
# in the same grid box.

# approximation error omitting the NAs from outside the convex hull                        
 stats( log10(abs(c(gHat$z - gHat1$z))) )
 
 image.plot(gHat$x, gHat$y,  (gHat$z - gHat1$z) )
 points( x,  pch=".", cex=.5)
 world( add=TRUE )
 
 }
 
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
