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
\name{offGridWeights}
\alias{offGridWeights}
\alias{addMarginsGridList}
\alias{mKrigFastPredictSetup}
\alias{offGridWeights1D}
\alias{offGridWeights2D}
\alias{findGridBox}

\title{
Utilities for fast spatial prediction. 
}
\description{
Based on a stationary Gaussian process model these functions support fast prediction onto
a grid using a sparse matrix approximation. They also allow for fast 
prediction to off-grid values (aka interpoltation) from an equally spaced rectangular grid and  using a spatial model. The sparsity comes about because only a fixed number of neighboring grid points
(NNSize) are used in the prediction.  The prediction variance for off-grid location is also give in the returned object. These function are used as the basis for approximate conditional simulation for large
spatial datasets asnd also for fast spatial prediction from irregular locations onto a grid.

}
\usage{
offGridWeights(s, gridList, NNSize = 2, mKrigObject = NULL, Covariance = NULL,
   covArgs = NULL, aRange = NULL, sigma2 = NULL, giveWarnings = TRUE,
   debug=FALSE)
   
   offGridWeights1D(s, gridList, np = 2, mKrigObject = NULL,
   Covariance = NULL,
   covArgs = NULL, aRange = NULL, sigma2 = NULL, giveWarnings = TRUE,
   debug=FALSE )
   
   offGridWeights2D(s, gridList, NNSize = 2, mKrigObject = NULL,
   Covariance = NULL,
   covArgs = NULL, aRange = NULL, sigma2 = NULL, giveWarnings = TRUE, 
  debug = FALSE,findCov = TRUE)
   
addMarginsGridList( xObs, gridList, NNSize)

findGridBox(s, gridList) 

mKrigFastPredictSetup(mKrigObject,
                          gridList, 
                          NNSize,
			  giveWarnings = TRUE)
}



\arguments{ 
 
\item{aRange}{
  The range parameter.}
\item{covArgs}{
  If \code{mKrigObject} is not specified  a list giving any additional arguments for the covariance function. }

  \item{Covariance}{
  The stationary covariance function  (taking pairwise distances as its
  first argument.) }
 
 \item{debug}{If TRUE returns intermediate calculations and structures for debugging and checking.}
 
\item{findCov}{If TRUE the prediction variances (via the Kriging
  model and assumed covariance function) for the offgrid locations are
  found. If FALSE just the prediction weights are returned. This option is a bit faster. 
 This switch only is implemented for the 2 D case since the extra  computational demands for 1 D are modest.}
 
 \item{giveWarnings}{If TRUE will warn if two or more observations
 are in the same grid box. See details below.}
 
  \item{gridList}{
  A list as the gridList format ( x and y components) that describes the
  rectagular grid. The grid must have at least \code{NNSize} extra grid points beyond
  the range of the points in \code{s}. See
  \code{\link{grid.list}}. 
}

\item{mKrigObject}{
  The output object (Aka a list with some specfic components.) from either mKrig or spatialProcess. This has the information about the covariance function used to do the Kriging. 
  The following items are coded in place of not supplying this object.  See the example below for more details.
}

\item{NNSize}{
  Number of nearest neighbor grid points to use for prediction.  \code{NNSize = 1}
  will use the 4 grid points that bound the offgrid point.  \code{NNSize = 2}
  will be a 4X4 subgrid that define a 3X3 set of boxes
  and the middle grid box containing the offgrid point.
  In general there will be  \code{(2*NNSize)^2} neighboring points used. 

}

  \item{s}{
  Off grid spatial locations
}

 \item{sigma2}{
  Marginal variance of the process. 
}
 
  \item{xObs}{
    Off grid spatial locations}
}

\details{
\strong{offGridWeights} This function creates the interpolation weights taking advantage of some
efficiency in the covariance function being stationary, use
of a fixed configuration of nearest neighbors, and Kriging
predictions from a rectangular grid. 

The returned matrix is in spam sparse matrix format. See
example below for the "one-liner" to make the prediction
once the weights are computed. Although created primarily
for (approximate) conditional simulation of a spatial process this
function is also useful for interpolating to off grid
locations from a rectangular field. It is also used for
fast, but approximate prediction for Kriging with a
stationary covariance. 

The function \code{offGridWeights} is a simple wrapper to
call either the 1D or 2D functions 

In most cases one would not use these approximations for a 1D problem.
However,
the 1D algorithm is included as a separate function for testing and also
because this is easier to read and understand the conversion between the
Kriging weights for each point and the sparse matrix encoding of them.


The interpolation errors are also computed based on the nearest neighbor
predictions. This is returned as a sparse matrix in the component SE. 
If all observations are in different grid boxes then \code{SE} is
diagonal
and agrees with the square root of the component
\code{predctionVariance} but
if multiple observations are in the same grid box then SE has blocks of
upper
triangular matrices that can be used to simulate the prediction error
dependence among observations in the same grid box. 
Explicitly if \code{obj} is the output object and there are \code{nObs}
observations then 

\preformatted{error <- obj$SE\%*\% rnorm( nObs)} 
will simulate a prediction error that includes the dependence. Note that
in
the case that there all observations are in separate grid boxes this code line is the same as 
\preformatted{error <- sqrt(obj$predictionVariance)*rnorm( nObs)}
It is always true that the prediction variance is given by 
\code{ diag( obj$SE\%*\% t( obj$SE))}.

The user is also referred to the testing scripts
\code{offGridWeights.test.R}  and 
\code{offGridWeights.testNEW.R}in \code{tests} where the 
Kriging predictions and standard errors are computed explicitly and 
tested against the sparse
matrix computation. This is helpful in defining exactly what
is being computed. 

Returned value below pertains to the offGridWeights function. 

\strong{addMarginsGridList} This is a supporting function that adds 
extra grid points for a gridList so that every irregular point has a
complete number of nearest neighbors.

\strong{findGridBox}This is a handy function that finds the indices of
the lower left grid box that contains the points in \code{s}. If the
point is not contained within the range of the grid and NA is returned 
fo the index. 
This function assumes that the grid points are equally spaced.}


\value{

\item{B}{A sparse matrix that is of dimension mXn with m the number of
locations (rows) in \code{s} and n being the total number of grid points.
\code{n = length(gridList$x)*length(gridList$y) }
}

\item{predictionVariance}{A vector of length as the rows of \code{s}
with the Kriging prediction variance based on the nearest neighbor
prediction and the specified covariance function. }
 \item{SE}{A sparse matrix that can be used to simulate dependence among
prediction errors for observations in the same grid box. 
(See explanation above.)}
 
}

\references{ 
Bailey, Maggie D., Soutir Bandyopadhyay, and Douglas Nychka. "Adapting conditional simulation using circulant embedding for irregularly spaced spatial data." Stat 11.1 (2022): e446.
}

\author{ Douglas Nychka and Maggie Bailey
}

\seealso{
\link{interp.surface} , \link{mKrigFastPredict}

}
\examples{

# an M by M  grid
M<- 400
xGrid<- seq( -1, 1, length.out=M)
gridList<- list( x= xGrid,
                 y= xGrid
                 )
NNSize<- 3 
 n<- 100
# sample n locations but avoid margins so that every point has a full set of nearest neighbors

set.seed(123)
s<- matrix(    runif(n*2, xGrid[(NNSize+1)],xGrid[(M-NNSize)]),
               n, 2 )

                  
obj<- offGridWeights( s, gridList, NNSize=3,
                   Covariance="Matern",
                   aRange = .1, sigma2= 1.0,
                   covArgs= list( smoothness=1.0)
                   )
		   
# make the predictions  by obj$B%*%c(y)
# where y is the matrix of values on the grid
 
# try it out on a simulated  Matern field  
CEobj<- circulantEmbeddingSetup( gridList,  
                  cov.args=list(
                  Covariance="Matern",
                   aRange = .1,
                    smoothness=1.0)
                    )
 set.seed( 333)                   
Z<- circulantEmbedding(CEobj)

#
# Note that grid values are "unrolled" as a vector
# for multiplication
# predOffGrid<- obj$B%*% c( Z)

predOffGrid<- obj$B\%*\% c( Z)

set.panel( 1,2)
zr<- range( c(Z))
image.plot(gridList$x, gridList$y, Z, zlim=zr)
bubblePlot( s[,1],s[,2], z= predOffGrid , size=.5,
highlight=FALSE, zlim=zr)
 set.panel()
 
}
\keyword{spatial}
