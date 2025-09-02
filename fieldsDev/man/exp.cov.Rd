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

\name{Covariance functions} 
\alias{Exp.cov}
\alias{Exp.simple.cov}
\alias{Rad.cov}
\alias{Rad.simple.cov}
\alias{stationary.cov}
\alias{stationary.taper.cov}
\alias{wendland.cov}
\alias{cubic.cov}
\alias{Tps.cov}
\alias{Paciorek.cov}
\title{
  Exponential  family, radial basis 
functions,cubic spline,  compactly  supported Wendland family, 
stationary covariances and non-stationary covariances. }
\description{
Given two sets of locations these functions compute the cross covariance matrix for
some covariance families. In addition these functions can take advantage 
of spareness, implement more efficient multiplcation of the 
cross covariance by a vector or matrix and also return a marginal 
variance to be consistent with calls by the Krig function.  

\code{stationary.cov} and \code{Exp.cov} have additional arguments for 
precomputed distance matrices and for calculating only the upper triangle 
and diagonal of the output covariance matrix to save time.  Also, they 
support using the \code{rdist} function with \code{compact=TRUE} or input 
distance matrices in compact form, where only the upper triangle of the 
distance matrix is used to save time.

Note: These functions have been been renamed from the previous fields functions
using 'Exp' in place of 'exp' to avoid conflict with the generic exponential 
function (\code{exp(...)})in R.
}
\usage{
Exp.cov(x1, x2 = NULL, aRange = 1, p = 1, distMat = NA, C =
                 NA, marginal = FALSE, onlyUpper = FALSE, theta = NULL,
                 ...)
                 
Exp.simple.cov(x1, x2, aRange =1, C=NA,marginal=FALSE, theta=NULL)

Rad.cov(x1, x2, p = 1, m=NA, with.log = TRUE, with.constant = TRUE, 
               C=NA,marginal=FALSE, derivative=0)

cubic.cov(x1, x2, aRange = 1, C=NA, marginal=FALSE, theta=NULL) 

Rad.simple.cov(x1, x2, p=1, with.log = TRUE, with.constant = TRUE, 
               C = NA, marginal=FALSE)

stationary.cov(x1, x2=NULL, Covariance = "Exponential", Distance = "rdist", 
  Dist.args = NULL, aRange = 1, V = NULL, C = NA, marginal = FALSE, 
  derivative = 0, distMat = NA, onlyUpper = FALSE, theta=NULL, ...)

stationary.taper.cov(x1, x2, Covariance="Exponential", 
           Taper="Wendland", 
           Dist.args=NULL, Taper.args=NULL, 
           aRange=1.0,V=NULL, C=NA, marginal=FALSE,
           spam.format=TRUE,verbose=FALSE, theta=NULL,...)
	   
Tps.cov(x1, x2 = NULL, cardinalX, m=2,
                   C = NA, aRange=NA,
                   marginal = FALSE
          ) 	   

wendland.cov(x1, x2, aRange = 1, V=NULL, k = 2, C = NA, 
             marginal =FALSE,Dist.args = list(method = "euclidean"), 
             spam.format = TRUE, derivative = 0, verbose=FALSE, theta=NULL)
Paciorek.cov(x1,
                           x2 = NULL,
                           Distance = "rdist",
                           Dist.args = NULL,
                           aRangeObj = 1,
                           rhoObj = NULL,
                           C = NA,
                           marginal = FALSE,
                           smoothness = .5)             
}
\arguments{
\item{x1}{ Matrix of first set of locations where each row gives the
coordinates of a particular point.}

\item{x2}{ Matrix of second set of locations where each row gives the
coordinatesof a particular point. If this is missing x1 is used. }

\item{aRange}{ Range (or scale) parameter. This should be a scalar (use
the V argument for other scaling options). Any distance calculated for
a covariance function is divided by aRange before the covariance function
is evaluated. For \code{Tps.cov} this argument is ignored.}

\item{aRangeObj}{ A fit object that defines the Range (or scale) parameter for arbitray locations using the generic \code{predict} function.}

\item{theta}{ Old version of the aRange parameter. If passed will be copied to aRange.}


\item{ cardinalX}{Locations added to the thin plate plate radial function to make it positive definite (See Details below).
}

\item{C}{ A vector or matrix with the same rows  as the number of rows of x2.
If specified the covariance matrix will be multiplied by this vector/matrix.}

\item{Covariance}{Character string that is the name of the covariance
shape function for the distance between locations. Choices in fields
are \code{Exponential}, \code{Matern}}

\item{derivative}{ If nonzero evaluates the partials of the
covariance function at locations x1. This must be used with  the "C" option
and is mainly called from within a predict function.  The partial
derivative is taken with respect to \code{x1}.  }

\item{Distance}{Character string that is the name of the distance
function to use. Choices in fields are \code{rdist},
\code{rdist.earth}}

\item{distMat}{
If the distance matrix between \code{x1} and \code{x2} has already been 
computed, it can be passed via this argument so it won't need to be 
recomputed.
}

\item{Dist.args}{ A list of optional arguments to pass to the Distance
function.}

\item{k}{The order of the Wendland covariance function.  See help on
 Wendland.}

\item{marginal}{If TRUE returns just the diagonal elements of the
covariance matrix using the \code{x1} locations. In this case this is
just 1.0. The marginal argument will trivial for this function is a
required argument and capability for all covariance functions used
with Krig.}



\item{m}{For the radial basis function p = 2m-d, with d being the
dimension of the locations, is the exponent applied to
the distance between locations. (m is a common way of parametrizing this
exponent.) Equivalently in Tps.cov the order of the spline. (See Details
section below).}

\item{onlyUpper}{
For internal use only, not meant to be set by the user.  Automatically 
set to \code{TRUE} by \code{mKrigMLEJoint} or \code{mKrigMLEGrid} if 
\code{lambda.profile} is set to \code{TRUE}, but set to \code{FALSE} 
for the final parameter fit so output is compatible with rest of 
\code{fields}.

If \code{TRUE}, only the upper triangle and diagonal of the covariance 
matrix is computed to save time (although if a non-compact distance 
matrix is used, the onlyUpper argument is set to FALSE). 
If \code{FALSE}, 
the entire covariance matrix is computed.  In general, it should 
only be set to \code{TRUE} for \code{mKrigMLEJoint} and
\code{mKrigMLEGrid}, 
and the default is set to \code{FALSE} so it is compatible with all of
\code{fields}.
}

\item{p}{ Exponent in the exponential covariance family.  p=1 gives an
exponential and p=2 gives a Gaussian.  Default is the exponential
form.  For the radial basis function this is the exponent applied to
the distance between locations. } 

\item{rhoObj}{
A fit object that defines a component of the marginal variance
(rho) parameter for
arbitray locations using the generic \code{predict} function. Note that in fields the complete marginal variance is \code{sigma2*rho} where \code{sigma2} can be estimated in \code{spatialProcess}.
 }
 
\item{smoothness}{For the Matern family the smoothnes of the process (aka "nu" in formulas). }

\item{spam.format}{If TRUE returns matrix in sparse matrix format
implemented in the spam package. If FALSE just returns a full
matrix. }

\item{Taper}{Character string that is the name of the taper function
to use. Choices in fields are listed in help(taper).}

\item{Taper.args}{ A list of optional arguments to pass to the Taper
function. \code{aRange} should always be the name for the range (or
scale) paremeter.}

\item{V}{ A matrix that describes the inverse linear transformation of
the coordinates before distances are found.  In R code this
transformation is: \code{x1 \%*\% t(solve(V))} Default is NULL, that
is the transformation is just dividing distance by the scalar value
\code{aRange}. See Details below.  If one has a vector of "aRange's"
that are the scaling for each coordinate then just express this as
\code{V = diag(aRange)} in the call to this function.}

\item{verbose}{If TRUE prints out some useful information for
debugging.}

\item{with.constant}{ If TRUE includes complicated constant for radial
 basis functions.  See the function \code{radbad.constant} for more
 details. The default is TRUE, include the constant. Without the usual
 constant the lambda used here will differ by a constant from spline
 estimators ( e.g.  cubic smoothing splines) that use the
 constant. Also a negative value for the constant may be necessary to
 make the radial basis positive definite as opposed to negative
 definite. }

\item{with.log}{If TRUE include a log term for even dimensions.  This
is needed to be a thin plate spline of integer order. }


\item{\dots}{ Any other arguments that will be passed to the
covariance function. e.g. \code{smoothness} for the Matern.}  }

\value{ If the argument C is NULL the cross covariance matrix is
returned.  In general if nrow(x1)=m and nrow(x2)=n then the returned
matrix will be mXn.  Moreover, if x1 is equal to x2 then this is the
covariance matrix for this set of locations.


If C is a vector of length n, then returned value is the
multiplication of the cross covariance matrix with this vector.
 
} \details{ For purposes of illustration, the function
\code{Exp.cov.simple} is provided in fields as a simple example and
implements the R code discussed below.  List this function out as a
way to see the standard set of arguments that fields uses to define a
covariance function.  It can also serve as a template for creating new
covariance functions for the \code{Krig} and \code{mKrig}
functions. Also see the higher level function \code{stationary.cov} to
mix and match different covariance shapes and distance functions.

A common scaling for stationary covariances: If \code{x1} and
 \code{x2} are matrices where \code{nrow(x1)=m} and \code{nrow(x2)=n}
 then this function will return a mXn matrix where the (i,j) element
 is the covariance between the locations \code{x1[i,]} and
 \code{x2[j,]}. The exponential covariance function is computed as
 exp( -(D.ij)) where D.ij is a distance between \code{x1[i,]} and
 \code{x2[j,]} but having first been scaled by aRange. Specifically if
 \code{aRange} is a matrix to represent a linear transformation of the
 coordinates, then let \code{u= x1\%*\% t(solve( aRange))} and \code{v=
 x2\%*\% t(solve(aRange))}.  Form the mXn distance matrix with
 elements:

\code{D[i,j] = sqrt( sum( ( u[i,] - v[j,])**2 ) )}.

and the cross covariance matrix is found by \code{exp(-D)}.  The
tapered form (ignoring scaling parameters) is a matrix with i,j entry
\code{exp(-D[i,j])*T(D[i,j])}. With T being a positive definite
tapering function that is also assumed to be zero beyond 1.

Note that if aRange is a scalar then this defines an isotropic
covariance function and the functional form is essentially
\code{exp(-D/aRange)}.

Implementation: The function \code{r.dist} is a useful FIELDS function
that finds the cross Euclidean distance matrix (D defined above) for
two sets of locations. Thus in compact R code we have

  exp(-rdist(u, v))

Note that this function must also support two other kinds of calls:

If marginal is TRUE then just the diagonal elements are returned (in R
code \code{diag( exp(-rdist(u,u)) )}).

If C is passed then the returned value is \code{ exp(-rdist(u, v))
\%*\% C}.


Some details on particular covariance functions:

\describe{


\item{Stationary covariance \code{stationary.cov}:}{Here the
computation is to apply the function Covariance to the distances found
by the Distance function.  For example

\code{Exp.cov(x1,x2, aRange=MyTheta)}

and

\code{stationary.cov( x1,x2, aRange=MyTheta, Distance= "rdist",
Covariance="Exponential")}

are the same. This also the same as

\code{stationary.cov( x1,x2, aRange=MyTheta, Distance= "rdist",
Covariance="Matern",smoothness=.5)}.  }

\item{Radial basis functions (\code{Rad.cov}:}{The
functional form is Constant* rdist(u, v)**p for odd dimensions and
Constant* rdist(u,v)**p * log( rdist(u,v) ) For an m th order thin plate
spline in d dimensions p= 2*m-d and must be positive. The constant,
depending on m and d, is coded in the fields function
\code{radbas.constant}. This form is only a generalized covariance
function -- it is only positive definite when restricted to linear
subspace.  See \code{Rad.simple.cov} for a coding of the radial basis
functions in R code.}

\item{Tps.cov}{This covariance can be used in a standard
"Kriging" computation to give a thin-plate spline (TPS). This is
useful because one can use the high
level function \code{spatialProcess} and supporting functions for
the returned
object, including conditional simulation.
The standard computation for a TPS uses the radial basis functions
as given
in \code{Rad.cov} and uses a QR decomposition based a polynoimial
matrix to
reduce the dimension of the radial basis function and yield a
positive definite
matrix. This reduced matrix is then used in the regular
compuations to find
the spatial estimate. The function \code{Krig} and specifically 
\code{Tps}
implements this algoritm.  The interested reader should look at
\code{Krig.engine.default} and specifically at the \code{TMatrix}
polynomial matrix and resulting reduced positive definite matrix
\code{tempM}. The difficulty with this approach is that is not
amenable to taking advantage of sparsity in the covariance matrix.

An alternative that is suggested by Grace Wahba in \emph{Spline
models for
observational
data} is to augment the radial basis functions with a low rank set
of
functions based on a low order polynomial evaluated at a set of
points. The set of locatios for this modifications are called
\emph{cardinal point}s due the to
property mentioned below.  This is
implemented in \code{Tps.cov} leading to a full rank  
(non-stationary!)
covariance function. If the fixed part of the spatial model also
includes
this same order polynomial then the final result gives a TPS and
is invariant
to the choice of cardinal points. To streamline using this
covariance when it isspecified in the \code{spatialProcess}
function the cardinal points will choosen automaitcally based on
the observation locations and  the spline order, \code{m} using a
space filling design.
A simple example  with fixed smoothing parameter, \code{lambda <- .1} may help


\preformatted{
data( "ozone2")
s<- ozone2$lon.lat
y<- ozone2$y[16,]
 
fitTps1<- spatialProcess( s,y, cov.function= "Tps.cov", lambda=.1)               
 }
 
 and
 compare the results to the standard algorithm.

\preformatted{
fitTps2<- Tps( s,y, scale.type ="unscaled", lambda=.1)

stats( abs(fitTps2$fitted.values - fitTps1$fitted.values)) 
}
 
 Here the default choice for the order is 2 and in two dimensions implies a
 linear polynomial. The arguments filled in by default are shown below
 
\preformatted{
 fitTps1$args
$cardinalX
        [,1]   [,2]
[1,] -85.289 40.981
[2,] -90.160 38.330
[3,] -91.229 43.812
attr(,"scaled:scale")
[1] 1 1
attr(,"scaled:center")
[1] 0 0

$aRange
[1] NA

fitTps1$mKrig.args
[[1]]
NULL

$m
[1] 2

$collapseFixedEffect
[1] TRUE

$find.trA
[1] TRUE

}

\code{ cardinalX} are the cardinal points chosen using a space
filling
criterion. These are attached to the covariance arguments list so
they are used
in the subsequent computations for this fit (such as  predict,
predictSE, sim.spatialProcess).

In general, if \code{d} is the dimension of the locations and
\code{m} the order of the spline one must have \code{2*m-d >0}.
The polynomial will have \code{ choose( m+d-1,d)} terms and so
this many cardinal points need to be specified. As mentioned above
these are chosen in a reasonable way if \code{spatialProcess} is
used.  
}

\item{Stationary tapered covariance \code{stationary.taper.cov}
:}{The
resulting cross covariance is the direct or Shure product of the
tapering function and the covariance. In R code given location
matrices, \code{x1} and \code{x2} and using Euclidean distance.
 
\code{Covariance(rdist( x1, x2)/aRange)*Taper( rdist( x1,
x2)/Taper.args$aRange)}

By convention, the \code{Taper} function is assumed to be identically
zero outside the interval [0,1]. Some efficiency is introduced within
the function to search for pairs of locations that are nonzero with
respect to the Taper. This is done by the SPAM function
\code{nearest.dist}.  This search may find more nonzero pairs than
dimensioned internally and SPAM will try to increase the space. One
can also reset the SPAM options to avoid these warnings.  For
spam.format TRUE the multiplication with the \code{C} argument is done
with the spam sparse multiplication routines through the "overloading"
of the \code{\%*\%} operator.  }


\item{Nonstationary covariance function, Paciorek.cov}{
This implements the nonstationary model developed by Chris Paciorek and
Mark Schervish that allows for a varying range parameter over space and
also a varying marginal variance. This is still experimental and spatialProcess has not been generalized to fit the parameter surfaces.
It can, however, be used to evaluate the model at fixed parameter surfaces. See the last example below. 

This covariance works by specifying a object aRangeObj such that the
generic call  \code{predict(aRangeObj, loc)} will evaluate the aRange
function at the locations \code{loc}. This object can be as simple a
fit to local estimated aRange parameters using a fields function such
as Tps or spatialProcess. More specific applications one can create a
special predict function. For example suppose log aRange follows a
linear model in the spatial coordinates and these coefficients are the
vector \code{Afit}. Define a class and a predict function. 

\preformatted{ 

aRangeObj<- list(coef=Afit)
class(aRangeObj)<- "myclass"

predict.myclass<- function( aRangeObj, x){
aRange<- exp(cbind( 1,x) \%*\% aRangeObj$coef)
return( aRange)
}
}

Now use \code{spatialProcess} with this object and make sure
\code{predict.myclass} is also loaded.

}
A similar strategy will also work for a varying marginal variance by
creating \code{sigmaObj} and if needed a companion predict method. 

}

About the FORTRAN: The actual function \code{Exp.cov} and 
\code{Rad.cov} call FORTRAN to 
make the evaluation more efficient this is especially important when the 
C argument is supplied. So unfortunately the actual production code in 
Exp.cov is not as crisp as the R code sketched above. See  
\code{Rad.simple.cov} for a R coding of the radial basis functions.  

}

\seealso{
 Krig, rdist, rdist.earth, gauss.cov, Exp.image.cov, Exponential, Matern, 
Wendland.cov, mKrig} 

\examples{
# exponential covariance matrix ( marginal variance =1) for the ozone
#locations 
out<- Exp.cov( ChicagoO3$x, aRange=100)

# out is a 20X20 matrix

out2<- Exp.cov( ChicagoO3$x[6:20,],ChicagoO3$x[1:2,], aRange=100)
# out2 is 15X2 matrix 

# Kriging fit where the nugget variance is found by GCV 
# Matern covariance shape with range of 100.
# 

fit<- Krig( ChicagoO3$x, ChicagoO3$y, Covariance="Matern", aRange=100,smoothness=2)

data( ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]
# Omit the NAs
good<- !is.na( y)
x<- x[good,]
y<- y[good]


# example of calling the taper version directly 
# Note that default covariance is exponential and default taper is 
# Wendland (k=2).

stationary.taper.cov( x[1:3,],x[1:10,] , aRange=1.5, Taper.args= list(k=2,aRange=2.0,
                       dimension=2) )-> temp
# temp is now a tapered 3X10 cross covariance matrix in sparse format. 

 is.spam( temp)  # evaluates to TRUE

# should be identical to
# the direct matrix product

 temp2<- Exp.cov( x[1:3,],x[1:10,], aRange=1.5) * Wendland(rdist(x[1:3,],x[1:10,]), 
                      aRange= 2.0, k=2, dimension=2)
 test.for.zero(  as.matrix(temp), temp2)

# Testing that the "V" option works as advertized ...
x1<- x[1:20,]
x2<- x[1:10,]

V<- matrix( c(2,1,0,4), 2,2)
Vi<- solve( V)

u1<- t(Vi\%*\% t(x1))
u2<- t(Vi\%*\% t(x2))

look<- exp(-1*rdist(u1,u2))
look2<- stationary.cov( x1,x2, V= V)
test.for.zero( look, look2)


# Here is an example of how the cross covariance multiply works
# and lots of options on the arguments


 Ctest<- rnorm(10)
 
 temp<- stationary.cov( x,x[1:10,], C= Ctest, 
        Covariance= "Wendland", 
            k=2, dimension=2, aRange=1.5 )

# do multiply explicitly

 temp2<- stationary.cov( x,x[1:10,],
        Covariance= "Wendland",
            k=2, dimension=2, aRange=1.5 )\%*\% Ctest

 test.for.zero( temp, temp2)


# use the tapered stationary version 
# cov.args is part of the argument list passed to stationary.taper.cov
# within Krig. 
# This example needs the spam package.
# 

\dontrun{

Krig(x,y, cov.function = "stationary.taper.cov", aRange=1.5,
      cov.args= list(Taper.args= list(k=2, dimension=2,aRange=2.0) )
           ) -> out2 
# NOTE: Wendland is the default taper here. 
}

# BTW  this is very similar to 
\dontrun{
 Krig(x,y, aRange= 1.5)-> out
}

##################################################
#### nonstationary covariance Paciorek.cov
##################################################
\dontrun{
M<- 20
gridList<- list(x=seq( 0,1,length.out=M),
                y=seq( 0,1,length.out=M))
sGrid<- make.surface.grid(gridList )
# An aRange surface
aRangeObj<- list(coef= c( 1,4,0))
class(aRangeObj)<- "myclass"

predict.myclass<- function( aRangeObj, x){
aRange<- exp(cbind( 1,x) \%*\% aRangeObj$coef)
return( aRange)
}

covMatrix<- Paciorek.cov( sGrid, sGrid, aRangeObj=aRangeObj)
# examine correlation surface between selected locations and the  full grid. 
set.panel( 2,2)
{
imagePlot( as.surface(sGrid, covMatrix[,10]))
imagePlot( as.surface(sGrid, covMatrix[,205]))
imagePlot( as.surface(sGrid, covMatrix[,305]))
imagePlot( as.surface(sGrid, covMatrix[,390]))

}

# simulation of the field
set.seed(222)
n<- nrow( sGrid)
f<- t(chol(covMatrix)) \%*\% rnorm(M^2)
set.panel()
imagePlot( as.surface(sGrid,f))
y<- f + .05*rnorm(n)
fitP<- spatialProcess( sGrid, y, cov.function="Paciorek.cov",
           cov.args= list(aRangeObj = aRangeObj ) )
# check estimated tau and sigma 
fitP$summary

# fitted surface
surface( fitP)

}

}
\keyword{spatial}
% docclass is function

