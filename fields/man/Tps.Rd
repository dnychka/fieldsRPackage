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
\name{Tps}
\alias{Tps}
\alias{fastTps}


\title{
  Thin plate spline regression 
}
\description{
Fits a thin plate spline surface to irregularly spaced data. The 
smoothing parameter is chosen  as a default by generalized cross-validation. The assumed 
model is additive  Y = f(X) +e  where f(X) is a d dimensional surface.

}
\usage{
Tps(x, Y, m = NULL, p = NULL, scale.type = "range", lon.lat = FALSE,
         miles = TRUE, method = "GCV", GCV = TRUE, ...)

fastTps(x, Y, m = NULL, p = NULL, aRange, lon.lat = FALSE,
                 find.trA = FALSE, REML = FALSE,theta=NULL,  ...)
		 
        
}
\arguments{
%To be helpful, a more complete list of arguments are described that are the 
%same as those for the Krig function.  
\item{x}{
Matrix of independent variables. Each row is a location or a set of 
independent covariates. 
}
\item{Y}{
Vector of dependent variables. 
}
\item{m}{
A polynomial function of degree (m-1) will be  
included in the model as the drift (or spatial trend) component. 
Default is the value such that 2m-d is greater than zero where d is the 
dimension of x. 
}

\item{p}{
Polynomial power for Wendland radial basis functions. Default is 2m-d 
where d is the dimension of x.  
}

\item{scale.type}{
The independent variables and knots are scaled to the specified 
scale.type. 
By default the scale type is "range", whereby 
the locations are transformed  
 to the interval (0,1) by forming (x-min(x))/range(x) for each x. 
Scale type of "user" allows specification of an x.center and
x.scale by the 
user. The default for "user" is mean 0 and standard deviation 1. Scale 
type of "unscaled" does not scale the data.  
}
\item{aRange}{The tapering range that is passed to the Wendland compactly 
supported covariance. The covariance (i.e. the radial basis function) is 
zero beyond range aRange. The larger aRange the closer this model will
approximate the standard thin plate spline.}

\item{lon.lat}{If TRUE locations are interpreted as lognitude and
latitude and great circle distance is used to find distances among
locations. The aRange scale parameter for \code{fast.Tps} (setting the
compact support of the Wendland function) in this case is in units of
miles (see example and caution below). }
 
 \item{method}{ Determines what "smoothing" parameter should be
     used. The default is to estimate standard GCV Other choices are:
     GCV.model, GCV.one, RMSE, pure error and REML. The differences
     are explained in the Krig help file.}
     
 \item{GCV}{If TRUE the decompositions are done to efficiently evaluate
 the estimate, GCV function and likelihood at multiple values
 of lambda. }    
 
\item{miles}{If TRUE great circle distances are in miles if FALSE
distances are in kilometers}


\item{find.trA}{If TRUE will estimate the effective degrees of freedom
using a simple Monte Carlo method (random trace). This will add to the
computational burden by approximately \code{NtrA} solutions of the
linear system but
the cholesky decomposition is reused.}

\item{REML}{If TRUE find the MLE for lambda using restricted maximum likelihood instead of the full version. }

\item{theta}{ Same as aRange.}
  
\item{\dots}{For \code{Tps} any argument that is valid for the
\code{Krig} function. Some of the main ones are listed below.

For \code{fastTps} any argument that is suitable for the \code{mKrig}
function see help on mKrig for these choices. The most common would be
\code{lambda} fixing the value of this parameter (tau^2/sigma^2), \code{Z} linear covariates or 
\code{ mKrig.args= list( m=1)} setting the regression model to be just a constant function.

Arguments for Tps:
\describe{
  \item{lambda}{ Smoothing parameter that is the ratio of the error
    variance (tau**2) to the scale parameter of the covariance
    function. If omitted this is estimated by GCV. }

  \item{Z}{Linear covariates to be included in fixed part of the model
    that are distinct from the default low order polynomial in
    \code{x}}

  \item{df}{ The effective number of parameters for the fitted
    surface. Conversely, N- df, where N is the total number of
    observations is the degrees of freedom associated with the
    residuals.  This is an alternative to specifying lambda and much
    more interpretable.}

  \item{cost}{ Cost value used in GCV criterion. Corresponds to a
    penalty for increased number of parameters. The default is 1.0 and
    corresponds to the usual GCV.}

  \item{weights}{ Weights are proportional to the reciprocal variance
    of the measurement error. The default is no weighting i.e. vector
    of unit weights. }

  \item{nstep.cv}{ Number of grid points for minimum GCV search. }

  \item{x.center}{ Centering values are subtracted from each column of
     the x matrix.  Must have scale.type="user".}

  \item{x.scale}{ Scale values that divided into each column after
    centering.  Must have scale.type="user".}

  \item{sigma}{Scale factor for covariance. }

  \item{tau2}{ Variance of errors or if weights are not equal to 1
     the variance is tau**2/weight.}

 

   \item{verbose}{ If true will print out all kinds of intermediate
      stuff.  }

   \item{mean.obj}{ Object to predict the mean of the spatial
      process. }

   \item{sd.obj}{ Object to predict the marginal standard deviation of
      the spatial process. }

    \item{null.function}{An R function that creates the matrices for
       the null space model.  The default is fields.mkpoly, an R
       function that creates a polynomial regression matrix with all
       terms up to degree m-1. (See Details) }

     \item{offset}{ The offset to be used in the GCV
        criterion. Default is 0. This would be used when Krig/Tps is
        part of a backfitting algorithm and the offset has to be
        included to reflect other model degrees of freedom.  }
}
}


}

\value{
A list of class Krig. This includes the 
fitted values, the predicted surface evaluated at the 
observation locations, and the residuals. The results of the grid 
search minimizing the generalized cross validation function are
returned in gcv.grid. Note that the GCV/REML optimization is 
done even if lambda or df is given. 
Please see the documentation on Krig for details of the returned 
arguments.  
}
\details{


\strong{Overview}

This is the classic nonparametric curve/surface estimate pioneered in statistics by Grace Wahba. The computational algorithm is based around a QR decomposition followed by an eigen decomposition on a reduced matrix derived from the spline radial basis functions. This insures a stable computation -- basically bombproof -- but is not the fastest. See the function \link{mKrig} and examples for fitting a thin plate spline using \link{spatialProcess} that uses a different linear algebra approach. This is implemented in the "fast" version albeit approximate version to exploit sparse matrices. 

This function also works for just a single dimension and reproduces the well known cubic smoothing spline for \code{ m ==2} Finally we note that a thin plate spline is a limiting case of a Gaussian process estimate as the range parameter in the Matern family increases to infinity. 
(Kriging). 

A "fast" version of this function uses a compactly supported Wendland covariance and sparse linear algebra for handling larger datta sets. Although a good approximation to Tps for sufficiently large aRange its actual form is very different from the textbook thin-plate definition. The user will see that \code{fastTps} is largely a wrapper for a call to \code{spatialProcess} with the Wendland covariance function. 

\strong{Background on the thin plate spline function, Tps}
A thin plate spline is the result of minimizing the residual sum of 
squares subject to a constraint that the function have a certain 
level of smoothness (or roughness penalty). Roughness is 
quantified by the integral of squared m-th order derivatives. For one 
dimension and m=2 the roughness penalty is the integrated square of 
the second derivative of the function. For two dimensions the 
roughness penalty is the integral of  

        (Dxx(f))**22 + 2(Dxy(f))**2 + (Dyy(f))**22 

(where Duv denotes the second partial derivative with respect to u 
and v.) Besides controlling the order of the derivatives, the value of 
m also determines the base polynomial that is fit to the data. 
The degree of this polynomial is (m-1). 

The smoothing parameter controls the amount that the data is 
smoothed. In the usual form this is denoted by lambda, the Lagrange 
multiplier of the minimization problem. Although this is an awkward 
scale, lambda = 0 corresponds to no smoothness constraints and the data 
is interpolated.  lambda=infinity corresponds to just fitting the 
polynomial base model by ordinary least squares.  

This estimator is implemented by passing the right generalized covariance
function based on radial basis functions to the more general function
Krig.  One advantage of this implementation is that once a Tps/Krig object
is created the estimator can be found rapidly for other data and smoothing
parameters provided the locations remain unchanged. This makes simulation
within R efficient (see example below). Tps does not currently support the
knots argument where one can use a reduced set of basis functions. This is
mainly to simplify the code and a good alternative using knots would be to use a
valid covariance from the Matern family and a large range parameter. 

\strong{Using a great circle distance function}
 The option to use great circle distance
 to define the radial basis functions (\code{lon.lat=TRUE}) is very useful
 for small geographic domains where the spherical geometry is well approximated by a plane. However, for large domains the spherical distortion be large enough that the radial basis functions no longer define a positive definite system and Tps will report a numerical error. An alternative is to switch to a three
dimensional thin plate spline with the locations being the direction cosines. This will 
give approximate great circle distances for locations that are close and also the numerical methods will always have a positive definite matrices. There are other radial basis functions that are specific to a spherical geometry but these are not implemented in fields. 

Here is an example using this idea for \code{RMprecip} and also some 
examples of building grids and evaluating the Tps results on them:
\preformatted{
# a useful function:
  dircos<- function(x1){
             coslat1 <- cos((x1[, 2] * pi)/180)
             sinlat1 <- sin((x1[, 2] * pi)/180)
             coslon1 <- cos((x1[, 1] * pi)/180)
             sinlon1 <- sin((x1[, 1] * pi)/180)
             cbind(coslon1*coslat1, sinlon1*coslat1, sinlat1)}
# fit in 3-d to direction cosines
  out<- Tps(dircos(RMprecip$x),RMprecip$y)
  xg<-make.surface.grid(fields.x.to.grid(RMprecip$x))
  fhat<- predict( out, dircos(xg))
# coerce to image format from prediction vector and grid points.
  out.p<- as.surface( xg, fhat)
  surface( out.p)
# compare to the automatic
  out0<- Tps(RMprecip$x,RMprecip$y, lon.lat=TRUE)
  surface(out0)
}

The function \code{fastTps} is really a convenient wrapper function that
calls \code{spatialProcess} with a suitable Wendland covariance 
function. This means one can use all the additional functions for
prediction and simulation built for the \code{spatialProcess} and
\code{mKrig} objects.  Some care needs to exercised in specifying
the support \code{aRange} -- a Goldilocks problem -- where aRange
is large enough so that every location has a reasonable number
(say 10 or more ) of neighboring locations that have non-zero
covariances but also the number of neighbors is not so large that
the sparsity of the covariance matrix is compromised. 
To figure out the neighborhood pattern are the spatial locations one can use the function \code{nearest.dist}, sparse matrix format,  and \code{table} function. 
\preformatted{
set.seed(222)
s<- cbind( runif(1e4),runif(1e4))
look<- nearest.dist(s,s, delta = .03)
look<- spam2spind(look)

stats( table( look$ind[,1]))
}
Here one has asummary of the number of nearest neighbors within a 
distance of .03 for these (randomly generated) locations.
I see a minimum of 7 for at least one location so aRange should be larger than .03. Trial and error with different deltas can lead to a better choice.  

 Note that unlike Tps the locations are not 
scaled to unit range and this can cause havoc in smoothing problems with
variables in very different units. 
So rescaling the locations \code{ x<- scale(x)} 
is a good idea for putting the variables on a common scale
for smoothing.  A conservative rule of thumb is to make \code{aRange}
large enough so that about 50 nearest neighbors are within this distance
for every observation location. 

This function does have the potential to approximate estimates of Tps 
for very large spatial data sets. See \code{wendland.cov} and help on 
the SPAM package, and the friendly spind format for more background.
Also, the function \code{predictSurface.fastTps} has been made more
efficient for the 
case of k=2 and m=2 -- a common choice for parameters.

 }

\section{References}{
See "Nonparametric Regression and Generalized Linear Models"  

by Green and Silverman. 
See "Additive Models" by Hastie and Tibshirani. 

See Wahba, Grace. "Spline models for observational data." Society for industrial and applied mathematics, 1990.

}
\seealso{
\code{\link{Krig}},
\code{\link{mKrig}},
\code{\link{spatialProcess}},
\code{\link{Tps.cov}}
\code{\link{sim.spatialProcess}},
\code{\link{summary.Krig}},
\code{\link{predict.Krig}},
\code{\link{predictSE.Krig}},
\code{\link{predictSurface}},
\code{\link{predictSurface.fastTps}},
\code{\link{plot.Krig}},
\code{\link{surface.Krig}}, 
\code{\link{sreg}}

}
\examples{
#2-d example 
data(ozone2)

x<- ozone2$lon.lat
y<- ozone2$y[16,]

fit<- Tps(x,y)  # fits a surface to ozone measurements. 

set.panel(2,2)
plot(fit) # four diagnostic plots of  fit and residuals. 
set.panel()

# summary of fit and estiamtes of lambda the smoothing parameter
summary(fit)

surface( fit) # Quick image/contour plot of GCV surface.

# NOTE: the predict function is quite flexible:

     look<- predict( fit, lambda=2.0)
#  evaluates the estimate at lambda =2.0  _not_ the GCV estimate
#  it does so very efficiently from the Krig fit object.

     look<- predict( fit, df=7.5)
#  evaluates the estimate at the lambda values such that 
#  the effective degrees of freedom is 7.5
 

# compare this to fitting a thin plate spline with 
# lambda chosen so that there are 7.5 effective 
# degrees of freedom in estimate
# Note that the GCV function is still computed and minimized
# but the lambda values used correpsonds to 7.5 df.

fit1<- Tps(x, y,df=7.5)

set.panel(2,2)
plot(fit1) # four diagnostic plots of  fit and residuals.
          # GCV function (lower left) has vertical line at 7.5 df.
set.panel()

# The basic matrix decompositions are the same for 
# both fit and fit1 objects. 

# predict( fit1) is the same as predict( fit, df=7.5)
# predict( fit1, lambda= fit$lambda) is the same as predict(fit) 


# predict onto a grid that matches the ranges of the data.  

out.p<-predictSurface( fit)
imagePlot( out.p) 

# the surface function (e.g. surface( fit))  essentially combines 
# the two steps above

# predict at different effective 
# number of parameters 
out.p<-predictSurface( fit, df=10)

\dontrun{
# predicting on a grid along with a covariate
  data( COmonthlyMet)	
# predicting average daily minimum temps for spring in Colorado
# NOTE to create an  4km  elevation grid: 
# data(PRISMelevation); CO.elev1 <- crop.image(PRISMelevation, CO.loc )
# then use same grid for the predictions: CO.Grid1<- CO.elev1[c("x","y")]
  obj<- Tps( CO.loc, CO.tmin.MAM.climate, Z= CO.elev)
  out.p<-predictSurface( obj,
              CO.Grid, ZGrid= CO.elevGrid)
  imagePlot( out.p)        
  US(add=TRUE, col="grey")
  contour( CO.elevGrid, add=TRUE, levels=c(2000), col="black")
}
\dontrun{
#A 1-d example  with confidence intervals
  out<-Tps( rat.diet$t, rat.diet$trt) # lambda found by GCV 
  out
  plot( out$x, out$y)
  xgrid<- seq(  min( out$x), max( out$x),,100)
  fhat<- predict( out,xgrid)
  lines( xgrid, fhat,)
  SE<- predictSE( out, xgrid)
  lines( xgrid,fhat + 1.96* SE, col="red", lty=2)
  lines(xgrid, fhat - 1.96*SE, col="red", lty=2)

# 
# compare to the ( much faster) B spline algorithm 
#  sreg(rat.diet$t, rat.diet$trt) 

# Here is a 1-d example with 95 percent  CIs  where sreg would not 
# work:
#  sreg would give the right estimate here but not the right CI's
  x<- seq( 0,1,,8)
  y<- sin(3*x)
  out<-Tps( x, y) # lambda found by GCV 
  plot( out$x, out$y)
  xgrid<- seq(  min( out$x), max( out$x),,100)
  fhat<- predict( out,xgrid)
  lines( xgrid, fhat, lwd=2)
  SE<- predictSE( out, xgrid)
  lines( xgrid,fhat + 1.96* SE, col="red", lty=2)
  lines(xgrid, fhat - 1.96*SE, col="red", lty=2)
}

# More involved example adding a covariate to the fixed part of model
\dontrun{
set.panel( 1,3)
# without elevation covariate
  out0<-Tps( RMprecip$x,RMprecip$y)
  surface( out0)
  US( add=TRUE, col="grey")

# with elevation covariate
  out<- Tps( RMprecip$x,RMprecip$y, Z=RMprecip$elev)
  
# NOTE: out$d[4] is the estimated elevation coefficient
# it is easy to get the smooth surface separate from the elevation.
  out.p<-predictSurface( out, drop.Z=TRUE)
  surface( out.p)
  US( add=TRUE, col="grey")
  
# and if the estimate is of high resolution and you get by with 
# a simple discretizing -- does not work in this case!
  quilt.plot( out$x, out$fitted.values)
  
#
# the exact way to do this is evaluate the estimate
# on a grid where you also have elevations 
# An elevation DEM from the PRISM climate data product (4km resolution)
  data(RMelevation)  
  grid.list<- list( x=RMelevation$x, y= RMelevation$y)
  fit.full<- predictSurface( out, grid.list, ZGrid= RMelevation)
  
# this is the linear fixed part of the second spatial model:
# lon,lat and elevation
  fit.fixed<- predictSurface( out, grid.list, just.fixed=TRUE,
                   ZGrid= RMelevation)
                   
# This is the smooth part but also with the linear lon lat terms. 
  fit.smooth<-predictSurface( out, grid.list, drop.Z=TRUE)
  
#
  set.panel( 3,1)
  
  fit0<- predictSurface( out0, grid.list)
  image.plot( fit0)
  title(" first spatial model (w/o elevation)")
  image.plot( fit.fixed)
  title(" fixed part of second model (lon,lat,elev linear model)")
  US( add=TRUE)
  image.plot( fit.full)
  title("full prediction second model")
  set.panel()
}
### 
### fast Tps
# m=2   p= 2m-d= 2
#
# Note: aRange = 3 degrees is a very generous taper range. 
# Use some trial aRange value with rdist.nearest to determine a
# a useful taper. Some empirical studies suggest that in the 
# interpolation case in 2 d the taper should be large enough to 
# about 20 non zero nearest neighbors for every location.

  out2<- fastTps( RMprecip$x,RMprecip$y,m=2, aRange=4, 
                      profileLambda=FALSE)

# note that fastTps produces a object of classes  spatialProcess and mKrig
# so one can use all the 
# the overloaded functions that are defined for these classes.
# predict, predictSE, plot, sim.spatialProcess 
# summary of what happened note estimate of effective degrees of 
# freedom
# profiling on lambda has been turned off to make this run quickly
# but it is suggested that one examines the the profile likelihood over lambda
  
  print( out2)

\dontrun{
set.panel( 1,2)
surface( out2)

#
# now use great circle distance for this smooth 
# Here "aRange" for the taper support is the great circle distance in degrees latitude.
# Typically for data analysis it more convenient to think in degrees. A degree of
# latitude is about 68 miles (111 km).
#
out3<- fastTps( RMprecip$x,RMprecip$y,m=2, lon.lat=TRUE, aRange= 4*68 )
print( out3)  # note the effective degrees of freedom is different.
surface(out3)

set.panel()
}

\dontrun{
#
# simulation reusing Tps/Krig object
#
fit<- Tps( rat.diet$t, rat.diet$trt)
true<- fit$fitted.values
N<-  length( fit$y)
temp<- matrix(  NA, ncol=50, nrow=N)
tau<- fit$tauHat.GCV
for (  k in 1:50){
ysim<- true + tau* rnorm(N) 
temp[,k]<- predict(fit, y= ysim)
}
matplot( fit$x, temp, type="l")

}
# 
#4-d example 
fit<- Tps(BD[,1:4],BD$lnya,scale.type="range") 

# plots fitted surface and contours 
# default is to hold 3rd and 4th fixed at median values 

surface(fit)   



}
\keyword{smooth}
% docclass is function
% Converted by Sd2Rd version 1.21.
