%#
%# fields  is a package for analysis of spatial data written for
%# the R software environment.
%# Copyright (C) 2021 Colorado School of Mines
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

\name{fields}
\alias{fields-package}
\alias{fields}
\title{
fields - tools for spatial data 
}
\description{
 \code{fields} is a collection of functions for curve and function
 fitting with an emphasis on spatial data and spatial statistics. It was developed over 20+ years to provide easy to use but sophisticated tools for analyzing
 spatial data, particularly that encountered in the environmental sciences. For the impatient users, jump to the examples below to see how easy this is use. 
 Please send bugs and questions to Doug Nychka, nychka@mines.edu. Positive
 comments are also welcome!

The
 major methods implemented include cubic and thin plate splines, 
 universal
 Kriging and Kriging for large data sets. A more modern framework  for Kriging is 
 spatial process estimation with covariance parameters determined by maximum likelihood and the uncertainty derived from assumptions of a Gaussian process. 
 Throughout we try to include reasonable defaults in functions that reflect
 our experience with analyzing spatial data. For example, the Matern covariance function is the default choice for the main spatial method. 
 
 A key feature of this package is any
 covariance function implemented in R code can be used for spatial prediction through the spatial functions.  
 Another important feature is that fields will take advantage of compactly
 supported covariance functions in a seamless way through 
the spam package. 
See \code{library( help=fields)} for a listing of all the
fields contents. We also  recommend the thoughtful vignette created by Ashton Weins, Mitchell Krock, and Emma Lilly (\code{fieldsVignette.pdf}) in the 
\href{https://github.com/dnychka/fieldsRPackage}{fields github repository}.

\code{fields} strives to have readable and tutorial code. Take a look at the 
source code for  \code{mKrig} to see how things work 
"under the hood". E.g. how a linear algebra computation is
overloaded to handle sparse matrices and how an output  object is built up sequentially throughout a computation.

The \code{fields} source code is liberally commented. Unfortunately on loading
this package,  R will
strip comments from the source for efficiency.   You can go to 
 \href{https://cran.r-project.org/package=fields}{CRAN fields page} to download
 the latest  "tarball" ( aka Package Source) and unzip to get code with comments.
 We also keep the most recent version of this package  at the 
\href{https://github.com/dnychka/fieldsRPackage}{fields github repository}.
and for commented source go to the 
the subdirectory \code{fields/src}. This may be a more recent version than what is posted to CRAN. 
}

\details{

\strong{Major methods} 
\itemize{ 

\item \code{\link{spatialProcess}} An easy to use method that fits a spatial process model
( e.g. Kriging) but also estimates the key spatial parameters:  nugget variance, sill variance and range by maximum likelihood. Default covariance model is a Matern covariance function. This function and related functions called by this are the core methods in fields and have much flexibility. 

\code{\link{spatialProcess}} allows one to supply a covariance
function that is
written in native R code.  See (\code{\link{stationary.cov}}) that
includes
several families of covariances  including the
Matern and several distance metrics including great circle distance
. 
\code{\link{sim.spatialProcess}} and \code{\link{simLocal.spatialProcess}} provide  "one liners"" for conditional simulation of the fitted surface. 
 
\item \code{\link{Tps}} Thin Plate spline
regression including GCV and REML estimates for the smoothing parameter. For moderate size data sets as a first look we use \code{Tps} all the time. See also \code{\link{fastTps}} for an 
approximate method to handle very large numbers of spatial locations 

\item \code{\link{sreg}} , \code{\link{splint}}   Fast 1-D cubic smoothing 
splines and interpolating  splines, a workhorse  algorithm for more EDA and more complicated methods.  

\item \code{\link{mKrig}} (micro Krig) 
    Efficient Universal Kriging and Gaussian process function,
    that can take advantage of sparse covariance
    
functions and is the core algorithm called by optimization
functions and for spatial predictio. 

\item \code{\link{QTps}} A easy to use extension of thin plate splines for quantile and robust surface fitting.

\item \code{\link{mKrigMLEGrid}} and \code{\link{mKrigMLEJoint}} for maximum likelihood estimates of covariance parameters. These functions also 
handle replicate fields, assumed to be independent realizations, at the same locations and can also take any covariate function function written in R following the \code{fields} format 
}  % end itemize


\strong{Other noteworthy functions}
\itemize{

\item \code{\link{vgram}} and \code{\link{vgram.matrix}} find variograms for spatial data (and
with temporal replications.

\item \code{\link{cover.design}} Generates space-filling designs where the distance 
function is expresed in R code.

\item Many
convenient functions for working with image data and rationally (well,
maybe reasonably) creating and  placing a color scale on  plots: 
\code{\link{as.image}}, \code{\link{imagePlot}}, \code{bubblePlot}, \code{
\link{drape.plot}}, \code{\link{quilt.plot}}
\code{\link{add.image}}, \code{\link{crop.image}}, \code{\link{half.image}},
\code{\link{average.image}}, 
\code{\link{designer.colors}}, \code{\link{color.scale}},
\code{\link{in.poly}}  
See also  \code{\link{grid.list}} for how fields works with grids and
\code{\link{US}}
and \code{\link{world}} for adding a map quickly. 
}


\strong{ Generic functions that support the methods} 

\code{plot} - diagnostic plots of fit \cr
\code{summary}- statistical summary of fit \cr
\code{print}- shorter version of summary \cr
\code{\link{surface}}- graphical display of fitted surface \cr
\code{predict}- evaluation fit at arbitrary points \cr
\code{\link{predictSE}}- prediction standard errors at arbitrary points. \cr
\code{\link{sim.rf}}- Simulate a random fields on a 2-d grid.

} % end details

\section{Getting Started}{

 Try some of the examples from help files for 
\code{\link{spatialProcess}} or \code{\link{Tps}}.

}

\section{\strong{Some fields datasets}}{
\itemize{
\item \code{\link{CO2}} Global satelite CO2 concentrations (simulated field)
\item \code{\link{COmonthlyMet}} Monthly mean temperatures and precip for Colorado
\item \code{\link{glacier}} An elevation dataset of a glacier also used by the applied math community to test interpolation methods. 
\item \code{\link{lennon}} Image of John Lennon
\item \code{\link{NorthAmericanRainfall}} 50+ year average and trend for summer rainfall at 
1700+ stations. 
\item \code{\link{ozone2}}  Daily max 8 hour ozone concentrations for the US midwest 
for summer 1987.
\item \code{\link{PRISMelevation}} Digital elevations for the 
   continental US at approximately 4km resolution
\item \code{\link{rat.diet}} Small paired study on rat food
intake over time.
\item \code{\link{RCMexample}} Regional climate model output
\item \code{\link{RMelevation}} Digital elevations for the Rocky
Mountain Empire

\item \code{\link{WorldBankCO2}} Demographic and carbon emission data 
 for 75 countries and for 1999. 
}
}


\section{\strong{DISCLAIMER:}}{
 The
authors can not guarantee the correctness of any function or program in
this package. 
}

\examples{
# some air quality data, daily surface ozone measurements for 
# the Midwest:

data(ozone2)
s<-ozone2$lon.lat
y<- ozone2$y[16,] # June 18, 1987 

# quick plot of spatial data with map
bubblePlot( s,y)
US( add=TRUE) # add US map

# fitting a thin plate spline surface (always a good place to 
# start). Here the  default smoothing (aka lambda) found by cross-validation
  fit0<- Tps(s,y)
# fits a GCV thin plate smoothing spline surface to ozone measurements.
# Hey, it does not get any easier than this!

  summary(fit0) #diagnostic summary of the fit 
  set.panel(2,2)
  plot(fit0) # four diagnostic plots of fit and residuals.

# quick plot of predicted surface
  set.panel()
  surface(fit0) # contour/image plot of the fitted surface
# see also predictSurface for more control over the evaluation grid
#
  US( add=TRUE, col="magenta", lwd=2) # US map overlaid
  title("Daily max 8 hour ozone in PPB,  June 18th, 1987")

####
  fit2<- spatialProcess( s,y)
# a "Kriging" model. The covariance defaults to a Matern 
# with smoothness 1.0.
# the nugget, sill and range parameters are found by maximum likelihood
# summary, plot, and surface also work for \code{fit2} !

  surface(fit2) # contour/image plot of the fitted surface
  US( add=TRUE, col="magenta", lwd=2) # US map overlaid
  title("Daily max 8 hour ozone in PPB,  June 18th, 1987")
\dontrun{
# And 20 approximate conditional draws of the spatial field on a grid
# with uncertainty in the 120PPB contour 
   look<- simLocal.spatialProcess(fit2, M=20)
for( k in 1:20){
contour( look$x, look$y, look$z[,,k], add=TRUE, level=c(120),
  col="white",  drawlabels=FALSE)
}

}  

}

\keyword{datasets}
