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

\name{sim.spatialProcess}
\alias{sim.Krig}
\alias{simSpatialData}
\alias{sim.spatialProcess}
\alias{simLocal.spatialProcess}
\alias{checkPredictGrid}
\alias{makePredictionGridList}
\alias{makeSimulationGrid}

\title{Unconditional and conditional simulation of a spatial process}
\description{
Generates exact (or approximate) random draws from the unconditional or
conditional 
distribution of a spatial process given specific observations. 
Draws from the conditional distribution, known as conditional simulation in geostatistics,  is a 
useful way to characterize the uncertainty in the predicted process from 
data.  Note that exact simulation is limted by the number of locations but there are approximate strategies to handle simulation for large grids of locations.  
 }
\usage{

simSpatialData(object,  M = 1, verbose = FALSE)
sim.spatialProcess(object, xp,  M = 1, verbose = FALSE, ...)
sim.Krig(object, xp, M = 1, verbose = FALSE, ...)

simLocal.spatialProcess(mKrigObject, predictionGridList = NULL,
                 simulationGridList = NULL, gridRefinement = 1, np = 2,
                 M = 1, nx = 80, ny = 80, verbose = FALSE, delta =
                 NULL, giveWarnings = TRUE, fast = FALSE, NNSize = 5,
                 ...)
             
checkPredictGrid(predictionGridList) 

makePredictionGridList(mKrigObject, nx, ny, np)

makeSimulationGrid(predictionGridList, gridRefinement)   
} 
 
\arguments{
 
 \item{delta}{If the covariance has compact support the simulation method can 
take advantage of this. This is the amount of buffer added for the simulation domain in the circulant embedding method. 
A minimum size would be \code{aRange} for the Wendland but a multiple of this maybe needed to obtain a positive definite 
circulant covariance function.  }

 
  \item{fast}{ If TRUE will use approximate,
                fast spatial prediction on grids. }
  
   \item{gridRefinement}{Amount to increase the number of grid points
  for the simulation grid.}


\item{giveWarnings}{If true will warn when more than one observation is in a grid box. This is instead of giving an error and stopping.}

  \item{mKrigObject}{An mKrig Object (or spatialProcess object)}

 
  \item{M}{Number of draws from conditional distribution.}
   \item{np}{Degree of nearest neighbors to use. Default \code{np=2} uses 16
   points in a 4X4 grid for prediction of the off grid point. }
  \item{NNSize}{Degree of neighborhood use for fast prediction.
  ( See \code{ predictSurface.mKrig}  with \code{fast= TRUE} ) } 
  \item{nx}{ Number of grid points in prediction locations for x coordinate.}
  \item{ny}{ Number of grid points in  prediction locations for x coordinate.}
  

 

  \item{object}{The spatial fit object.}
  
   \item{predictionGridList}{A grid list specifying the grid locations
   for the conditional samples. }
 

  \item{simulationGridList}{ A \code{gridlist} describing grid for
simulation. If missing this is created from the range of the
locations, \code{nx}, \code{ny}, \code{gridRefinement}, and
\code{gridExpansion}
or from the range and and \code{nxSimulation}, \code{nySimulation}.}

 
  \item{xp}{Same as predictionPoints above.}

 \item{\dots}{Any other arguments to be passed to the predict function.
 	Usually this is the \code{Z} or \code{drop.Z} argument when there are 
 	additional covariates in the fixed part of the model.
 	(See example below.) }

  \item{verbose}{If true prints out intermediate information. }
}

\details{

These functions generate samples from an unconditional or conditional multivariate (spatial) 
distribution, or an approximate one. The \strong{unconditional} simulation function,
\code{simSpatialData},
is a handy way to generate synthetic observations from a fitted model.
Typically one would use these for a parametric bootstrap.
The functions that simulate \strong{conditional} distributions are much
more involved
in their coding. They are useful
for describing the uncertainty in predictions using the estimated spatial 
process under Gaussian assumptions. An important assumption throughout 
these functions is that all covariance parameters are fixed at their 
estimated or prescribed values from the passed object. Although these functions might
be coded up easily by the users these
versions have the advantage that they take the \code{mKrig}, \code{spatialProcess} or
\code{Krig} objects as a way to specify the
model in an unambiguous way. 

Given a spatial process  h(x)= P(x) + g(x) observed at 

Y.k =  Z(x.k)d + P(x.k) + g(x.k) + e.k

where P(x) is a low order, fixed polynomial and g(x) a Gaussian spatial 
process and Z(x.k) is a vector of covariates that are also indexed by space (such as elevation).
Z(x.k)d is a linear combination of the the covariates with the parameter vector d being a
component of the fixed part of the 
model and estimated in the usual way by generalized least squares.

With Y= Y.1, ..., Y.N,
the goal is to sample the conditional distribution of the process. 
 
[h(x) | Y ]  or the full prediction Z(x)d + h(x)

For fixed a covariance this is just a multivariate normal sampling
problem.  \code{sim.Krig.standard} samples this conditional process at
the points \code{xp} and is exact for fixed covariance parameters.
\code{sim.Krig.grid} also assumes fixed covariance parameters and does
approximate sampling on a grid.

The outline of the algorithm is

0) Find the spatial prediction at the unobserved locations based on
 the actual data. Call this h.hat(x) and this is the conditional mean. 

1) Generate a realization that includes both prediction and observation locations from the unconditional spatial process and from this process
simluate synthetic observations. 

2) Use the spatial prediction model ( using the true covariance) to
estimate the spatial process at unobserved locations.

3) Find the difference between the simulated process and its
prediction based on synthetic observations. Call this e(x).

4) h.hat(x) + e(x) is a draw from [h(x) | Y ].

The approximations for this simulation come in at step 1).
Here the field at the observation locations is
approximated using a local conditional simulation from
the nearest grid points. 

NOTE: A fixed part in the model is handled easily by
simply
making the prediction from the synthetic observations
that have mean zero but include estimation of the
fixed part as part of the prediction. 
Because the
regression estimates are unbaised, this gives a
valid draw from the correct mulitvariate distribution
even though the synthetic observations do not include
a fixed part. 


\code{sim.spatialProcess} Follows this algorithm exactly. For the case of an
addtional covariate  this of course needs to be included. For a model
with covariates use \code{drop.Z=TRUE} for the function to ignore prediction
using the covariate and generate conditional samples for just the spatial
process and any low order polynomial. Finally, it should be noted that this
function will also work with an \code{mKrig} object because the essential
prediction information in the mKrig and spatialProcess objects are the same.
The naming is through convenience. 

\code{sim.Krig}  Also follows this algorithm exactly but for the older \code{Krig} object.  Note the inclusion of
drop.Z=TRUE or FALSE will determine whether the conditional simulation 
includes the covariates Z or not. (See example below.) 
\code{simLocal.spatialProcess}{This function is
designed for conditional simulation for a large
prediction grid and for a large number of observation
s. The approximation will be accurate for fine grids
that separate clusters of observation locations. E.g.
multiple observations in a single gridbox are treated
exactly. If observation location are separated by a
grid box then due to the screening effect the
approximation error will be negliable, especially if a
nugget component  (tau ) is present. 

The 1D version of this function will not be much more
efficient than an exact computation. However, it is
included as easy to read source code and for checking
(see examples.)

See the utility function \code{\link{offGridWeights}} for the function that creates weights used to generate the (approximate) conditional sample.  The functions 
\code{checkPredictGrid}, 
\code{makePredictionGridList} and 
\code{makeSimulationGrid} are utilities to setup the grids when not specified. 

}


}
\value{
\code{sim.Krig and sim.spatialProcess}
 a matrix with rows indexed by the locations in \code{xp} and columns being the 
 \code{M} independent draws.


\code{simLocal.spatialProcess} a list with \code{predictionPoints} being the
locations where the field has been simulated.If these have been created from a 
grid list that information is stored in the \code{attributes} of \code{predictionPoints}. 
 \code{Ensemble} is a
matrix where rows index the simulated values of the field and columns
are the different draws, \code{call} is the calling sequence. Not that if \code{predictionPoints}
has been omitted in the call or is created beforehand using \code{make.surface.grid} it is 
easy to reformat the results into an image format for ploting using \code{as.surface}.
e.g. if \code{simOut} is the output object then to plot the 3rd draw:

\preformatted{
     imageObject<- as.surface(simOut$PredictionGrid, simOut$Ensemble[,3] )
     image.plot( imageObject)
}




}

\author{Doug Nychka}
\seealso{ sim.rf, Krig, spatialProcess}
\examples{
## 10 member ensemble for the O3 data

\dontrun{
data( "ozone2")
fitObject<- spatialProcess( ozone2$lon.lat, ozone2$y[16,],
                              smoothness=.5)

nx<- 65
ny<- 55

xGridList<- fields.x.to.grid( fitObject$x, nx=nx, ny=ny)

xGrid<- make.surface.grid( xGridList)

allTime0<- system.time(
 look0<- sim.spatialProcess(fitObject, xp=xGrid, M=5)
)
print( allTime0)

## Local simulation with extra refinement of the grid for embedding
## and same grid size for prediction 
## this runs much faster compared to exact method above 
## as nx, ny are increased  e.g. nx= 128, ny=128 is dramatic difference 

allTime1<- system.time(
  look<- simLocal.spatialProcess(fitObject, M=5,nx=nx, ny=ny,
                     gridRefinement = 3,
                     np=3)
) 
print( allTime1)
print( look$timing)

allTime2<- system.time(
  look<- simLocal.spatialProcess(fitObject, M=5,nx=nx, ny=ny,
                     gridRefinement = 3,
                     np=3, 
                     fast=TRUE)
) 
print( allTime2)
print( look$timing)
}

\dontrun{
## A simple example for setting up a bootstrap 
## M below should be
## set to a much larger sample size, however, ( e.g. M <- 200) for better
## statistics

data( ozone2)
obj<- spatialProcess( ozone2$lon.lat,ozone2$y[16,] )
aHat<- obj$summary["aRange"]
lambdaHat<- obj$summary["lambda"]

######### boot strap 
# create M independent copies of the observation vector
# here we just grab the model information from the 
# spatialProcess object above.
#  
# However, one could just create the list 
#  obj<- list( x= ozone2$lon.lat,
#       cov.function.name="stationary.cov",
#     summary= c( tau= 9.47, sigma2= 499.79, aRange= .700),
#    cov.args= list( Covariance="Matern", smoothness=1.0),
#     weights= rep( 1, nrow(ozone2$lon.lat) )
# )
# Here summary component  has the parameters
# tau, sigma2 and aRange
# and cov.args component has the remaining ones.

set.seed(223)
M<- 25
ySynthetic<- simSpatialData( obj, M)

bootSummary<- NULL

for(  k in 1:M){
cat( k, " ")
# here the MLEs are found using the easy top level level wrapper
# see mKrigMLEJoint for a more efficient strategy
  newSummary<- spatialProcess(obj$x,ySynthetic[,k],
                    cov.params.start= list(
			                   aRange = aHat,
			                  lambda = lambdaHat)
                               )$summary
  bootSummary<- rbind( bootSummary, newSummary)
  }
cat( fill= TRUE)
# the results and 95% confidence interval  
  stats( bootSummary )

  obj$summary
  tmpBoot<- bootSummary[,c("lambda", "aRange") ]
  confidenceInterval <- apply(tmpBoot, 2,
                               quantile, probs=c(0.025,0.975) )
# compare to estimates used as the "true" parameters			       
  obj$summary[2:5] 
  print( t(confidenceInterval) )
# compare to confidence interval using large sample theory  
  print( obj$CITable)
}

\dontrun{
# conditional simulation with covariates
# colorado climate example
  data(COmonthlyMet)
  fit1E<- spatialProcess(CO.loc,CO.tmin.MAM.climate, Z=CO.elev   )
# conditional simulation at missing data
  good<- !is.na(CO.tmin.MAM.climate ) 
  infill<- sim.spatialProcess( fit1E, xp=CO.loc[!good,], 
                Z= CO.elev[!good], M= 10)
# get an elevation grid  ... NGRID<- 50 gives a nicer image but takes longer 
 NGRID <- 25  
 # get elevations on a grid  
   COGrid<- list( x=seq( -109.5, -100.5, ,NGRID), y= seq(36, 41.75,,NGRID) )
   COGridPoints<- make.surface.grid( COGrid)
 # elevations are a bilinear interpolation from the 4km
 # Rocky Mountain elevation fields data set.   
   data( RMelevation)
   COElevGrid<- interp.surface( RMelevation, COGridPoints )
# NOTE call to sim.spatialProcess treats the grid points as just a matrix
# of locations the plot has to "reshape" these into a grid 
# to use with image.plot 
   SEout<- sim.spatialProcess( fit1E, xp=COGridPoints,  Z= COElevGrid, M= 30)
# for just the smooth surface in lon/lat
#  SEout<- sim.spatialProcess( fit1E, xp=COGridPoints,  drop.Z=TRUE, M= 30)
# in practice M should be larger to reduce Monte Carlo error.      
   surSE<- apply( SEout, 2, sd )
   image.plot( as.surface( COGridPoints, surSE)) 
   points( fit1E$x, col="magenta", pch=16) 
   
}

### Approximate conditional simulation 
\dontrun{
 # create larger lon/lat grid  
   NGRID <- 200  
   COGrid<- list( x=seq( -109.7, -100.5, ,NGRID),
                    y= seq(36, 41.75,,NGRID) )
 # interpolation elevations to this grid. 
 # This took about 40 seconds
   COElevGrid<- interp.surface.grid( RMelevation, COGrid )
  system.time( 
    SEout0<- simLocal.spatialProcess( fit1E,COGrid , 
    ZGrid= COElevGrid$z,
    M= 10)
    )
   
  
}
### Approximate conditional simulation and with approximate prediction
### increase  np and NNSize to improve approximations 
### This takes about 8 seconds of course one would want more thatn 10 reps
### to estimate the SE. Use drop.Z=TRUE to just get the spatial surface without ### the fixed part 
\dontrun{
system.time(
   SEout2<- simLocal.spatialProcess( fit1E, COGrid , 
                   ZGrid= COElevGrid$z, np = 2,
   fast= TRUE, NNSize=5, M= 10)
   )
   
  look <- apply( SEout2$z,c(1,2), sd)
  imagePlot(SEout2$x, SEout2$y, look, col=viridis(256) )
  points( fit1E$x, pch=16, cex=.5, col="magenta")
  title("Monte Carlo prediction SE")
}


#### example using Krig object and exact conditional simulation.
\dontrun{
data( ozone2)
set.seed( 399)
# fit to day 16 from Midwest ozone data set.
  out<- Krig( ozone2$lon.lat, ozone2$y[16,], Covariance="Matern", 
            aRange=1.0,smoothness=1.0, na.rm=TRUE)

# NOTE aRange =1.0 is not the best choice but 
# allows the sim.rf circulant embedding algorithm to 
# work without increasing the domain.

#six missing data locations
 xp<-  ozone2$lon.lat[ is.na(ozone2$y[16,]),]

# 5 draws from process at xp given the data 
# this is an exact calculation
 sim.Krig( out,xp, M=5)-> sim.out

}

\dontrun{
## testing the local method on a 1D case.

set.seed(124) 
# 10 observations -- massive dataset!
s<-  cbind(runif(  10, 5,45)) 
y<-  cbind(runif(  10))
aRange<- 10 
obj<- mKrig( s, y, aRange=aRange,Covariance="Matern",
             smoothness=1.0,
             lambda=.01,tau=sqrt(.01),
                      m=0)
#                    
# M should be much larger for an accurate check on code
#

gridList<- list( x= seq(0,50, length.out=15))

look<- simLocal.spatialProcess(obj,
                               np=3,
                               predictionGridList = gridList,
                               gridRefinement = 3,
                               M=50, extrap=TRUE) 
                               
simSE<- apply(look$z, 1, sd )

checkSE<- predictSE( obj, xnew= cbind(look$x ), drop.Z=TRUE )
# percentage error from true SE at each location.
stats( 100*abs(1- simSE/checkSE) )
    
# Maggie plot  
plot( look$x, checkSE, type="b", col="blue",
xlab="Location", ylab="Prediction SE")
rug( look$x, col="blue", lwd=3)
points( look$x, simSE, col="orange3", pch=16)
xline( s, col="grey", lwd=2)
title("Exact (blue) and Monte Carlo (orange) 
for the prediction SE based on observations (grey) ")

}





}
\keyword{spatial}
% at least one, from doc/KEYWORDS
