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

\name{predictSE}
\alias{predictSE}
\alias{predictSE.Krig}
\alias{predictSE.mKrig}
\alias{predictSEUsingKrigA}
\title{
  Standard errors of predictions for Krig spatial process estimate 
}
\description{
Finds the standard error ( or covariance) of prediction based on a linear
combination of
the observed data. The linear combination is usually the "Best Linear
Unbiased Estimate" (BLUE) found from the Kriging equations. 
This statistical computation is done under the assumption that the
covariance function is known.
}
\usage{
predictSE(object, ...) 
\method{predictSE}{Krig}(object, x = NULL, cov = FALSE, verbose = FALSE,...)
\method{predictSE}{mKrig}(object, xnew = NULL, XMat = NULL, verbose = FALSE, drop.XMat
                 = FALSE, ...)
 }
\arguments{
\item{drop.XMat}{If FALSE find standard error without including the additional spatial covariates
described by \code{XMat}. If TRUE find full standard error with spatial covariates if they are part of the model.}
\item{object}{ A fitted object that can be used to find prediction standard errors. This is usually from fitting a spatial model to data. e.g. a Krig or mKrig object.  
}
\item{xnew}{
Points to compute the predict standard error or the prediction 
cross covariance matrix.  
}
\item{x}{
Same as \code{xnew} -- points to compute the predict standard error or the prediction 
cross covariance matrix.  
}
\item{cov}{
If TRUE the full covariance matrix for the predicted values is returned.
Make sure this will not be big if this option is used. ( e.g. 50X50 grid
will return a matrix that is 2500X2500!)  If FALSE just the marginal
standard deviations of the predicted values are returned. Default is
FALSE -- of course.
}

\item{verbose}{If TRUE will print out various information for debugging.}

\item{\dots}{
These additional arguments passed to the predictSE function.
}
\item{XMat}{Additional matrix of spatial covariates used for prediction. These are used to 
determine the additional covariance contributed in teh fixed part of the model.}

}

\value{
A vector of standard errors for the predicted values of the Kriging fit. 
}
\details{
The predictions are represented as a linear combination of the dependent 
variable, Y. Call this LY.  Based on this representation the conditional 
variance is the same as the expected value of (P(x) + XMat(X) - LY)**2.
where 
P(x)+XMat(x) is the value of the surface at x and LY is the linear 
combination that estimates this point. Finding this expected value is 
straight forward given the unbiasedness of LY for P(x) and the covariance 
for XMat and Y. 

In these calculations it is assumed that the covariance parameters are fixed.
This is an approximation since in most cases they have been estimated from the 
data. It should also be noted that if one assumes a Gaussian field and known 
parameters in the covariance, the usual Kriging estimate is the 
conditional mean of the field given the data. This function finds the 
conditional standard deviations (or full covariance matrix) of the 
fields given the data.  

 There are two useful extensions supported by this function.  Adding the 
variance to the estimate of the spatial mean if this is a correlation 
model. (See help file for Krig) and calculating the variances under 
covariance misspecification.  The function \code{predictSE.KrigA} uses 
the smoother matrix ( A(lambda) ) to find the standard errors or 
covariances directly from the linear combination of the spatial 
predictor. Currently this is also the calculation in 
\code{predictSE.Krig} although a shortcut is used 
\code{predictSE.mKrig} for mKrig objects.

}
\seealso{
spatialProcess, mKrig, Krig, predict.Krig, predictSurfaceSE, 
}
\examples{
# 
 data(ozone2)
 sOzone<- ozone2$lon.lat
 yOzone<- ozone2$y[16,]
# Matern smoothness =1.0 sigma2 and tau found by ML
  fit0<- spatialProcess( sOzone,yOzone, aRange=10)
# prediction SEs at  miss data locations
  missingObs<- is.na(yOzone)
  missingSE<- predictSE(fit0, sOzone[missingObs,]) 
# NOTE based on dispatching predictSE.mKrig is called. 
# 95% CI for  ozone surface
   fHat<- predict(fit0, sOzone[missingObs,]) 
   cbind(fHat - 1.96*missingSE, fHat + 1.96*missingSE )
# 95% CI for  ozone _observation_ 
   fHat<- predict(fit0, sOzone[missingObs,])
   tau2<- fit0$summary["tau"]^2
   cbind(fHat - 1.96*sqrt(missingSE + tau2) ,
         fHat + 1.96*sqrt(missingSE + tau2) )

# SE on  a  grid of locations  
  xg<-make.surface.grid( 
  list(x=seq(-91,-83,,30), 
       y=seq( 37, 45,,30))
  )     
  gridSE<- predictSE(fit0,xg)   # std errors of predictions 
  
#at the grid points out is a vector of length 900 
#reshape the grid points into a 30X30 matrix etc.  

   out.p<-as.surface( xg, gridSE) 
   surface( out.p, type="C") 
   points( sOzone[!missingObs,], col="magenta", pch=16)

# this is equivalent to the higher level functions
# (We are lazy and the grid will be a bit different here. See the 
# help file!)
 out<- predictSurfaceSE( fit0, nx=30, ny=30, extra=TRUE) 
 imagePlot( out, col=viridis(256))
 contour( out, add=TRUE, levels=c( 10,20, 20))
 points( sOzone[!missingObs,], col="magenta", pch=16)


}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
