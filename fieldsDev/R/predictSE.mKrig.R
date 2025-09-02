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
predictSE.mKrig<- function(object, xnew = NULL, 
                              Z = NULL, verbose = FALSE,
                              drop.Z = FALSE, ...) {
  #
  # name of covariance function
  call.name <- object$cov.function.name
  # check for collapsed fixed effects
  
  # if(drop.Z){
  #   stop("  Sorry drop.Z not supported")
  # }
  #
  # default is to predict at data x's
  if (is.null(xnew)) {
    xnew <- object$x
  }
  if ( is.null(Z)& !is.null( object$Z) ) {
    Z <- object$Z
  }
  xnew <- as.matrix(xnew)
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
  }
  if (verbose) {
    cat("Passed locations", fill=TRUE)
    print(xnew)
    cat("Assumed covariates", fill=TRUE)
    print(Z)
  }
  objectSummary<- object$summary
  lambda <- objectSummary["lambda"]
  sigma2 <- objectSummary["sigma2"]
  tau2 <- objectSummary["tau"]^2
  if (verbose) {
    print(c(lambda, sigma2, tau2))
  }
  k0 <- do.call(call.name, c(object$args, list(x1 = object$x, 
                                               x2 = xnew)))
  #
  # old form based on the predict function
  #   temp1 <-  sigma2*(t0%*% object$Omega %*%t(t0)) -
  #          sigma2*predict( object, y= k0, x=x) -
  #          sigma2*predict( object, y= k0, x=x, just.fixed=TRUE)
  # alternative formula using the beta and c.coef coefficients directly.
  # collapseFixedEffect=FALSE because
  # we want the "fixed effect" computation
  # to be done separately for each column of k0
  hold <- mKrig.coef(object, y = k0, collapseFixedEffect=FALSE)
  #print( hold$beta)
  # temp0 are marginal variances -- trival in the stationary case!
  temp0<- sigma2 * do.call(call.name,
                            c(object$args, 
                              list(x1 = xnew, marginal = TRUE))) 
  temp <- temp0 - sigma2 *colSums((k0) * hold$c.coef)
  # add contribution from fixed effect covariates 
  if( !drop.Z){
    # fixed effects matrix includes both spatial drift and covariates.
    t0 <- t(cbind(fields.mkpoly(xnew, m = object$m), Z))
    temp1 <- sigma2 * (colSums(t0 * (object$Omega %*% t0)) 
                       - 2 * colSums(t0 * hold$beta))
    #print( colSums(t0 * (object$Omega %*% t0)))
    temp <- temp + temp1
  }
  # return square root as the standard error in units of observations.
  return(sqrt(temp))
}
