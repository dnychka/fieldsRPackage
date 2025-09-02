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
summary.spatialProcess <- function(object, ...) {
# output list
  outObject<- list()
  digits<- 4
  if (is.matrix(object$residuals)) {
    n <- nrow(object$residuals)
    nData <- ncol(object$residuals)
  }
  else {
    n <- length(object$residuals)
    nData <- 1
  }
  
  c1 <- "Number of Observations:"
  c2 <- n
  
  if (nData > 1) {
    c1 <- c(c1, "Number of data sets fit:")
    c2 <- c(c2, nData)
  }
  
  c1 <- c(c1, "Degree of polynomial in fixed part: ")
  
  
  if(object$m !=0 ){
    c2 <- c(c2, object$m - 1)
  }
  else{
    c2 <- c(c2, NA)
  }
  c1 <- c(c1, "Total number of parameters in fixed part: ")
  c2 <- c(c2, object$nt)
  
  if (object$nZ > 0) {
    c1 <- c(c1, "Number of additional covariates (Z)")
    c2 <- c(c2, object$nZ)
  }
  
  if ( !object$simpleKriging ){
    if (!is.na(object$gamma[1])) {
    c1 <- c(c1, "Number of common covariates (ZCommon)")
    c2 <- c(c2, length(object$gamma))
    }
      }
  
   
  c1 <- c(c1, "sigma Process stan. dev: ")
  c2 <- c(c2, signif( sqrt(object$MLESummary["sigma2"]), digits))
  
  c1 <- c(c1, "tau  Nugget stan. dev:")
  c2 <- c(c2, signif(object$MLESummary["tau"], digits))
  
  c1 <- c(c1, "lambda   tau^2/sigma^2: ")
  c2 <- c(c2, signif(object$MLESummary["lambda"], digits))
  
  c1 <- c(c1, "aRange parameter (in units of distance): ")
  c2 <- c(c2, signif(object$MLESummary["aRange"], digits))
  
  if (!is.na(object$eff.df)) {
    c1 <- c(c1, "Approx.  degrees of freedom for curve")
    c2 <- c(c2, signif(object$eff.df, digits))
    if (length(object$trA.info) < object$np) {
      c1 <- c(c1, "   Standard Error of df estimate: ")
      c2 <- c(c2, signif(sd(object$trA.info)/sqrt(length(object$trA.info)), 
                         digits))
    }
  }
  
  c1 <- c(c1, "log Likelihood: ")
  c2 <- c(c2, object$summary["lnProfileLike.FULL"])
  c1 <- c(c1, "log Likelihood REML: ")
  c2 <- c(c2, object$summary["lnProfileREML.FULL"])
  
  summaryStuff <-  cbind(c1, c2)
  dimnames(summaryStuff) <- list(rep("",
                                     dim(summaryStuff)[1]),
                                 rep("", dim(summaryStuff)[2]))
  ###########
  outObject$summaryTable <- summaryStuff
  outObject$collapseFixedEffect <- object$collapseFixedEffect
  outObject$CITable <- object$CITable
  ###########
  if (!is.null(object$MLEInfo)) {
    outObject$MLEpars <-  names(object$MLEInfo$pars.MLE)
    outObject$MLESummary <- object$summary
  }
  else{
    outObject$MLEpars <- NA
    outObject$MLESummary <- object$summary
  }
  ########### information for SE for fixed effects
  if (!is.null(object$beta)){
    if (outObject$collapseFixedEffect | (nData == 1)) {
      outObject$fixedEffectsCov <- object$fixedEffectsCov
      SE <- sqrt(diag(outObject$fixedEffectsCov))
      beta <-  object$beta[, 1]
      SEbeta <- SE[1:length(beta)]
      pValue <- pnorm(abs(beta / SEbeta), lower.tail = FALSE) * 2
      outObject$fixedEffectsTable <- cbind(signif(beta, digits),
                                           signif(SEbeta, digits),
                                           signif(pValue, digits))
      # get row names
      if (is.null(object$fixedEffectNames)) {
        outObject$fixedEffectNames <- paste0("d", 1:(object$nt))
      }
      else{
        outObject$fixedEffectNames <- object$fixedEffectNames
      }
      dimnames(outObject$fixedEffectsTable) <-
        list(outObject$fixedEffectNames,
             c("estimate", "SE", "pValue"))
      
    }
    else{
      # if more that one realization just summarize the coefficients
      outObject$fixedEffectsTable <- stats(t(object$beta))
    }
  }
  if (is.null(object$beta)) {
    outObject$fixedEffectsTable <- NA
     }
  
  if (!is.null(object$gamma)) {
    gamma <- object$gamma
    
    if( object$collapseFixedEffect){
    SE <- sqrt(diag(object$fixedEffectsCov))
    SEgamma <- SE[(1:length(gamma)) + length(beta)]
    }
    else{
      SEgamma <- sqrt(diag(object$fixedEffectsCovCommon))
    }
    
    pValue <- pnorm(abs(gamma / SEgamma), lower.tail = FALSE) * 2
    outObject$fixedEffectsTableCommon <-
      cbind(signif(gamma, digits),
            signif(SEgamma, digits),
            signif(pValue, digits))
    tmp <- paste0("gamma", 1:length(gamma))
    dimnames(outObject$fixedEffectsTableCommon) <- list(tmp,
                                                        c("estimate", "SE", "pValue"))
  }
  
  if (is.null(object$gamma)) {
    outObject$fixedEffectsTableCommon <- NA
  }
  
  
  #####################
  outObject$nData <- nData
  outObject$call <- object$call
  outObject$cov.function <- object$cov.function
  outObject$args <- object$args
  outObject$nonzero.entries <- object$nonzero.entries
  outObject$MLEInfo <- object$MLEInfo
  
  class(outObject) <- "spatialProcessSummary"
  
  return(outObject)
  
}
