#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2022 Colorado School of Mines
# 1500 Illinois St., Golden, CO 80401
# Contact: Douglas Nychka,  douglasnychka@gmail.edu,
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
print.spatialProcessSummary <- function(x, digits = 4, ...) {
  cat("CALL:\n")
  dput(x$call)
  cat("\n")
  cat(" SUMMARY OF MODEL FIT:\n")
  print(x$summaryTable, quote = FALSE)
  cat("\n")
  if( !is.na(x$fixedEffectsTable[1])){
  cat(" ESTIMATED COEFFICIENTS FOR FIXED PART:", fill = TRUE)
  cat("\n")
  print(x$fixedEffectsTable, quote = FALSE)
  cat("\n")
  }
  if( !is.na(x$fixedEffectsTableCommon[1])){
  cat(" ESTIMATED COEFFICIENTS FOR COMMON COVARIATES:", fill = TRUE)
  cat("\n")
  print(x$fixedEffectsTableCommon, quote = FALSE)
  cat("\n")
  }
  cat(" COVARIANCE MODEL:", x$cov.function, fill = TRUE)
  if (x$cov.function == "stationary.cov") {
    covName<- ifelse(is.null(x$args$Covariance), 
                     "Exponential", x$args$Covariance)
    cat("  Covariance function: ", covName, fill = TRUE)
  }
  if (!is.null(x$args)) {
    cat("   Non-default covariance arguments and their values ", 
        fill = TRUE)
    nlist <- as.character(names(x$args))
    NL <- length(nlist)
    for (k in 1:NL) {
      covItem<- x$args[[k]]
      cat( nlist[k], ":",
           fill = TRUE)
      if (object.size(covItem) > 10) {
        print(covItem)
      }
      else {
        cat("Too large to print value, size > 10 ...", 
            fill = TRUE)
      }
    }
  }
  cat( "Nonzero entries in covariance matrix ", x$nonzero.entries, fill=TRUE)
  
  cat("\n")
  cat("SUMMARY FROM Max. Likelihood ESTIMATION:", fill=TRUE)
  
  if( !is.null( x$MLEInfo) ){
    cat("Parameters found from optim: ", fill=TRUE )
    print( x$MLESummary[x$MLEpars] )
    cat("Approx. confidence intervals for MLE(s) ", fill=TRUE )
    print( x$CITable)
    
    cat("\n")
    cat(" Note: MLEs for  tau and sigma found analytically from lambda", fill=TRUE)
    cat("\n")
  }
  else{
    cat(" lambda and range are supplied : ", fill=TRUE)
    cat(" tau and sigma  analytically  derived from lambda", fill=TRUE)
    cat("\n")
  }
  cat("Summary from estimation:", fill=TRUE)
  print( x$MLESummary)
 
  invisible(x)
}
