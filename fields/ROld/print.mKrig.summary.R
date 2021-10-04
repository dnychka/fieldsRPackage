# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2018
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
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
print.mKrigSummary <- function (x, digits = 4, ...) 
{
  cat("Call:\n")
  dput(x$call)
  print(x$summaryTable, quote = FALSE)
  
  # fixed effects are reported differently when fields are replicated.    
  nData<- x$nData
  cat(" ", fill = TRUE)
  if( nData == 1 | x$collapseFixedEffect ){
    cat(" ", fill = TRUE)
    cat("Summary of fixed effects", fill = TRUE)
    print( x$fixedEffectsTable)
  }
  else {
    cat("Estimated fixed effects found separately for each replicate field", 
        fill = TRUE)
  }
  cat(" ", fill = TRUE)
  cat("Covariance Model:", x$cov.function, fill = TRUE)
  if (x$cov.function == "stationary.cov") {
    cat("   Covariance function:  ", ifelse(is.null(x$args$Covariance), 
                                            "Exponential", x$args$Covariance), fill = TRUE)
  }
  if (!is.null(x$args)) {
    cat("   Non-default covariance arguments and their values ", 
        fill = TRUE)
    nlist <- as.character(names(x$args))
    NL <- length(nlist)
    for (k in 1:NL) {
      cat("   Argument:", nlist[k], " ")
      if (object.size(x$args[[k]]) <= 1024) {
        cat("has the value(s): ", fill = TRUE)
        print(x$args[[k]])
      }
      else {
        cat("too large to print value, size > 1K ...", 
            fill = TRUE)
      }
    }
  }
  invisible(x)
}