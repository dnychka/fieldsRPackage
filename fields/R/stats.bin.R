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
"stats.bin" <- function(x, y, N = 10, breaks = NULL, prettyBins=FALSE) {
    out <- list()
    
    if (is.null(breaks)) {
      if(prettyBins){
        breaks <- pretty(x, N)
      }
      else{
        breaks<- seq( min( x), max(x), length.out=N)
      }
    }
    NBIN <- length(breaks) - 1
    centers <- (breaks[1:NBIN] + breaks[2:(NBIN + 1)])/2
    test <- describe()
    obj <- matrix(NA, ncol = NBIN, nrow = length(test))
    dimnames(obj) <- list(test, format(1:NBIN))
    obj[, 1] <- describe(y[x <= breaks[2] & x >= breaks[1]])
    for (k in 2:NBIN) {
        obj[, k] <- describe(y[x <= breaks[k + 1] & x > breaks[k]])
    }
    list(centers = centers, breaks = breaks, stats = obj)
}
