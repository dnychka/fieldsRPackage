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
"describe" <- function(x) {
    lab <- c("N", "mean", "Std.Dev.", "min", "Q1", "median", 
        "Q3", "max", "missing values")
    if (missing(x)) {
        return(lab)
    }
    temp <- rep(0, length(lab))
    xt <- x[!is.na(x)]
    ix <- order(xt)
    n <- length(xt)
    if (!is.numeric(xt) || all(is.na(x))) {
        return(c(n, rep(NA, length(lab) - 2), length(x) - length(xt)))
    }
    if (n == 1) {
       tmp<- c(n, xt[1], NA, rep(xt[1], 5), length(x) - length(xt))
       names( tmp)<- lab
        
    }
    else {
        tmp<- c(n, mean(xt), sqrt(var(xt)), min(xt), quantile(xt, 
            c(0.25, 0.5, 0.75)), max(xt), 
            length(x) - length(xt))
      names( tmp)<- lab
      }
    
    return( tmp)
}
