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



bplot <- function(x, by, pos = NULL, at = pos, add = FALSE, 
    boxwex = 0.8, xlim = NULL, ...) {
    if (!missing(by)) {
        x <- split(c(x), as.factor(by))
    }
    if (!add & !is.null(at) & is.null(xlim)) {
        xlim <- range(at)
    }
    if (!is.null(at)) {
        boxwex <- boxwex * min(diff(sort(at)))
    }
    boxplot(x, at = at, xlim = xlim, add = add, boxwex = boxwex, 
        ...)
    
}

bplot.xy<- function (x, y, N = 10, breaks = NULL, 
                        plot = TRUE,axes=TRUE,  ...) 
{
    if(is.null(breaks)){
      breaks<- pretty(x, N, eps.correct = 1)
    }
    # N ignored if breaks are specified. 
    NBIN <- length(breaks) - 1
    centers <- (breaks[1:NBIN] + breaks[2:(NBIN + 1)])/2
    obj <- split(y, cut(x, breaks))
    if (length(obj) == 0) {
        stop("No points within breaks")
    }
    if (plot) {
        bplot(obj, at = centers, show.names = FALSE, axes = axes, 
              ...)
        if( axes){ 
            axis(1)
        }
    }
    else {
        return(list(centers = centers, breaks = breaks, boxplot.obj = boxplot(obj, 
                                                                              plot = FALSE)))
    }
}


