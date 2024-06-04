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
"Exponential" <- function(d, aRange=1.0,  phi = 1.0,
                          theta=NULL, range=NULL) {
    #
    # Matern covariance function transcribed from Stein's book page 31
    # nu==smoothness==.5, alpha ==  1/range
    # 
    # GeoR parameters map to kappa==smoothness and phi == range
    # to make old code from Fuentes and also the package SEHmodel work
    # phi is accepted as the marginal variance of the process (see below)
    # within fields this parameter is "sigma2" and is a standard notation
                               #
    # check for negative distances
    #
    # fill in other arguments as aRange.
    # Note that the theta argument has been depreciated.
    if( !is.null(theta)){  
     aRange<- theta
    }
    if( !is.null(range)){  
     aRange<- range
    }
    if (any(d < 0)) 
        stop("distance argument must be nonnegative")
    #
    return(phi*exp(-d/aRange))
}
