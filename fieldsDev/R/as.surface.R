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
as.surface <-
  function(obj,
           z,
           location = NULL,
           order.variables = "xy") {
    #
    if (is.list(obj)) {
      grid.list <- obj
    }
    if (is.matrix(obj)) {
      grid.list <- attr(obj, "grid.list")
    }
    #
    #  OK now have a grid, parse this to figure
    #  nx and ny the x and y sequences and extract names
    #  if a 1D grid then hold$ny is NA
    
    hold <- parse.grid.list(grid.list, order.variables = "xy")
    
    grid1D <- is.na(hold$ny)
    
    # now fill in the "z" values -- the actual grid values
    if (!grid1D) {
      if (!is.null(location)) {
        # location is a logical or two column matrix of indices
        # with the position of the z values
        temp <- rep(NA,  hold$ny * hold$nx)
        temp[location] <- z
        temp <- matrix(temp, ncol = hold$ny, nrow = hold$nx)
      }
      else{
        if (length(c(z)) != hold$ny * hold$nx) {
          print(length(c(z)))
          print(c(hold$ny * hold$nx))
          stop("z does not match number of grid points")
        }
        temp <- matrix(z,    ncol = hold$ny, nrow = hold$nx)
      }
    }
    else{
      temp <- matrix(z,    ncol = 1, nrow = hold$nx)
    }
    
    #
    # note that coercing z to a matrix is just reformatting
    # using the standard ordering.
    #
    # output list is all the grid stuff and the matrix z.
    return( c(hold,
              list(z = temp),
              list(grid1D = grid1D))
            )
  }
