#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2021 Colorado School of Mines
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
envelopePlot <- function(x1, y1, x2 = x1, y2,
                         col ="thistle1" , lineCol = "thistle3", ...) {
#  sort the curves -- just in case they are passed out of order
    ind<- order( x1)
    x1<- x1[ind]
    y1<- y1[ind]
    ind<- order( x2)
    x2<- x2[ind]
    y2<- y2[ind]
   
  polygon(c(x1, rev(x2)), c(y1, rev(y2)), col = col, border = NA, ...)
  lines(x1, y1, lwd = 3, col = lineCol)
  lines(x2, y2, lwd = 3, col = lineCol)
}
