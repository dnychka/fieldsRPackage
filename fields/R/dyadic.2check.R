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
dyadic.2check <- function(m, n, cut.p = 2) {
    # checks that n is of the form
    # n=p*2^m where p <= cut.p
    m2 <- as.integer(m)
    n2 <- as.integer(n)
    while ((n2 > cut.p) & (m2 > cut.p)) {
        if ((m2%%2 != 0) | (n2%%2 != 0)) {
            cat(n, "and", m, "must equal p*2^L where p is less than or equal to ", 
                cut.p, fill = TRUE)
            return(FALSE)
        }
        m2 <- m2/2
        n2 <- n2/2
    }
    return(TRUE)
}
