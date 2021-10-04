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
"arrow.plot" <- function(a1, a2, u = NA, v = NA, arrow.ex = 0.05, 
    xpd = TRUE, true.angle = FALSE, arrowfun = arrows, ...) {
    if (is.matrix(a1)) {
        x <- a1[, 1]
        y <- a1[, 2]
    }
    else {
        x <- a1
        y <- a2
    }
    if (is.matrix(a2)) {
        u <- a2[, 1]
        v <- a2[, 2]
    }
    ucord <- par()$usr
    arrow.ex <- arrow.ex * min(ucord[2] - ucord[1], ucord[4] - 
        ucord[3])
    if (true.angle) {
        pin <- par()$pin
        r1 <- (ucord[2] - ucord[1])/(pin[1])
        r2 <- (ucord[4] - ucord[3])/(pin[2])
    }
    else {
        r1 <- r2 <- 1
    }
    u <- u * r1
    v <- v * r2
    maxr <- max(sqrt(u^2 + v^2))
    u <- (arrow.ex * u)/maxr
    v <- (arrow.ex * v)/maxr
    invisible()
    old.xpd <- par()$xpd
    par(xpd = xpd)
    arrowfun(x, y, x + u, y + v, ...)
    par(xpd = old.xpd)
}
