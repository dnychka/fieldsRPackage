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
"summaryGCV.sreg" <- function(object, lambda, cost = 1, 
    nstep.cv = 20, offset = 0, verbose = TRUE, ...) {
    out <- object
    tauHat.pure.error <- out$tauHat.pure.error
    pure.ss <- out$pure.ss
    nt <- 2
    np <- out$np
    N <- out$N
    out$cost <- cost
    out$offset <- offset
    lambda.est <- rep(NA, 6)
    names(lambda.est) <- c("lambda", "trA", "GCV", "GCV.one", 
        "GCV.model", "tauHat")
    #
    # fill in stuff for this  lambda
    lambda.est[1] <- lambda
    temp <- sreg.fit(lambda, out)
    lambda.est[2] <- temp$trace
    lambda.est[3] <- temp$gcv
    lambda.est[4] <- temp$gcv.one
    if (!is.na(tauHat.pure.error)) {
        lambda.est[5] <- temp$gcv.model
    }
    lambda.est[6] <- temp$tauHat
    lambda.est
}
