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
"summaryGCV.Krig" <- function(object, lambda, cost = 1, 
    verbose = FALSE, offset = 0, y = NULL, ...) {
    out <- object
    nt <- out$nt
    np <- out$np
    N <- out$N
    D <- out$matrices$D
    if (is.null(y)) {
        u <- out$matrices$u
        tauHat.pure.error <- out$tauHat.pure.error
        pure.ss <- out$pure.ss
    }
    else {
        out2 <- Krig.coef(out, y)
        u <- out2$u
        tauHat.pure.error <- out2$tauHat.pure.error
        pure.ss <- out2$pure.ss
    }
    info <- list(matrices = list(D = D, u = u), N = N, nt = nt, 
        cost = cost, pure.ss = pure.ss, tauHat.pure.error = tauHat.pure.error, 
        offset = offset)
    if (verbose) {
        print(info)
    }
    lambda.est <- rep(NA, 6)
    names(lambda.est) <- c("lambda", "trA", "GCV", "GCV.one", 
        "GCV.model", "tauHat")
    lambda.est[1] <- lambda
    lambda.est[2] <- Krig.ftrace(lambda, D)
    lambda.est[3] <- Krig.fgcv(lambda, info)
    lambda.est[4] <- Krig.fgcv.one(lambda, info)
    if (!is.na(tauHat.pure.error)) {
        lambda.est[5] <- Krig.fgcv.model(lambda, info)
    }
    lambda.est[6] <- sqrt(Krig.fs2hat(lambda, info))
    lambda.est
}
