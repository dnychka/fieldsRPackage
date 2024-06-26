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
vgram.matrix <- function(dat, R = NULL, dx = NULL, dy = NULL) {

# a useful function for matching shifted indices
# (the kind of internal function Dorit does not like!)
    SI <- function(ntemp, delta) {
        n1 <- 1:ntemp
        n2 <- n1 + delta
        good <- (n2 >= 1) & (n2 <= ntemp)
        return(cbind(n1[good], n2[good]))
     }
#  
# check if dat is in $x $y $z image format
# if so extract z and find dx and dy
if( is.list( dat)){
    if(is.null(dat$z)){
        stop("need a z matrix in dat")
    }
    if( !is.null(dat$x)& is.null(dx)){
        dx<- dat$x[2]-dat$x[1]
        if( any( diff( dat$x)-dx > 10 *.Machine$double.eps) ){
            stop(" x grid values are not equal")
        }
       
    }
    if( !is.null(dat$y)& is.null(dy)){
        dy<- dat$y[2]-dat$y[1]
        if( any( diff( dat$y)-dy > 10 *.Machine$double.eps) )
            {
            stop(" x grid values are not equal")
        }
    }
    dat <- dat$z
}
# if dat is not an image defaults for dx and dy are 1.0    
    if( is.null( dx)){
        dx<- 1
    }
    if( is.null( dy)){
        dy<- 1
    }
    # if R not specified go out to 5 grid boxes in distance
    # 
    if( is.null(R)){
        R<-  5* max( c(dx, dy))
    }
    
    M<- nrow(dat)
    N<- ncol( dat)
    # create all possible separations for a grid up to a distance R
    m <- min( c(round(R/dx),M) )
    n <- min( c(round(R/dy),N) )
    #
    # all relavent combinations:  note that negative increments are
    # needed as well as  positive ones
    ind <- rbind(as.matrix(expand.grid(0, 1:n)), as.matrix(expand.grid(1:m, 
        0)), as.matrix(expand.grid(c(-(m:1), 1:m), 1:n)))
    # distances - only take those within a distance R.
    # and trim everything to this bound
    d <- sqrt((dx * ind[, 1])^2 + (dy * ind[, 2])^2)
    good <- (d > 0) & (d <= R)
    ind <- ind[good, ]
    d <- d[good]
    ind <- ind[order(d), ]
    d <- sort(d)
    #
    # arrays to hold statistics
    nbin <- nrow(ind)
    holdVG <- rep(NA, nbin)
    holdRVG <- rep(NA, nbin)
    holdN <- rep(0, nbin)
    # loop over each separation
    for (k in 1:nbin) {
        # indices for original and shifted image that are within array bounds
        MM <- SI(M, ind[k, 1])
        NN <- SI(N, ind[k, 2]) 
        if( length(MM)>0 & length(NN)>0){     
        # find differences       
          BigDiff <- (dat[MM[, 1], NN[, 1] ] - dat[MM[, 2], NN[,2] ] )  
        # standard and the Cressie robust version.
        # modified to handle NAs
          holdVG[k] <-  mean(0.5 * (BigDiff)^2, na.rm = TRUE)
          holdRVG[k] <- mean(abs(BigDiff)^0.5, na.rm = TRUE)
          holdN[k]   <-   sum( !is.na(BigDiff) ) 
        }
    }
    # finish robust estimate Cressie (1993) formula 2.4.12
    holdRVG <- 0.5 * (holdRVG^4)/(0.457 + 0.494 * holdN)
    # collapsed variogram to common distances this what one would look
    # at under the stationary case.
    top <- tapply(holdVG * holdN, d, FUN = "sum")
    bottom <- tapply(holdN, d, FUN = "sum")
    dcollapsed <- as.numeric(names(bottom))
    vgram <- top/bottom
    #  wipe out pesky row names
    dimnames(vgram) <- NULL
    out <- list(vgram = vgram, d = dcollapsed, ind = ind, d.full = d, 
        vgram.full = holdVG, robust.vgram = holdRVG, N = holdN, 
        dx = dx, dy = dy)
    class(out) <- "vgram.matrix"
    return(out)
}
