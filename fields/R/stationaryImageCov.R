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
stationaryImageCov <- function(
                               
                               ind1=NULL, 
                               ind2=NULL,
                               Y, 
                               covObject = NULL,
                               gridList=NULL, 
                              mKrigObject = NULL,
    setup = FALSE,
    M = NULL,
    cov.function="stationary.cov",
    delta=NULL, 
    cov.args=NULL, 
    ...) {
    
   
    if (is.null(covObject)) {
  # sort out where the covariance model is 
      if( is.null(mKrigObject)){
        cov.args<-c( cov.args, list(...))
      }
    covObject<- circulantEmbeddingSetup( 
                                  grid = gridList, 
                                     M = M, 
                           mKrigObject = mKrigObject, 
                          cov.function = "stationary.cov",
                              cov.args = cov.args,
                                 delta = NULL,
                                giveWarnings=FALSE)
                                         
    if( setup){
    return(covObject)
    }
    }
    temp <- matrix(0, nrow = covObject$M[1], ncol = covObject$M[2])
    if (is.null(ind1)) {
        temp[1:covObject$m[1], 1:covObject$m[2]] <- Y
        indexSubset<- list( x= 1:covObject$m[1], y=1:covObject$m[2] )
        return( 
        Re(fft(fft(temp) * covObject$wght, inverse = TRUE)[
          indexSubset$x,indexSubset$y])
        )
    }
    else {
        if (missing(ind2)) {
            temp[ind1] <- Y
        }
        else {
            temp[ind2] <- Y
        }
        #
        # as promised this is a single clean step
        #
        Re(fft(fft(temp) * covObject$wght, inverse = TRUE)[ind1])
    }
}
