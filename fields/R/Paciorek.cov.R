#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2022 Colorado School of Mines
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
"Paciorek.cov" <- function(x1,
                           x2 = NULL,
                           Distance = "rdist",
                           Dist.args = NULL,
                           aRangeObj = 1,
                           sigma2Obj = NULL,
                           sigma2 = 1,
                           C = NA,
                           marginal = FALSE,
                           smoothness = .5)
{
  # get covariance function arguments from call
  if (class(aRangeObj)[1] == "numeric") {
    class(aRangeObj) <- "constant"
  }
  
  # coerce x1 and x2 to matrices
  if (is.data.frame(x1))
    x1 <- as.matrix(x1)
  if (!is.matrix(x1))
    x1 <- matrix(c(x1), ncol = 1)
  
  if (is.null(x2))
    x2 = x1
  if (is.data.frame(x2))
    x2 <- as.matrix(x2)
  if (!is.matrix(x2) & !is.null(x2))
    x2 <- matrix(c(x2), ncol = 1)
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  #
  if (!marginal) {
    aRangex1 <- c(predict(aRangeObj, x1))
    aRangex2 <- c(predict(aRangeObj, x2))
    # this is (Sigma_i + Sigma_j )/2
    aRange2Matrix <- (outer(aRangex1 ^ 2, aRangex2 ^ 2, '+') / 2)
    
    distMat <- do.call(Distance, c(list(x1 = x1,
                                        x2 = x2), Dist.args))
    
    dimX<- ncol( x1) # adjust for dimension.
    #  detS1= (determinant Sigma_i) ^1/4 
    detS1<- (aRangex1)^(2*dimX/4)
    detS2<- (aRangex2)^(2*dimX/4)
    # detS12 = (determinant (Sigma_i + Sigma_j)/2 )  ^1/2
    detS12 <- aRange2Matrix^(dimX/2)
    normCov <-outer(detS1, detS2, "*")/detS12
  
    covMat <- normCov * Matern(
                         distMat / sqrt(aRange2Matrix), 
                       smoothness = smoothness) 
    

    
    if (!is.null(sigma2Obj)) {
      cat( "sigma Obj used", fill=TRUE)
      sigmax1 <- c(sqrt(predict(sigma2Obj, x1)))
      sigmax2 <- c(sqrt(predict(sigma2Obj, x2)))
      covMat <-  t(t(sigmax1 * covMat) * sigmax2)
    }
    else{
      covMat <-  sigma2 * (covMat)
    }
    
    if (is.na(C[1])) {
      # distMat is a full matrix
      return(covMat)
    }
    # or multiply cross covariance by C
    # as coded below this is not particularly efficient of memory
    else{
      return(covMat %*% C)
    }
  }
  else{
    if (!is.null(sigma2Obj)) {
      marginalVariance <- predict(sigma2Obj, x1)
    }
    else{
      marginalVariance <- rep(sigma2, ncol(x1))
    }
    return(marginalVariance)
  }
  
}
