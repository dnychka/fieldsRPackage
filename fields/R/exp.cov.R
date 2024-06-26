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
"Exp.cov" <- function(x1, x2=NULL, aRange = 1, p=1, 
                      distMat = NA, C = NA, marginal = FALSE, onlyUpper=FALSE,
                      theta=NULL, ...)  {
  # theta argument has been depreciated.
  if( !is.null( theta)){
    aRange<- theta
  }
  
  if (!is.matrix(x1)) 
    x1 <- as.matrix(x1)
  if (is.null(x2)) 
    x2 <- x1
  if (!is.matrix(x2)) 
    x2 <- as.matrix(x2)
  if (length(aRange) > 1)
    stop("Non-scalar aRange as input to Exp.cov is depracated.  Use the V argument in stationary.cov or scale
         the input locations beforehand.")
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  # scale the coordinates by aRange if distance matrix isn't precomputed
  # a more general scaling by a matrix is done in stationary.cov
  if(is.na(distMat[1]) || !is.na(C[1])) {
    x1 <- x1*(1/aRange)
    x2 <- x2*(1/aRange)
  }
  #
  # there are three main possible actions listed below:
  #
  # if no cross covariance matrix and marginal variance not desired
  if (is.na(C[1]) && !marginal) {
    
    #compute distance matrix if necessary
    if(is.na(distMat[1]))
      distMat = rdist(x1, x2, compact=TRUE)
    else
      distMat = distMat*(1/aRange)
    
    #only exponentiate by p if p != 1
    if(p != 1)
      distMat = distMat^p
    
    if(inherits(distMat, "dist")) {
      #distMat is in compact form, so evaluate over all distMat and convert to matrix form
      
      if(onlyUpper)
        return(compactToMat(exp(-distMat), diagVal=1))
      else
        #if onlyUpper==FALSE, fill in lower triangle of covariance matrix as well
        return(compactToMat(exp(-distMat), diagVal=1, lower.tri=TRUE))
      
    }
    else {
      #distMat is an actual matrix
      
      #only evaluate upper triangle of covariance matrix if possible
      if(onlyUpper && nrow(distMat) == ncol(distMat))
        return(ExponentialUpper(distMat))
      else
        return(exp(-distMat))
    }
    
  }
  #
  # multiply cross covariance matrix by C
  # in this case implemented in C
  #
  else if (!is.na(C[1])) {
    return(.Call("multebC", 
                 nd = as.integer(d),
                 x1 = as.double(x1), 
                 n1 = as.integer(n1),
                 x2 = as.double(x2),
                 n2 = as.integer(n2), 
                 par = as.double(p),
                 c = as.double(C),
                 work = as.double(rep(0, n2)) , PACKAGE="fields")
           )
  }
  #
  # return marginal variance ( 1.0 in this case)
  else if (marginal) {
    return(rep(1, nrow(x1)))
  }
  
  #not possible to reach this point
}
