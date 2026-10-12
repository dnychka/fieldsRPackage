#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2026 Colorado School of Mines
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

Tps.cov<-function (x1, x2 = NULL, cardinalX, mDegree=2,
                   C = NA, aRange=NA,
                   marginal = FALSE
          ) 
{
  m<- mDegree
  if (is.null(x2) ) {
    x2 <- x1
   }
  
  if( !is.na(aRange)){
    stop("Tps should  be passed an aRange parameter that is NA
         to indicate this not used in the covariance.")
  }
  
  PCard<- fields.mkpoly(cardinalX, m=m)
  PCoef<- solve(PCard, diag(1, ncol(PCard)) )
  # Note by construction  PCard%*%PCoef = I
  P1<- fields.mkpoly(x1, m=m)%*%PCoef
  P2<- fields.mkpoly(x2, m=m)%*%PCoef
  
  if ( !marginal) {
   
   
    # a computation that is less efficient but easier to read:
    # bigE <- Rad.cov( x1, x2, m=m)
    # K1<- Rad.cov( x1,cardinalX, m=m)%*%t(P2)
    # K2<- P1%*% t(Rad.cov( x2,cardinalX, m=m))
    # K3<- P1%*% Rad.cov(cardinalX,cardinalX, m=m)%*%t(P2)
    # covMatrix<- bigE - K1 - K2 + K3 + P1%*%t( P2)
    # and ... if C is passed  covMatrix%*%C
    # 
      if (is.na(C[1])) {
        bigE<- Rad.cov( x1,x2, m=m)
        K1<-       Rad.cov( x1,      cardinalX, m=m, C=t(P2))
        K2 <-   t( Rad.cov( x2,      cardinalX, m=m, C=t(P1)))
        K3<- P1%*% Rad.cov(cardinalX,cardinalX, m=m, C=t(P2))
        covMatrix<- bigE - K1 - K2 + K3 + P1%*%t( P2)
        return(covMatrix)
      }
      else{
        bigEC<- Rad.cov( x1,x2, m=m,C=C)
        K1C <-         Rad.cov( x1,      cardinalX, m=m, C=t(P2)%*%C)
        K2C <-  P1%*% Rad.cov(cardinalX,        x2, m=m, C=C) 
        K3C <-  P1%*% Rad.cov(cardinalX, cardinalX, m=m, C=t(P2)%*%C )
        K4C <- P1%*%(t( P2)%*%C)
        covMatrixC<- bigEC - K1C - K2C + K3C + K4C
        return( covMatrixC)
        # in more verbose form
        # bigE<- Rad.cov( x1,x2, m=m)
        # K1<-       Rad.cov( x1,      cardinalX, m=m, C=t(P2))
        # K2 <-   t( Rad.cov( x2,      cardinalX, m=m, C=t(P1)))
        # K3<- P1%*% Rad.cov(cardinalX,cardinalX, m=m, C=t(P2))
        # covMatrix<- bigE - K1 - K2 + K3 + P1%*%t( P2)
        # return(covMatrix%*%C)
      }
    }   
    else {
      # NOTE: diag ( bigE) is zero
      # K2 <- t( Rad.cov( x1,      cardinalX, m=m, C=t(P1)))
      # K1 <- t(K2)
      # K3 <- P1%*% Rad.cov(cardinalX,cardinalX, m=m, C=t(P1))
      # marginalVariance0<- diag(  - K2  - t( K1) + K3 + P1%*%t( P1) )
       DK2<- rowSums(Rad.cov( x1,      cardinalX, m=m)*P1)
       DK3<- colSums(Rad.cov(cardinalX,cardinalX, m=m, C=t(P1))*t(P1) )
       DK4<- rowSums(P1^2)
      marginalVariance<- (-2*DK2 + DK3 +DK4 )
      return( marginalVariance)
    }
}



