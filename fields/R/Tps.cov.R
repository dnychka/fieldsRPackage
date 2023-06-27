
Tps.cov<-function (x1, x2 = NULL, cardinalX, m=2,
                   C = NA, aRange=NA,
                   marginal = FALSE
          ) 
{
  
  if (is.null(x2) ) {
    x2 <- x1
   }
  
  if( !is.na(aRange)){
    stop("Tps should  be passed an aRange parameter that is NA
         to indicate this not used in the covariance.")
  }
  
  if ( !marginal) {
   
    PCard<- fields.mkpoly(cardinalX, m=m)
    PCoef<- solve(PCard, diag(1, ncol(PCard)) )
    # Note by consturction  PCard%*%PCoef = I
     P1<- fields.mkpoly(x1, m=m)%*%PCoef
     P2<- fields.mkpoly(x2, m=m)%*%PCoef
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
        # bigE<- Rad.cov( x1,x2, m=m)
        # K1<-       Rad.cov( x1,      cardinalX, m=m, C=t(P2))
        # K2 <-   t( Rad.cov( x2,      cardinalX, m=m, C=t(P1)))
        # K3<- P1%*% Rad.cov(cardinalX,cardinalX, m=m, C=t(P2))
        # covMatrix<- bigE - K1 - K2 + K3 + P1%*%t( P2)
        # return(covMatrix%*%C)
      }
    }   
    else {
      return(rep(1, nrow(x1)))
    }
}



