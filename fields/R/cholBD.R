cholBD <- function(A) {
  ## obtain R such that R'R = A.
  # Where A is banded matrix contained in R.
  #
  #    The band storage scheme is illustrated by the following example, when
  #    N = 6, KD = 2, and UPLO = 'U':
  #  
  #    On entry:                       On exit:
  #  
  #        *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
  #        *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
  #       a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
  #  
  
  # CHARACTER          UPLO
  # INTEGER            INFO, KD, LDAB, N
  # *     ..
  # *     .. Array Arguments ..
  # DOUBLE PRECISION   AB( LDAB, * )
  #SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
  # UPLO set to 'U'
  n<- ncol( A)
  KD<- nrow( A)-1 # number of off diagonal bands
 out<- .Fortran( "bandchol", 
                 N = as.integer(n),
                 KD = as.integer(KD),
                  A = as.double(c(A)),
                 LDAB= as.integer(KD+1),
                 info=as.integer(0)
                 )
  if( out$info!=0){
    print( out$info)
    stop("error in call to bandChol")
  }
 return( out)
} ## bandchol

