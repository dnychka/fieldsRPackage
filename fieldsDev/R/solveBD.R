
  solveBD<- function(U,B){
#SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
#          UPLO='U'
#          N is INTEGER
#          The order of the matrix A.  N >= 0.
#          KD is INTEGER
#          The number of superdiagonals of the matrix A if UPLO = 'U',
#          NRHS  the number of columns
#          of the matrix B. 
#          AB is DOUBLE PRECISION array, dimension (LDAB,N)
#          The triangular factor U  from the Cholesky factorization
#          A = U**T*U  of the band matrix A, stored in the
#          first KD+1 rows of the array. 
#          LDAB is INTEGER
#          The leading dimension of the array AB.  LDAB >= KD+1.
# 
  #          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
  #          On entry, the right hand side matrix B.
  #          On exit, the solution matrix X.
  # 
  #          LDB is the leading dimension of the array B.  LDB >= max(1,N).
  #          INFO is INTEGER
  #          = 0:  successful exit
  #          < 0:  if INFO = -i, the i-th argument had an illegal value
  # 
  n<- ncol(U)
  LDAB<- nrow( U)
  LDB<- nrow(B)
  if( n !=LDB){
    print( c( n, LDB))
    stop("dimension of cholesky decomp and 
         row dimension of  B do not match.")
  }
  NRHS= ncol(B)
  KD<- nrow( U)-1
  out<- .Fortran( "bandsolve",
                  N =  as.integer(n), 
                  KD = as.integer( KD),
                  NRHS= as.integer(NRHS),
                  AB= as.double(U),
                  LDAB= as.integer(LDAB),
                  B=as.double(B), 
                  LDB= as.integer( LDB),
                  INFO= as.integer(0),
                  PACKAGE = "fields"
              )
  return(out)
  }