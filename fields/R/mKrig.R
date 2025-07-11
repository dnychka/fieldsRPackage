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
mKrig <- function(x, y, weights=rep(1, nrow(x)), Z = NULL, ZCommon=NULL,
                  cov.function="stationary.cov", 
                  cov.args = NULL, lambda = NA, m = 2, 
                  chol.args = NULL, find.trA = TRUE, NtrA = 20, 
                  iseed = NA, na.rm=FALSE, 
                  collapseFixedEffect = TRUE, 
                  tau=NA, sigma2=NA, verbose=FALSE, simpleKriging=FALSE, 
                  ...) {
  # pull extra covariance arguments from ...  and overwrite
  # any arguments already named in cov.args
  ind<- match( names( cov.args), names(list(...) ) )
  cov.args = c(cov.args[is.na(ind)], list(...))
  #
  #If cov.args$find.trA is true, set onlyUpper to FALSE (onlyUpper doesn't
  #play nice with predict.mKrig, called by mKrig.trace)
  #
  # next function also omits NAs from x,y,weights, and Z  if na.rm=TRUE.
  object<- mKrigCheckXY( x, y, weights, Z, ZCommon,na.rm = na.rm)
  # as the computation progresses additional components are 
  # added to the object list and by the end this is the returned 
  # object of class mKrig.
  
  if(find.trA == TRUE && supportsArg(cov.function, "onlyUpper"))
    cov.args$onlyUpper= FALSE
  if(find.trA == TRUE && supportsArg(cov.function, "distMat"))
    cov.args$distMat= NA
 
  if( !is.na(tau)|!is.na(sigma2)){
    fixedParameters<- TRUE
# work through the 3 cases for sigma2 and tau  
# note that for 2 of these also need lambda
    if( !is.na(tau)&!is.na(sigma2)){
    lambda<- tau^2/sigma2}
    if( is.na(tau)){
      tau <- sqrt( lambda*sigma2)
    }
    if( is.na(sigma2)){
      sigma2 <- tau^2/lambda
    }
  }
  else{
    fixedParameters<- FALSE
  }
  
  object$fixedParameters<- fixedParameters
  
  # check for duplicate x's.
  # stop if there are any
  if (any(duplicated(cat.matrix(x)))) {
    stop(" locations are not unique see help(mKrig) to
         collapse to unique locations and weighted obs ")
  }
  # create fixed part of model as m-1 order polynomial
  # NOTE: if m==0 then fields.mkpoly returns a NULL to 
  # indicate no polynomial part.
  if( simpleKriging){
    m<- 0
  }
  Tmatrix <- cbind(fields.mkpoly(object$x, m), object$Z)
  # set some dimensions
    np <- nrow(object$x)
  if( is.null(Tmatrix) ){
    nt<- 0
    }
  else{
    nt<- ncol(Tmatrix) 
  }
  if( is.null(object$Z)){
    nZ<- 0
  }
  else{
    nZ<- ncol(object$Z)
  }
  if( simpleKriging & (nZ>0)){
    stop("simpleKriging is TRUE but some covariates
         have been passed to mKrig")
  }
  ind.drift <- c(rep(TRUE, (nt - nZ)), rep(FALSE, nZ)) 
  # as a place holder for reduced rank Kriging, distinguish between
  # observations locations and  the locations to evaluate covariance.
  # (this is will also allow predict.mKrig to handle a Krig object)
  object$knots <- object$x
  # covariance matrix at observation locations
  # NOTE: if cov.function is a sparse constuct then Mc will be sparse.
  # see e.g. wendland.cov
  covArgsFull<- c(cov.args=list(cov.args), 
                  list(x1 = object$knots,
                       x2 = object$knots)
  )
                    
  if( verbose){
    cat("********************************",  fill=TRUE)
    cat("**** Main call to chol in mKrig", fill=TRUE)
    cat("***cov.function: ",  fill=TRUE)
    print(cov.function)
    cat("***cov.args:  ",  fill=TRUE)
    cat( names( cov.args), fill=TRUE, sep=",")
  }
  
  Mc <- do.call(cov.function, c(cov.args, list(x1 = object$knots, 
                                               x2 = object$knots)))

  #
  # decide how to handle the pivoting.
  # one wants to do pivoting if the matrix is sparse.
  # if Mc is not a matrix assume that it is in sparse format.
  #
  sparse.flag <- !is.matrix(Mc)
  #
  # set arguments that are passed to cholesky
  #
  if (is.null(chol.args)) {
    chol.args <- list(pivot = sparse.flag)
  }
  else {
    chol.args <- chol.args
  }
  # quantify sparsity of Mc for the mKrig object
  if( verbose){
    cat("**** In mKrig", fill=TRUE)
    cat("**** sparse flag",  sparse.flag, fill=TRUE)
  }
  nzero <- ifelse(sparse.flag, length(Mc@entries), np^2)
  # add diagonal matrix that is the observation error Variance
  # NOTE: diag must be a overloaded function to handle sparse format.
  if (lambda != 0) {
    if(! sparse.flag)
        invisible(.Call("addToDiagC", Mc, as.double(lambda/object$weights), nrow(Mc), PACKAGE="fields")
                  )
    else
      diag(Mc) = diag(Mc) + lambda/object$weights
  }
  #  MARK LINE Mc
  # At this point Mc is proportional to the covariance matrix of the
  # observation vector, y.
  #
  # cholesky decoposition of Mc
  # do.call used to supply other arguments to the function
  # especially for sparse applications.
  # If chol.args is NULL then this is the same as
  #              Mc<-chol(Mc), chol.args))
  Mc <- do.call("chol", c(list(x = Mc), chol.args))

  lnDetCov <- 2 * sum(log(diag(Mc)))
  
  #
  # start linear algebra to find estimates and likelihood
  # Note that all these expressions make sense if y is a matrix
  # of several data sets and one is solving for the coefficients
  # of all of these at once. In this case beta and c.coef are matrices
  #
  if( !is.null(Tmatrix) | !is.null(ZCommon) ){
    # Efficent way to multply inverse of Mc times the Tmatrix  
    TStar<- forwardsolve(Mc, x = Tmatrix, k=ncol(Mc), transpose = TRUE, upper.tri = TRUE)
    qr.VT <- qr(TStar)
    # GLS covariance matrix for fixed part.
    Rinv <- solve(qr.R(qr.VT))
    Omega <- Rinv %*% t(Rinv)
    lnDetOmega<- sum( log( diag(Rinv)^2 ) ) 
    # now do generalized least squares for beta
    YStar<- forwardsolve(Mc, transpose = TRUE, 
                         object$y, upper.tri = TRUE)
  }
  ##########################################
  ### only the T and Z covariate fixed parts
  ##########################################
  if( !is.null(Tmatrix) & is.null(ZCommon) ){
    beta <- as.matrix(qr.coef(qr.VT, YStar))
    if (collapseFixedEffect) {
      # use a common estimate of fixed effects across all replicates      
      betaMeans <- rowMeans(beta)
      beta <- matrix(betaMeans, ncol = ncol(beta), 
                       nrow = nrow(beta))
    }
#    
#  Omega from above  is  solve(t(Tmatrix)%*%solve( Sigma)%*%Tmatrix)
#   where Sigma = cov.function( x,x) + lambda/object$weights    
#   proportional to fixed effects covariance matrix.
#   for the GLS estimates of
#   the fixed linear part of the model. 
#    
#  SEdcoef = diag( Omega) * sigma2.MLE.FULL
# 
# if fixed effects are pooled across replicate fields then
# adjust the Omega matrix to reflect a mean estimate.
    if (collapseFixedEffect) {
      Omega <- Omega/ ncol(beta)
    }
# set ZCommon  parameters to NA   
  gamma<- NA
  OmegaZCommon<- NULL
  
#   GLS residual now used to evaluate likelihood
  resid<-  object$y - Tmatrix %*% beta
  
  }
  if( !is.null(ZCommon) ){
    if( is.null(T)){
      stop("need a fixed part matrix  (m>0, T and/or Z) to add ZCommon")
    }
    # check dimensions
     n<- nrow(object$y)
     M<- ncol( object$y)
     if( M ==1){
       stop("No replications just use the Z argument!")
     }
     N<- n*M
     if( nrow( ZCommon)!= N){
       stop("dimension of ZCommon should be nObs X nReps")
     }
     
     ZCStar<- matrix( NA, N, ncol(ZCommon))
     
     for( k in 1:ncol(ZCommon) ) {
      ZCtmp<- matrix(ZCommon[,k],n,M )
      temp<- forwardsolve(Mc, 
                          x = ZCtmp,
                          k=ncol(Mc),
                          transpose = TRUE,
                          upper.tri = TRUE)
       ZCStar[,k]<- c(temp)
     }
     testZ<- matrix( NA, N, ncol(ZCommon))
  # Project ZCommon onto subspace orthogonal to  TStar  
     for( k in 1:ncol(ZCommon) ) {
       temp<-  qr.resid(qr.VT, matrix( ZCStar[,k],n,M) )
       testZ[,k]<- c(temp)
     }
     testY<- c( qr.resid(qr.VT,  YStar) )
     ##########################################
     # ifelse block on collapseFixedEffects
     ##########################################
     if (!collapseFixedEffect) {
     
     gamma<- lsfit(  testZ,testY, intercept=FALSE)$coefficients
     OmegaZCommon<-  solve( t( testZ)%*%testZ )
     
     # find correct beta having adjusted by gamma   
     tmpResid<- YStar -  matrix(ZCStar%*%gamma,n,M)
     beta<- qr.coef(qr.VT, tmpResid )
     }
     else{
       # collapse beta fit -- common fixed effect parameters in T
       # but need to adjust for ZCommon that might not be the same at 
       # each realization.
       UStar<- cbind( rep(1,M)%x%TStar , ZCStar)
       allPar<- lsfit( UStar,c(YStar), intercept=FALSE)$coefficients
       # extract beta and gamma 
       nBeta<- ncol(TStar)
       nZC<- ncol(ZCommon)
       betaCommon<- allPar[1: nBeta]
       gamma<- allPar[(1: nZC)+ nBeta]
       # repeat beta for all realizations to have consitent format with 
       # collapseFixedEffects =FALSE
       beta <- matrix(betaCommon, ncol = M, 
                      nrow = nBeta)
       # correct Omega matrices
       Omega<- solve( t(UStar)%*%UStar)
       OmegaZCommon<- Omega[(1: nZC)+ nBeta,(1: nZC)+ nBeta]
     }
     #   GLS residual now used to evaluate likelihood   
     resid<- object$y - Tmatrix%*%beta - matrix(ZCommon%*%gamma,n,M)
  }
  if( is.null(Tmatrix)){
# much is set to NULL because no fixed part of model    
    nt<- 0
    resid<- object$y
    Rinv<- NULL
    Omega<- NULL
    qr.VT<- NULL
    beta<- NULL
    lnDetOmega <- 0
    # set ZCommon  parameters to NULL  
    gamma<- NULL
    OmegaZCommon<- NULL
  }
  # and now find c.
  #  the coefficents for the spatial part.
  # if there is also a linear fixed part  resid are the residuals from the 
  # GLS regression.
  
  c.coef <- as.matrix(forwardsolve(Mc, transpose = TRUE,
                                   resid, upper.tri = TRUE))
  # save intermediate result this is   t(y- T beta)( M^{-1}) ( y- T beta)
  quad.form <- c(colSums(as.matrix(c.coef^2)))
  
  # compute full likelihood if 2 out three covariance parameters are given
  if(fixedParameters){
    lnLike<- lnProfileLike <- (-quad.form/(2*sigma2) - log(2 * pi) * (np/2) - (np/2) * 
                                  log(sigma2) - (1/2) * lnDetCov )
    lnLikeREML<- lnLike + (1/2) * lnDetOmega
    lnLike.FULL<- sum( lnLike)
    lnLikeREML.FULL<- sum(lnLikeREML)
  }
  else{
    lnLike<-NA
    lnLike.FULL<-NA
    lnLikeREML<-NA
    lnLikeREML.FULL<-NA
  }
  
  # find c coefficients
  c.coef <- as.matrix(backsolve(Mc, c.coef))
  # find the residuals directly from solution
  # to avoid a call to predict
  object$residuals <- lambda * c.coef/object$weights
  object$fitted.values <- object$y - object$residuals
  # MLE estimate of sigma and tau
  #    sigmahat <- c(colSums(as.matrix(c.coef * y)))/(np - nt)
  # NOTE if y is a matrix then each of these are vectors of parameters.
  sigma2.MLE <- (quad.form/np)
  #sigma2hat <- c(colSums(as.matrix(c.coef * object$y)))/np
  tau.MLE <- sqrt(lambda * sigma2.MLE)
  # the  log profile likehood with  sigma2.MLE  and  dhat substituted
  # leaving a profile for just lambda.
  # NOTE if y is a matrix then this is a vector of log profile
  # likelihood values.
  lnProfileLike <- (-np/2 - log(2 * pi) * (np/2) - (np/2) * 
                      log(sigma2.MLE) - (1/2) * lnDetCov)
  # see section 4.2 handbook of spatial statistics (Zimmerman Chapter)
  # for this amazing shortcut to get the REML version 
  lnProfileREML <-  lnProfileLike + (1/2) * lnDetOmega
  # following FULL means combine the estimates across all replicate fields 
  # mean for MLE is justified as it is assumed locations and weights the same across 
  # replicates. 
  sigma2.MLE.FULL <- mean(sigma2.MLE)
  tau.MLE.FULL <- sqrt(lambda * sigma2.MLE.FULL)
  # if y is a matrix then compute the combined likelihood
  # under the assumption that the columns of y are replicated
  # fields
  lnProfileLike.FULL <- sum((-np/2 - log(2 * pi) * (np/2) - 
                               (np/2) * log(sigma2.MLE.FULL) 
                               - (1/2) * lnDetCov)
                            )
  lnProfileREML.FULL <- sum((-np/2 - log(2 * pi) * (np/2) - 
                               (np/2) * log(sigma2.MLE.FULL) 
                             - (1/2) * lnDetCov
                             + (1/2) * lnDetOmega  )
                            )
  
  #
  # return coefficients and include lambda as a check because
  # results are meaningless for other values of lambda
  # returned list is an 'object' of class mKrig (micro Krig)
  # also save the matrix decompositions so coefficients can be
  # recalculated for new y values.  Make sure onlyUpper and 
  # distMat are unset for compatibility with mKrig S3 functions
  #
  
   if(!is.null(cov.args$onlyUpper))
     cov.args$onlyUpper = FALSE
   if(!is.null(cov.args$distMat))
     cov.args$distMat = NA
 # build return object except for effective  degrees of freedom computation
 # and the summary vector 
  replicateInfo = list(
    lnProfileLike = lnProfileLike,
    lnProfileREML =  lnProfileREML,
    lnLike= lnLike,
    lnLikeREML= lnLikeREML, 
    tau.MLE = tau.MLE, 
    sigma2.MLE = sigma2.MLE,
    quad.form = quad.form
  )
   object2 <-  
     list(
       beta = beta,
       gamma = gamma,
       c.coef = c.coef,
       nt = nt,
       np = np,
       lambda.fixed = lambda,
       cov.function.name = cov.function,
       args = cov.args,
       m = m,
       chol.args = chol.args,
       call = match.call(),
       nonzero.entries = nzero,
       replicateInfo = replicateInfo,
       lnLike.FULL = lnLike.FULL,
       lnLikeREML.FULL = lnLikeREML.FULL,
       lnDetCov = lnDetCov,
       lnDetOmega = lnDetOmega,
       Omega = Omega,
       lnDetOmega = lnDetOmega,
       OmegaZCommon = OmegaZCommon,
       qr.VT = qr.VT,
       Mc = Mc,
       Tmatrix = Tmatrix,
       ind.drift = ind.drift,
       nZ = nZ,
       fixedEffectsCov = Omega * sigma2.MLE.FULL,
       fixedEffectsCovCommon = OmegaZCommon * sigma2.MLE.FULL,
       collapseFixedEffect = collapseFixedEffect,
       simpleKriging=simpleKriging
     )
   
    object<- c( object, object2)
  #
  
  #
  # estimate effective degrees of freedom using Monte Carlo trace method.
  # this is optional because not needed for predictions and likelihood
  # but necessary for GCV
   
  if (find.trA) {
    object3 <- mKrig.trace(object, iseed, NtrA)
    object$eff.df <- object3$eff.df
    object$trA.info <- object3$trA.info
    object$GCV <- (sum(object$residuals^2)/np)/(1 - object3$eff.df/np)^2
    if (NtrA < np) {
      object$GCV.info <- (sum(object$residuals^2)/np)/(1 - object3$trA.info/np)^2
    }
    else {
      object$GCV.info <- NA
    }
  }
  else {
    object$eff.df <- NA
    object$trA.info <- NA
    object$GCV <- NA
  }
   
  ################### compile summary vector of parameters
  summaryPars<- rep(NA,10)
  names( summaryPars) <- c( "lnProfileLike.FULL","lnProfileREML.FULL",
                            "lnLike.FULL","lnREML.FULL",
                            "lambda" ,
                            "tau","sigma2","aRange","eff.df","GCV")
  summaryPars["lnProfileLike.FULL"]<- lnProfileLike.FULL
  summaryPars["lnProfileREML.FULL"]<- lnProfileREML.FULL
  summaryPars["lnLike.FULL"]<- lnLike.FULL
  summaryPars["lnREML.FULL"]<- lnLikeREML.FULL
  
  
  if( fixedParameters){
    summaryPars["tau"]  <- tau
    summaryPars["sigma2"]<- sigma2
    }
  else{  
    summaryPars["tau"]  <- tau.MLE.FULL
    summaryPars["sigma2"]<- sigma2.MLE.FULL
  }
  
  summaryPars["lambda"]<- lambda
  summaryPars["aRange"] <-ifelse( !is.null(cov.args$aRange), 
                                     cov.args$aRange, NA)
  summaryPars["eff.df"] <- object$eff.df
  summaryPars["GCV"] <- object$GCV
  object$summary<- summaryPars
  
  
  
  if( verbose){
    cat( "****summary of object from mKrig", fill=TRUE)
    print( object$summary)
  }
  ########################
  ### add in some depreciated components so that LatticeKrig 8.4
  ### passes its tests.
  ########################
  
  object$rho.MLE<- sigma2.MLE
  object$rho.MLE.FULL<- sigma2.MLE.FULL
  object$lnProfileLike<- lnProfileLike
  object$lnProfileLike.FULL<- lnProfileLike.FULL
  object$quad.form <- quad.form
  object$rhohat<- sigma2.MLE.FULL
  object$d<- beta
    
  class(object) <- "mKrig"
 
  return(object)
}
