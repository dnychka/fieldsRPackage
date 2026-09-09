naturalSplineBasis <- function(sGrid,
                               sKnots,
                               degree = 3,
                               derivative = 0) {
  sKnots0<- c( rep( min(sKnots),degree),sort(sKnots),
               rep( max(sKnots),degree) )
  
  basis <- splineDesign(sKnots0, sGrid,
                        ord= degree+1,
                        derivs=derivative)
  return( basis )
  
}