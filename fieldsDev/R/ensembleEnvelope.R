ensembleEnvelope<- function(object=NULL,ensemble= NULL, alpha.level =.95,
                            ...){
  if( is.null( ensemble)){
  ensemble<- sim.spatialProcess(object,
                           ...)
  }
 dimEnsemble<-  c( 1)
  muHat<- apply( ensemble,dimEnsemble,mean, na.rm = TRUE)
  SE<- apply( ensemble,dimEnsemble, sd, na.rm = TRUE)
  Zscore<- (ensemble - muHat)/SE
  ZMax<- apply( abs(Zscore), 1, max)
  bound<- quantile( ZMax, 1-alpha.level, na.rm = TRUE)
  return(
  cbind( muHat - SE*bound, muHat + SE*bound )
  )
}