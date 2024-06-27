setwd("~/Dropbox/ComoSchool/Data/GHCNStuff")

library( fields)

 
 load("NAPrcpSubset.rda" )
 
 # if any of JJA missing then report as missing
 NAPrcpJJA<- apply(NAPrcpSubset[,,6:8], c(1,2), sum)
 s<- NAPrcpLocsSubset
 
# NorthAmericanRainfall
# 100th meridian story
M<- nrow(NAPrcpInfoSubset )
coef<-matrix( NA, M,2)
coefSE<- matrix( NA, M,2)

x<- 1971:2024
x<- x - 1997 # centered so intercept estimates at year=1997

for(  k in 1:M){
  cat(k, " ")
  y<- NAPrcpJJA[k,]
  ind<- !is.na(y)
  yTmp<- y[ind]
  xTmp<- x[ind]
  fit<- lm( yTmp~xTmp)
  look<- summary( fit)$coefficients
  coef[k,]<- look[,1]
  coefSE[k,]<- look[,2]
}


NorthAmericanRainfall2<- 
  list( longitude= NAPrcpLocsSubset[,1],
        latitude= NAPrcpLocsSubset[,2],
        precip= coef[,1],
        precipSE= coefSE[,1],
        elevation=  NAPrcpInfoSubset$elev,
        trend= coef[,2],
        trendSE= coefSE[,2]
  )
library(mapproj)
xStereo<- mapproject( NorthAmericanRainfall2$lon,NorthAmericanRainfall2$lat, projection="stereographic")
NorthAmericanRainfall2$x.s<- cbind( xStereo$x, xStereo$y)
NorthAmericanRainfall2$projection<- .Last.projection

save(NorthAmericanRainfall2,
     file="NorthAmericanRainfall2.rda")