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
##END HEADER


# test of sreg and related functions

suppressMessages(library(fields))
#options(echo=FALSE)

test.for.zero.flag<- 1



# simple covariance function for implementation
exp_cov <- function(dist){
  sigma2<-1 
  covariance <- sigma2* exp(-dist / 10) # 10 is arbitrary 
  return(covariance)
}

# -----------------------------
# Define grid and observations
# -----------------------------
 
m<- 10
n<- 11
nx<- m
ny<- n
M<- 15
dx<- 1
dy<- 1
sigma2<-2.0 
NNSize<-3 


# first a case where all obs in same grid box.
# addition of "dx" also tests that this works when grid is not just integers
# set dx=1 for the most basic case
dx<- .5 
s0<- rbind( 
            c(5.1,6.2),
            c(5.1,6.5),
            c( 5.85,6.45)
            )
 s0<- s0*dx
 
test0<-  offGridWeights( s0, list( x= (1:m)*dx, y=(1:n)*dx),
                                              aRange=10*dx, sigma2=sigma2, 
                                              Covariance="Exponential", 
                                              NNSize=2, 
                                              debug=TRUE) 
# explicit nearest neighbors in this case
sTmp<- cbind( rep(4:7,4), rep(5:8,each=4) )
sGrid<- sTmp*dx

# check that same grid being used by function
test.for.zero(sGrid, cbind(test0$gridX[,1], test0$gridY[,1]) )

S21  <- 2.0* exp( -rdist( s0, sGrid)/(10*dx))
S11  <- 2.0* exp( -rdist( sGrid , sGrid)/(10*dx) )
S22 <-  2.0* exp( -rdist( s0, s0)/(10*dx))
# local weights applied for prediction 
Btest<- S21%*% solve( S11)
# find indices for neigborhood
sIndex<- sTmp[,1] + (sTmp[,2]-1)*m
# Kriging weights
Bfull<-  spam2full(test0$B[,sIndex])
test.for.zero( Bfull, Btest)
# standard error matrix 
# note that transpsoe also taken so SEtest%*%t( SEtest) = cov matrix
SEtest<-  t(chol(S22 - S21%*% solve( S11)%*%t(S21) ))
SEfull <- spam2full(test0$SE)
test.for.zero( SEfull, SEtest)


# now test several observation locations

 dx<- .45
 s<- rbind( 
   c(5.1,6.2),
   c(7.1,7.2),
   c(5.1,6.5),
   c(8.5,4.4),
   c( 5.85,6.45),
   c(7.3,7.4)
 )
 s<- s * dx
# Note s0 from above is  s[c(1,3,5),]
ind1<- c(1,3,5)
 
 sTmp<- cbind( rep(4:7,4), rep(5:8,each=4) )
 sGrid<- sTmp*dx
 sIndex<- sTmp[,1] + (sTmp[,2]-1)*m
 
 S21<- 2.0* exp( -rdist( s[ind1,], sGrid)/(10*dx) )
 S11<- 2.0* exp( -rdist( sGrid , sGrid)/(10*dx) )
 S22<- 2.0* exp( -rdist( s[ind1,], s[ind1,])/(10*dx) )
 
 sparseObj<-  offGridWeights( s, list( x= (1:m)*dx, y=(1:n)*dx),
                                 aRange=(10*dx), sigma2=sigma2, 
                                 Covariance="Exponential", 
                                 NNSize=2, 
                                 debug=TRUE)
 
test.for.zero( sparseObj$Sigma21Star[ind1,], S21 )
test.for.zero( sparseObj$Sigma11Inv, solve(S11) )

Btest<- S21%*% solve( S11)
look2<- spam2full( sparseObj$B)
test.for.zero( Btest,look2[ind1, sIndex] )


SEfull<- spam2full( sparseObj$SE)
SE2full<- (SEfull)%*%t(SEfull)
test.for.zero(diag( SE2full), sparseObj$predictionVariance )

SEtest<- t(chol(S22 - S21%*%solve( S11)%*%t( S21) ))
test.for.zero(SEtest, SEfull[ind1, ind1] )

# check that debug FALSE also works

sparseObj1<-  offGridWeights( s, list( x= (1:m)*dx, y=(1:n)*dx),
                                aRange=(10*dx), sigma2=sigma2, 
                                Covariance="Exponential", 
                                NNSize=2, 
                                debug=FALSE)

test.for.zero( sparseObj$B, sparseObj1$B)
test.for.zero( sparseObj$SE, sparseObj1$SE)

cat("all done with off grid weight tests part 2", fill=TRUE)

