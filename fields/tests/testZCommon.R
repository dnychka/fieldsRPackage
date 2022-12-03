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


suppressMessages(library(fields))
options( echo=FALSE)
test.for.zero.flag<- 1

set.seed(122)

n<-6
M<- 3
N<- n*M
ZCommon<- matrix( round(100*runif( N*3)), N,3)

# data locs
#s<- cbind( 1:n)
s<- matrix( runif( n*2), n,2)
T<-  cbind( 1, s)
Sigma<- exp( -rdist(s,s)/2) 
beta<- rep( 5, ncol(T) )
gamma<- 1:3  

set.seed(111)
f<- t(chol(Sigma))%*% matrix(rnorm( n*M),n,M)
E<- matrix( .05*rnorm(n*M), n,M)
# the data obs
y0<- matrix( rep(T%*%beta, M),  n,M) + f +  E
y<- matrix( rep( T%*%beta, M) + ZCommon%*%gamma, n,M) + f +  E

# NOTE for fixed covariance parameters spatailProcess is a wrapper for a call to mKrig


look0<- spatialProcess(s,y0, aRange=2, lambda=.05^2, smoothness=.5)
look0B<- mKrig(s,y0, aRange=2, lambda=.05^2,
               collapseFixedEffect = TRUE)
             
look1<- spatialProcess(s,y0, aRange=2, lambda=.05^2, smoothness=.5,
                    collapseFixedEffect = FALSE)
look1B<- mKrig(s,y0, aRange=2, lambda=.05^2,
               collapseFixedEffect = FALSE)

test.for.zero( look1B$beta, look1$beta)



look2<- spatialProcess(s,y, aRange=2, lambda=.05^2, smoothness=.5,
                      ZCommon = ZCommon)
look2B<- mKrig(s,y, aRange=2, lambda=.05^2,
               ZCommon = ZCommon )

test.for.zero( look2B$beta, look2$beta)
test.for.zero( look2B$gamma, look2$gamma)

look3<- spatialProcess(s,y, aRange=2, lambda=.05^2,smoothness=.5,
                     collapseFixedEffect = FALSE,
                          ZCommon = ZCommon)

look3B<- mKrig(s,y, aRange=2, lambda=.05^2,
               collapseFixedEffect = FALSE,
               ZCommon = ZCommon )

# hand computations 

Gamma<- solve( chol(Sigma + diag(.05^2,n)) )
YStar<- c(t( Gamma)%*%y)
YStar0<- c(t( Gamma)%*%y0)
bigSigma<- diag(1,M)%x%(Sigma+ diag(.05^2,n))
bigGamma<- solve( chol(bigSigma) )
bigTA<- rep(1,M)%x%T # collapse ==TRUE
bigT<- diag(1,M)%x%T # collapse ==FALSE

# collapseFixedEffects = TRUE  
bigTAStar<- t( bigGamma)%*%bigTA
coefTestA<- lm( YStar0 ~ bigTAStar -1)$coefficients
test.for.zero(  look0$beta[,1] , coefTestA )


# collapseFixedEffects = FALSE 
bigTA<- diag(1,M)%x%T
bigTAStar<- t( bigGamma)%*%bigTA

coefTestA<- lm( YStar0 ~ bigTAStar -1)$coefficients
test.for.zero(  look1$beta , coefTestA )




# collapseFixedEffects =TRUE ZCommon
bigT<- rep(1,M)%x%T
U<- cbind(bigT,ZCommon )

UStar<- t( bigGamma)%*%U
coefTest<- lm( YStar ~ UStar -1)$coefficients

test.for.zero( c( look2$beta[,1], look2$gamma), coefTest,
               tag= "spatialProcess ZCommon collapse" )


# collapseFixedEffects = FALSE ZCommon
bigT<- diag(1,M)%x%T
U<- cbind(bigT,ZCommon )
UStar<- t( bigGamma)%*%U
coefTest<- lm( YStar ~ UStar -1)$coefficients

test.for.zero( c( look3$beta, look3$gamma), coefTest, tag= "spatialProcess ZCommon" )
test.for.zero( c( look3B$beta, look3B$gamma), coefTest,tag= "mKrig ZCommon"  )

####################################################
# generate a large set of  2D realizations
####################################################
set.seed(122)
options(echo=FALSE)
n<-100
M<- 1000
N<- n*M
ZCommon<- matrix( round(100*runif( N*3)), N,3)

# data locs
#s<- cbind( 1:n)
s<- matrix( runif( n*2), n,2)
T<-  cbind( 1, s)
Sigma<- exp( -rdist(s,s)/2) 
beta<- rep( 5, ncol(T) )
gamma<- 1:3  

set.seed(111)
f<- t(chol(Sigma))%*% matrix(rnorm( n*M),n,M)
E<- matrix( .05*rnorm(n*M), n,M)
# the data obs
y0<- matrix( rep(T%*%beta, M),  n,M) + f +  E
y<- matrix( rep( T%*%beta, M) + ZCommon%*%gamma, n,M) + f +  E

look<- spatialProcess(s,y0, aRange=2, lambda=.05^2, smoothness=.5,
                       ZCommon = ZCommon)
test.for.zero( look$beta[,1], beta, tol=2e-2,tag= "big data beta")
test.for.zero( look$summary["tau"], .05, tol=1e-2, tag="big data tau")
test.for.zero( look$summary["sigma2"],1, tol=1e-2, tag="big data sigma2")

look<- spatialProcess(s,y, aRange=2, lambda=.05^2, smoothness=.5,
                       ZCommon = ZCommon)
test.for.zero( look$beta[,1], beta, tol=2e-2,
               tag= "big data beta ZCommon")
test.for.zero( look$gamma, gamma, tol=1e-3,
               tag= "big data gamma ZCommon")
test.for.zero( look$summary["tau"], .05, tol=1e-2,
               tag="big data tau ZCommon")
test.for.zero( look$summary["sigma2"],1, tol=1e-2,
               tag="big data sigma2 ZCommon")



options(echo=TRUE)
