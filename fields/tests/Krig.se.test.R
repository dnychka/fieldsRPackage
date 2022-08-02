#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2021 Colorado School of Mines
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


suppressMessages(library(fields))

# tests of predictSE
# against direct linear algebra 

#options( echo=FALSE)



x0<- expand.grid( c(-8,-4,0,20,30), c(10,8,4,0))


out<- Krig( ChicagoO3$x, ChicagoO3$y, cov.function = "Exp.cov", aRange=50)


# direct calculation
Krig.Amatrix( out, x=x0)-> A
test.for.zero( A%*%ChicagoO3$y, predict( out, x0),tag="Amatrix vs. predict")

Sigma<- out$sigmahat*Exp.cov( ChicagoO3$x, ChicagoO3$x, aRange=50)
S0<- out$sigmahat*c(Exp.cov( x0, x0, aRange=50))
S1<- out$sigmahat*Exp.cov( out$x, x0, aRange=50)

#yhat= Ay
#var( f0 - yhat)=    var( f0) - 2 cov( f0,yhat)+  cov( yhat)

look<- S0 - t(S1)%*% t(A) - A%*%S1 +  
       A%*% ( Sigma + diag(out$tauHat.MLE**2/out$weightsM))%*% t(A)
#
#compare to 
# diagonal elements


test2<- predictSE( out, x= x0) 
test.for.zero( sqrt(diag(  look)), test2,tag="Marginal predictSE")

out2<- Krig( ChicagoO3$x, ChicagoO3$y, cov.function = "Exp.cov", aRange=50,
            lambda=out$lambda)

test2<- predictSE( out2, x= x0) 
test.for.zero( sqrt(diag(  look)), test2,tag="Marginal predictSE fixed ")

test<- predictSE( out, x= x0, cov=TRUE)
test.for.zero( look, test,tag="Full covariance predictSE")


# simulation based.

set.seed( 333)

sim.Krig( out, x0,M=4e3)-> test
 # columns are the realizations rows are locations

var(t(test))-> look

predictSE( out, x=x0)-> test2
mean( diag( look)/ test2**2)-> look2
test.for.zero(look2, 1.0, tol=1.5e-2, tag="Marginal standard Cond. Sim.")

predictSE( out, x=x0, cov=TRUE)-> test2

# multiply simulated values by inverse square root of covariance
# to make them white

eigen( test2, symmetric=TRUE)-> hold
hold$vectors%*% diag( 1/sqrt( hold$values))%*% t( hold$vectors)-> hold
cor(t(test)%*% hold)-> hold2
# off diagonal elements of correlations -- expected values are zero. 

abs(hold2[ col(hold2)> row( hold2)])-> hold3

test.for.zero(   mean(hold3), 0, relative=FALSE, tol=.02,
          tag="Full covariance standard Cond. Sim.")


# test of A matrix
#
# first create and check a gridded test case. 


data( ozone2)
as.image(ozone2$y[16,], x= ozone2$lon.lat, ny=24, nx=20, 
          na.rm=TRUE)-> dtemp
#
# A useful disctrtized version of ozone2 data
 
x<- dtemp$xd
y<- dtemp$z[ dtemp$ind]
weights<- dtemp$weights[ dtemp$ind]

Krig( x, y, Covariance="Matern", 
   aRange=1.0, smoothness=1.0, weights=weights) -> out



  set.seed(234)
  ind0<- cbind( sample( 1:20, 5), sample( 1:24, 5))

  x0<- cbind( dtemp$x[ind0[,1]], dtemp$y[ind0[,2]]) 

# an  inline check plot(out$x, cex=2); points( x0, col="red", pch="+",cex=2)

# direct calculation as backup ( also checks weighted case)

Krig.Amatrix( out, x=x0)-> A
test.for.zero( A%*%out$yM, predict( out, x0),tag="Amatrix vs. predict")

Sigma<- out$sigmahat*stationary.cov( 
out$xM, out$xM, aRange=1.0,smoothness=1.0, Covariance="Matern")

S0<- out$sigmahat*stationary.cov( 
x0, x0, aRange=1.0,smoothness=1.0, Covariance="Matern")

S1<- out$sigmahat*stationary.cov(
out$xM, x0, aRange=1.0,smoothness=1.0, Covariance="Matern")



#yhat= Ay
#var( f0 - yhat)=    var( f0) - 2 cov( f0,yhat)+  cov( yhat)
 
look<- S0 - t(S1)%*% t(A) - A%*%S1 +
       A%*% ( Sigma + diag(out$tauHat.MLE**2/out$weightsM) )%*% t(A)

test<- predictSE( out, x0, cov=TRUE)

test.for.zero( c( look), c( test), tag="Weighted case and exact for ozone2 full 
cov", tol=1e-8)


cat("all done testing predictSE.Krig ", fill=TRUE)
options( echo=TRUE)
