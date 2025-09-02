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
options(echo=FALSE)
test.for.zero.flag<- 1

data(ozone2)

x<- ozone2$lon.lat
y<- ozone2$y[16,]

temp<- Rad.cov( x,x, p=2)
temp2<- RadialBasis( rdist( x,x), M=2, dimension=2)

temp3<-  rdist( x,x)
temp3 <- ifelse( abs(temp3) < 1e-14, 0,log( temp3)*(temp3^2) )
temp3<- radbas.constant( 2,2)*temp3

test.for.zero( temp, temp2, tag="Tps radial basis function 2d")
test.for.zero( temp, temp3, tag="Tps radial basis function 2d")
test.for.zero( temp2,temp3, tag="Tps radial basis function 2d")


set.seed( 123)
xtemp<- matrix( runif( 40*3), ncol=3) 
temp<- Rad.cov( xtemp,xtemp, p= 2*4-3)
temp2<- RadialBasis( rdist( xtemp,xtemp), M=4, dimension=3)

temp3<-  rdist( xtemp,xtemp)
temp3 <- ifelse( abs(temp3) < 1e-14, 0, temp3^(2*4 -3) )
temp3<- radbas.constant( 4,3)*temp3

test.for.zero( temp, temp2, tag="Tps radial basis function 3d")
test.for.zero( temp, temp3, tag="Tps radial basis function 3d")
test.for.zero( temp2,temp3, tag="Tps radial basis function 3d")

#### testing multiplication of a vector
#### mainly to make the FORTRAN has been written correctly
#### after replacing the ddot call with an explicit  do loop
set.seed( 123)
C<- matrix( rnorm( 10*5),10,5 )
x<- matrix( runif( 10*2), 10,2)
temp3<- rdist( x,x)
K<-  ifelse( abs(temp3) < 1e-14, 0,log( temp3)*(temp3^2) )
K<- K * radbas.constant( 2,2)
test.for.zero( Rad.cov( x,x,m=2, C=C) , K%*%C, tol=1e-10)

set.seed( 123)
C<- matrix( rnorm( 10*5),10,5 )
x<- matrix( runif( 10*3), 10,3)
temp3<- rdist( x,x)
K<-  ifelse( abs(temp3) < 1e-14, 0,(temp3^(2*4-3)) )
K<- K * radbas.constant( 4,3)
test.for.zero( Rad.cov( x,x,m=4, C=C) , K%*%C,tol=1e-10)


#####  testing derivative formula

set.seed( 123)
C<- matrix( rnorm( 10*1),10,1 )
x<- matrix( runif( 10*2), 10,2)
temp0<-  Rad.cov( x,x, p=4, derivative=1, C=C)

eps<- 1e-6
temp1<- (
           Rad.cov( cbind(x[,1]+eps, x[,2]),x, p=4, derivative=0, C=C) 
         - Rad.cov( cbind(x[,1]-eps, x[,2]),x, p=4, derivative=0, C=C) )/ (2*eps)
temp2<- (
           Rad.cov( cbind(x[,1], x[,2]+eps),x, p=4, derivative=0, C=C) 
         - Rad.cov( cbind(x[,1], x[,2]-eps),x , p=4,derivative=0,C=C) )/ (2*eps)

test.for.zero( temp0[,1], temp1, tag=" der of Rad.cov", tol=1e-6)
test.for.zero( temp0[,2], temp2, tag=" der of Rad.cov", tol=1e-6)



# comparing  Rad.cov used by Tps with simpler function called 
# by stationary.cov
set.seed( 222)
x<- matrix( runif( 10*2), 10,2)
C<- matrix( rnorm( 10*3),10,3 ) 
temp<- Rad.cov( x,x, p=2, C=C)
temp2<- RadialBasis( rdist( x,x), M=2, dimension=2)%*%C
test.for.zero( temp, temp2)

#### Basic matrix form for Tps as sanity check
data("ozone2")
s<- ozone2$lon.lat
y<- ozone2$y[16,]

good<- !is.na( y)
s<- s[good,]
y<- y[good]
data(ozone2)


obj<-Tps( s,y, scale.type="unscaled", with.constant=FALSE)

# now work out the matrix expressions explicitly
lam.test<- obj$lambda
N<-length(y)

Tmatrix<- cbind( rep( 1,N), s)
D<- rdist( s,s)
R<- ifelse( D==0, 0, D**2 * log(D))
A<- rbind(
          cbind( R+diag(lam.test,N), Tmatrix),
          cbind( t(Tmatrix), matrix(0,3,3)))

 hold<-solve( A, c( y, rep(0,3)))
 c.coef<- hold[1:N]
 d.coef<- hold[ (1:3)+N]
 zhat<-  R%*%c.coef + Tmatrix%*% d.coef
  test.for.zero( zhat, obj$fitted.values, tag="Tps 2-d m=2 sanity check")
# out of sample prediction
snew<- rbind( c( -87,41),
              c( -81,44)
              )
T1<- cbind(  1, snew)
D<- rdist( snew,s)
R1<- ifelse( D==0, 0, D**2 * log(D))
z1<-  R1%*%c.coef + T1%*% d.coef
  test.for.zero( z1, predict( obj, x=snew), tag="Tps 2-d m=2 sanity predict")

#### test Tps verses Krig note scaling must be the same
   out<- Tps( s,y)
   out2<- Krig( s,y, Covariance="RadialBasis", 
           M=2, dimension=2, scale.type="range", method="GCV")
   test.for.zero( predict(out), predict(out2), tag="Tps vs.  Krig w/ GCV")

# test for fixed lambda
   test.for.zero( 
   predict(out,lambda=.1), predict(out2, lambda=.1),
   tag="Tps vs. radial basis w Krig")

#### testing derivative using predict function 
   set.seed( 233)
   x<- matrix( (rnorm( 1000)*2 -1), ncol=2)
   y<- (x[,1]**2 + 2*x[,1]*x[,2] -  x[,2]**2)/2

   out<- Tps( x, y, scale.type="unscaled")

   xg<- make.surface.grid( list(x=seq(-.7,.7,,10),  y=seq(-.7,.7,,10)) )
   test<- cbind( xg[,1] + xg[,2], xg[,1] - xg[,2])
#   test<- xg
   look<- predictDerivative.Krig( out, x= xg) 
   test.for.zero( look[,1], test[,1], tol=1e-3)
   test.for.zero( look[,2], test[,2], tol=1e-3)

############################################################
### testing Tps version using spatialProcess and Tps.cov
############################################################
   
   set.seed(222)
   n<- 50
   x1<- cbind( runif(n), runif(n))*100
   x2<-  cbind( runif(5), runif(5))
   #x2<- x1
   cardinalX<- cbind( runif(3), runif(3))
   m<- 2
   
   # simple check of marginal variances 
   look<- Tps.cov( x1,x1,cardinalX, m=m)
   look2<- Tps.cov( x1,cardinalX=cardinalX, m=m, marginal=TRUE)
   test.for.zero(diag(look), look2, tag="Tps.cov marginal" )
   
   
## comparing with the Tps function   
   data("ozone2")
   s<- ozone2$lon.lat
   y<- ozone2$y[16,]
   
   good<- !is.na( y)
   s<- s[good,]
   y<- y[good]
##### Tps used as benchmark
   
   out0<- Tps( s,y, scale.type ="unscaled", method="REML")
   lambdaHat<- out0$lambda.est[6,1]
   fHat<- predict( out0)
   
   cardinalX<- s[1:3,]
   out2<- mKrig( s,y, cov.function="Tps.cov",
                 cov.args= list( cardinalX=cardinalX,
                                 aRange=NA),
                 m=2,lambda=lambdaHat
   )
   
   # should be invariant to cardinal points and not need aRange =NA
   # defaults supplied by spatialProcessSetDefaults when detecting Tps.cov
   out3<- spatialProcess( s, y, cov.function="Tps.cov",
                             mKrig.args = list(m=2,  NtrA = 200),
                            lambda=lambdaHat )
   
   out4<- spatialProcess( s, y, cov.function="Tps.cov",
                          cov.args= list( 
                            aRange=NA),
                          #REML=TRUE,
                          verbose=FALSE, gridN=50
   )
   
   test.for.zero(fHat, predict( out2, s), tag="Tps vs spatialProcess" )
   # other parameters and likelihood
   test.for.zero(out0$tauHat.MLE, 
                 out3$summary["tau"],tag="Tps vs spatialProcess tau" )
   test.for.zero(fHat, predict( out3, s),
                 tag="Tps vs spatialProcess default" )
   # eff.df exact for spatialProcess because nTrA is larger than 
   # sample size.
   test.for.zero( out0$eff.df, out3$summary["eff.df"],
                  tag="Eff.df Tps vs spatialProcess")
   # look at prediction standard error computation.
   SE0<- predictSE( out0, s)
   SE3<- predictSE( out3,s)
   test.for.zero(SE0, SE3, tag="Tps vs spatialProcess SE")
# compare REML for Tps and spatialProcess  
   # test.for.zero(out0$lambda.est["REML", "-lnLike Prof"], 
   #               out4$summary["lnProfileREML.FULL"],
   #               tag="Tps vs spatialProcess  REML log like" )
   # test.for.zero(out0$lambda.est["REML", "lambda"], 
   #               out4$summary["lambda"],
   #               tag="Tps vs spatialProcess w REML lambda Hat" )
   # 
   gridList<- fields.x.to.grid( s, nx=20, ny=20)
   
   sGrid<- make.surface.grid( gridList)
   
   fHatGrid0<- predict( out0, sGrid)
   fHatGrid3<- predict( out3,sGrid)
   test.for.zero(fHatGrid0,fHatGrid3,
                 tag="Tps vs spatialProcess predictions on a grid")
   
   fHatGrid0SE<- predict( out0, sGrid)
   fHatGrid3SE<- predict( out3,sGrid)
   test.for.zero(fHatGrid0SE,fHatGrid3SE,
                 tag="Tps vs spatialProcess SE predictions on a grid")
   
options( echo=TRUE)
cat("all done testing Tps and spatialProcess with Tps.cov", fill=TRUE)
# 
# lambda0<- out0$lambda.est["REML", "lambda"]
# 
# Krig.flplike(out0, lambda0)

