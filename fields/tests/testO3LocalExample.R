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

data( "ozone2")

s<- ozone2$lon.lat
z<- ozone2$y[16,]
good<- !is.na(z)
s<- s[good,]
z<- z[good]

obj0<- spatialProcess(s,z, smoothness = 1.5)

sGrid2<- list( x= seq( min(s[,1]), max(s[,1]),length.out=45),
               y= seq( min(s[,2]), max(s[,2]),length.out=50)
)


# sanity check as to accuracy
#
fHat<-  predictSurface( obj0, gridList=sGrid2, extrap=TRUE)
fHat1<- predictSurface( obj0, gridList=sGrid2, fast =TRUE, NNSize=20, extrap=TRUE,
                        verbose=TRUE)
test.for.zero( fHat$z, fHat1$z, tol=1e-7, tag="First a sanity check for fast predict!")

look<- mKrigFastPredictSetup(obj0,
                                  sGrid2, 
                                  NNSize=4,
                                  giveWarnings=TRUE,
                                  verbose=FALSE
)




set.seed(222)
out1<- simLocal.spatialProcess(obj0, sGrid2, 
                                   NNSize=8,
                                   M=5,
                                   extrap=TRUE,
                                  verbose=FALSE) 
set.seed(222)
out2<- simLocal.spatialProcess(obj0, sGrid2, M=5, 
                               fast=FALSE,
                               NNSizePredict=8,
                               NNSize=8,
                               extrap=TRUE,
                               verbose=FALSE) 
# check that the same synthetic data sets are generated 
test.for.zero( out1$hHat, out2$hHat, 
               tag="hHat fast/exact", tol=.005)
# 
# test.for.zero( out1$outData, out2$outData, 
#                tag="syn data on fast/exact 5 reps")
# test.for.zero( out1$outhData, out2$outhData, 
#                tag="syn h  on fast/exact 5 reps")

# compare sample SEs
SE1<- apply( out1$z, c( 1,2), sd)
SE2<- apply( out2$z, c( 1,2), sd)

test.for.zero( SE2, SE1,
               tag="SEs local, fast/local, exact", tol=.01)

test.for.zero( out1$z, out2$z, tol=.025, tag="z on fast/exact 5 reps")

RMSE<- sqrt(mean(c(out1$z- out2$z)^2))
RRMSE<- sqrt(mean(c(out1$z- out2$z)^2))/sqrt(mean(c(out1$z)^2))

test.for.zero( RMSE, 0, relative=FALSE, tol=.0002, tag="RMSE z")

## unequal NNSize
set.seed(222)
out1<- simLocal.spatialProcess(obj0, sGrid2, 
                               NNSize=8,
                               M=5,
                               extrap=TRUE) 
set.seed(222)
out3<- simLocal.spatialProcess(obj0, sGrid2, M=5, 
                               fast=TRUE,
                               NNSizePredict=10,
                               NNSize=8,
                               extrap=TRUE) 
test.for.zero( out1$hHat, out3$hHat, tol= .001,
               tag="syn h  on fast/exact 5 reps")  


## compare to exact method
set.seed(222)
out4<- simLocal.spatialProcess(obj0, sGrid2, M=200, 
                               fast=TRUE,
                               NNSizePredict=10,
                               NNSize=8,
                               extrap=TRUE)
SEApprox<- apply(out4$z, c(1,2), sd)

set.seed(222)
xp<- make.surface.grid(sGrid2)
SE0<- predictSurfaceSE(obj0, grid.list=sGrid2,
                       extrap=TRUE)
outMonte<- sim.spatialProcess(obj0, xp, M=1000) 
SEMonte<- apply(outMonte, c(1), sd)


# test.for.zero(SEMonte, SE0$z, 
#               tag="Monte Carlo SE vs  Formula" )

RRMSE<- sqrt(mean( ((SEMonte- c(SE0$z))/c(SE0$z)) ^2))
test.for.zero( RRMSE,0, tol=.03, relative=FALSE, 
               tag="RRMSE Monte vs. exact")



# test.for.zero(SE0$z, SEApprox, 
#               tag="Exact SE vs local/fast" )

RRMSE<- sqrt(
  mean(
    ((c(SE0$z)- c(SEApprox))/c(SE0$z) )^2
    )
)
test.for.zero( RRMSE,0, tol=.05, relative=FALSE,
               tag="RRMSE local/fast vs. exact")

cat("all done!", fill=TRUE)
#
# below are some qq plots that confirm the accuracy of the code
#

# n<- ncol( outMonte)
# Chi2<- (SEMonte^2/c(SE0$z)^2)
# qqnorm( (Chi2 -1)* sqrt(n/2) )
# abline(0,1,col="red")

# n<- 200
# Chi2<- (c(SEApprox)^2/c(SE0$z)^2)
# qqnorm( (Chi2 -1)* sqrt(n/2) )
# abline(0,1,col="red")
# 
# bplot.xy( log10(SEMonte), log10(c(SE0$z)))
# abline( 0,1,col="red")
# 
# bplot.xy( log10(SEMonte), log10(c(SEApprox)))
# abline( 0,1,col="red")

# set.panel(1,2)
# imagePlot( as.surface(sGrid2, SEExact))
# points( s, col="magenta", cex=.2)
# imagePlot(out4$x, out4$y,  SEApprox)
# points( s, col="magenta", cex=.2)
# 
# imagePlot( as.surface(sGrid2, outExact[,1]))
