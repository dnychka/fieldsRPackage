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

# tests of SEs from the Omega matrix 
# using Monte Carlo

suppressMessages(library(fields))
options( echo=FALSE)

test.for.zero.flag<- TRUE
X<- ChicagoO3$x
n<- nrow( X)
tau2<- .1
sigma<- 2.0
processCov<-  sigma*Exp.cov( X,X,aRange=50)

cholCov<- chol( processCov + diag( tau2 , n ) )

nreps<- 1e5

set.seed( 111)
# extra random covariate
Z<- matrix( rnorm(n),n)
# the stochastic part 
E<- matrix( rnorm(n*nreps),n,nreps)
# NOTE: all fixed effects set to zero
Y<- t( cholCov)%*%E 
out<- mKrig( X,Y, aRange=50, Z=Z,
             collapseFixedEffect = FALSE,
              lambda=.1/2.0)
testCov<- var( t(out$beta) )

test.for.zero(testCov, sigma*out$Omega, tol=.05  )












