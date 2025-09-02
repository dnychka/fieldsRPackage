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
test.for.zero.flag<-1

set.seed( 245)

x<- runif(3)

coef<- runif( 5)
temp<- fields.evlpoly( x, coef)

temp2<- coef[1]

for(  k in (2:5) ){
temp2<- temp2 + coef[k]*x**(k-1)
}

test.for.zero( temp, temp2)


set.seed( 124)
x<-  matrix( runif(12), ncol=3)

fields.mkpoly(x, m=3)-> out

attr( out, "ptab")-> ptab

J<- nrow( ptab)

coef<- runif( J)
temp<- fields.evlpoly2( x, coef, ptab)

temp2<-out%*% coef

test.for.zero( temp,temp2)

fields.derivative.poly( x, m=3, coef)-> temp

fields.mkpoly( cbind( x[,1:2], x[,3]+1e-6), m=3)%*% coef-> temp2
fields.mkpoly( cbind( x[,1:2], x[,3]-1e-6), m=3)%*% coef-> temp3
temp2<- (temp2- temp3)/ 2e-6

test.for.zero( temp[,3], temp2)

cat("Done testing polynomial evaluation",fill=TRUE)

options( echo=FALSE)




