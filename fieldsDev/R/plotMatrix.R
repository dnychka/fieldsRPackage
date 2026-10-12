#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2026 Colorado School of Mines
# 1500 Illinois St., Golden, CO 80401
# Contact: Douglas Nychka,  douglasnychka@gmail.com,
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
plotMatrix<- function( A, M=5,N=5,...){
  
  n<- nrow(A)
  m<- ncol( A)
  
  obj<-  t((A)[n:1,] )
  
  imagePlot( obj,  axes=F, xlab="columns", ylab="rows",
             ...)
   nG<- round(pretty( 1:n, n=N) )
  axis(side=2,  at= (nG-1)/(n-1), labels= format( nG) )
  mG<-  seq( 1,m, round(m/M) )
  mG<- round( pretty( 1:m, n=M))
  axis(side=1,  at= 1- (mG-1)/(m-1),
       labels= format( mG), gap.axis=0, las=2 )
  box()
}