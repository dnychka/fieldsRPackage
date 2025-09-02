#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2024 Colorado School of Mines
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

###  gridList<- list( x= seq( 0,1,.2), y= seq( -1,1,.25 ))
###  set.seed(112)
###   s<- cbind( runif( 20, -.5, 1.5), runif( 20,-2.5,2.4))
###  findGridBox( s, gridList)

findGridBox<-function(s, gridList)                       {
  #
  # function assumes the grid points are the centers and 
  # equally spaced in X and Y.
  
  m<- length( gridList$x)
  n<- length( gridList$y)
  
  dx<- gridList$x[2]- gridList$x[1]
  dy<- gridList$y[2]- gridList$y[1]
  # create new grid  of m+1 by n+1 assuming gridList are the grid box centers. 
  xGrid0<- c((gridList$x - dx/2),  gridList$x[m] + dx/2)
  yGrid0<- c((gridList$y - dy/2),  gridList$y[n] + dy/2)
  
  # lower left corner of grid box indices containing the locations s
  s0<-  cbind( 
    trunc( (s[,1]- xGrid0[1] )/dx) + 1 ,
    trunc( (s[,2]- yGrid0[1] )/dy) + 1
  ) 
  # set points of s  outside of grid to NA
  indX<- ( s0[,1]<1 | s0[,1]>m )
  indY<- ( s0[,2]<1 | s0[,2]>n )
  s0[indX,1]<- NA
  s0[indY,2]<- NA
  if( any(indX | indY )){
    warning( "Some points outside range of grid")
  }
  
  # print( s0)
  # 
  # # grid centers
  # sGrid<- make.surface.grid( gridList)
  # xr<- range( c( xGrid0,s[,1] ) )
  # yr<- range( c( yGrid0,s[,2] ) )
  # plot(s, pch=16, col="magenta", cex=1.2, xlim=xr, ylim=yr )
  # points( sGrid, col="grey")
  # xline( xGrid0); yline( yGrid0)
  # out<- (!indX) & (!indY )
  # print( s0[out,])
  # points( cbind(xGrid0[s0[out,1]], yGrid0[s0[out,2]]), pch="+", cex=2,col="blue")
  # 
  # index  of locations when 2D array is unrolled
  # note points outside grid are NAs
  s0Index<- as.integer( s0[,1] + (s0[,2]-1)*m)
  
  return(  
    list( index= s0Index, gridIndex= s0,  m=m, n=n ) 
    )
  
}
