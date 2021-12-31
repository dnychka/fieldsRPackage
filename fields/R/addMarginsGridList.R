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
addMarginsGridList<- function( xObs, gridList, np){
  xMin<- min(xObs[,1])
  xMax<- max(xObs[,1]) 
  yMin<- min(xObs[,2])
  yMax<- max(xObs[,2]) 
  
  gridXMin<- gridList$x[1]
  gridYMin<- gridList$y[1]
  
  dx<- gridList$x[2] - gridList$x[1]
  nx<- length( gridList$x)
 
  ind1<-  floor(    (xMin -  gridXMin)/dx)
  ind2<-  ceiling(  (xMax -  gridXMin)/dx)
  if( (ind1 < 0)|  (ind2 > nx) ) {
    stop("locations outside of x grid" )
  }
  #indRangeX<- c( min( ind1 - np +1, 0), max( ind2 + np  ,nx ))
  # add extra grid box to avoid roundoff errors 
  indRangeX<- c( min( ind1 - np , 0), max( ind2 + np +1 ,nx ))
  xGrid<- (indRangeX[1]:(indRangeX[2] - 1)) *dx + gridXMin

  dy<- gridList$y[2] - gridList$y[1]  
  ny<- length( gridList$y)
  ind1<-  floor(    (yMin -  gridYMin)/dy)
  ind2<-  ceiling(  (yMax -  gridYMin)/dy)
  if( (ind1 < 0)|  (ind2 > ny) ) {
    stop("locations outside of y grid" )
  }
  #indRangeY<- c( min( ind1 - np +1 , 0), max( ind2 + np  , ny))
  indRangeY<- c( min( ind1 - np  , 0), max( ind2 + np +1 , ny))
  yGrid<- (indRangeY[1]:(indRangeY[2] - 1)) *dy + gridYMin
  
  # new gridList contains the orginal grid and extends it so that
  # there are np grid points beyond the min and max ranges for the coordinates. 
  gridListNew<- list( x = xGrid,
                      y = yGrid)
  indX<-  match( round(gridList$x,   13),
                 round(gridListNew$x,13) )
 # print( indX)
  
  
  indY<-  match(round(gridList$y,   13),
                round(gridListNew$y,13) )
 # print( indY)
  return( 
    list(gridListNew = gridListNew, gridList = gridList,
          indX= indX,
          indY= indY,
          np=np )
  )
}

#
#  uncomment to test 
#
#  xObs<- cbind( c( 1.2,4.9 ),c( 2.2, 5.5) )
#  gridList<-  list( x= seq( 1, 5,, 9), y= seq(1,6,, 21 ) )
#   marginInfo<- addMarginsGridList( xObs, gridList, 1)
#   gridListNew<- marginInfo$gridListNew
# test.for.zero(gridList$x, gridListNew$x[ marginInfo$indX] )
# test.for.zero(gridList$y, gridListNew$y[ marginInfo$indY] )
# sum( gridListNew$x <  min(xObs[,1] ) ) 
# sum( gridListNew$x > max(xObs[,1] ) )
# sum( gridListNew$y < min(xObs[,2] ) ) 
# sum( gridListNew$y >  max(xObs[,2] ) )
# 
# xObs<- cbind( c( 1,5 ),c( 1,6) )
# gridList<-  list( x= seq( 1, 5,, 9), y= seq(1,6,, 21 ) )
# marginInfo<- addMarginsGridList( xObs, gridList, 1)
# gridListNew<- marginInfo$gridListNew
# test.for.zero(gridList$x, gridListNew$x[ marginInfo$indX] )
# test.for.zero(gridList$y, gridListNew$y[ marginInfo$indY] )
# sum( gridListNew$x <  min(xObs[,1] ) ) 
# sum( gridListNew$x > max(xObs[,1] ) )
# sum( gridListNew$y < min(xObs[,2] ) ) 
# sum( gridListNew$y >  max(xObs[,2] ) )
# 
# xObs<- cbind( c( 1,5 ),c( 1,6) )
# gridList<-  list( x= seq( 1, 5,, 9), y= seq(1,6,, 21 ) )
# marginInfo<- addMarginsGridList( xObs, gridList, 4)
# gridListNew<- marginInfo$gridListNew
# test.for.zero(gridList$x, gridListNew$x[ marginInfo$indX] )
# test.for.zero(gridList$y, gridListNew$y[ marginInfo$indY] )
# sum( gridListNew$x <  min(xObs[,1] ) ) 
# sum( gridListNew$x > max(xObs[,1] ) )
# sum( gridListNew$y < min(xObs[,2] ) ) 
# sum( gridListNew$y >  max(xObs[,2] ) )
