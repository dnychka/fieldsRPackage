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
addMarginsGridList<- function( xObs, gridList, NNSize){

  np<- NNSize
# ranges and sizes for the inoput grid and the obs  
  xMin<- min(xObs[,1])
  xMax<- max(xObs[,1]) 
  yMin<- min(xObs[,2])
  yMax<- max(xObs[,2]) 
  
  gridXMin<- gridList$x[1]
  gridYMin<- gridList$y[1]
  
  dx<- gridList$x[2] - gridList$x[1]
  nx<- length( gridList$x)
  dy<- gridList$y[2] - gridList$y[1]  
  ny<- length( gridList$y)
  
 #### x grid block
  ind1<-  floor(    (xMin -  gridXMin)/dx)
  ind2<-  ceiling(  (xMax -  gridXMin)/dx)
  
  if( (ind1 < 0)|  (ind2 > nx) ) {
    stop("locations outside of x grid" )
  }
  indRangeX<- c( min( ind1 - np , 0), max( ind2 + np +1 ,nx ))
  xGrid<- (indRangeX[1]:(indRangeX[2] - 1)) *dx + gridXMin

  ind1<-  floor(    (yMin -  gridYMin)/dy)
  ind2<-  ceiling(  (yMax -  gridYMin)/dy)
  if( (ind1 < 0)|  (ind2 > ny) ) {
    stop("locations outside of y grid" )
  }
  #### y grid block 
  indRangeY<- c( min( ind1 - np  , 0), max( ind2 + np +1 , ny))
  yGrid<- (indRangeY[1]:(indRangeY[2] - 1)) *dy + gridYMin
  
  # new gridList contains the orginal grid
  # up to round off based on  integer steps in dx and dy and extends
  # it so that
  # there are np grid points beyond the min and max ranges for 
  # the original ranges. 
  gridListNew<- list( x = xGrid,
                      y = yGrid)
  
  # these indices are used to locate the original grid as a subset of the larger
  # one
  # the strange coding is due to the fact that the original grid and the larger one
  # may differ at the level of roundoff error. 
  indX<-  c( 
              which.min( abs(gridList$x[1]  - gridListNew$x)),
              which.min( abs(gridList$x[nx] - gridListNew$x))
             )
  indY<-  c( 
              which.min( abs(gridList$y[1]  - gridListNew$y)),
              which.min( abs(gridList$y[ny] - gridListNew$y))
  )
              
  return( 
    list(gridListNew = gridListNew, gridList = gridList,
          indX= indX,
          indY= indY,
          np=np )
  )
}

#