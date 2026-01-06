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
unrollXMatGrid<- function( grid.list, XMatGrid){
  if( is.null(XMatGrid)){
    return(XMatGrid)
  }
  if( is.list( XMatGrid) ){
     if( any(grid.list[[1]] != XMatGrid[[1]]) | any(grid.list[[2]] != XMatGrid[[2]]) ){
         stop("grid list does not match grid for covariates")
       }  
# wipe out the x and y components of XMatGrid 
  XMatGrid<- XMatGrid$z
  }
# check dimensions
    XMatdim<- dim( XMatGrid)
      nx<- length( grid.list[[1]])
      ny<- length( grid.list[[2]])
      if( (XMatdim[1] != nx) | (XMatdim[2] != ny) ){
         cat( "XMat:", XMatdim, "grid", nx, ny, fill=TRUE )
         stop( "Dimension of XMatGrid does not match dimensions of location grid list.")
      }
# reshape as a matrix where rows index locations.
# Note that this works whether XMatdim[3] exists or not! 
      return( matrix( c(XMatGrid),  nrow= XMatdim[1]*XMatdim[2] ))
 }
