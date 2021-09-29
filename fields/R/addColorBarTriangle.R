# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2018
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
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
addColorBarTriangle <- function(lowerColor=NULL,
                                upperColor=NULL,
                                horizontal=TRUE) {
usr<- par()$usr
plt<- par()$plt

par(xpd=TRUE)
if(horizontal){
  deltaY<-  (usr[4] - usr[3])
  unitX<-   (usr[2] - usr[1])/(plt[2] - plt[1])
  deltaX<-  unitX*(plt[4] - plt[3])
}
else{ 
  deltaX<-  (usr[2] - usr[1])
  unitY<-   (usr[4] - usr[3])/(plt[4] - plt[3])
  deltaY<-  unitY*(plt[2] - plt[1])
}
#
if( !is.null(upperColor) ){
  if(horizontal){
    triangleUpper<- rbind( 
                     c( usr[2],              usr[3] ),
                     c( usr[2]+ deltaX,      usr[3] + deltaY/2 ),
                     c( usr[2],              usr[4])
    )
  }
  else{
    triangleUpper<- rbind( c( usr[1],          usr[4] ),
                       c( usr[1] + deltaX/2,  usr[4] + deltaY ),
                       c( usr[2],          usr[4])
    )
  }
  polygon( triangleUpper, col=upperColor, border=upperColor)
  lines(triangleUpper)
}
#
if(!is.null(lowerColor)){
  if(horizontal){
    triangleLower<- rbind( c( usr[1],          usr[3] ),
                   c( usr[1]- deltaX,  usr[3] + deltaY/2 ),
                   c( usr[1],          usr[4])
    )
  }
  else{
    triangleLower<- rbind( c( usr[1],          usr[3] ),
                           c( usr[1]+ deltaX/2,  usr[3] - deltaY ),
                           c( usr[2],          usr[3])
    ) 
  }
  polygon( triangleLower, col=lowerColor, border=lowerColor)
  lines(triangleLower)
}
#
par(xpd=FALSE)
}
