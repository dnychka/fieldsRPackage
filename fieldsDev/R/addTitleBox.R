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
addTitleBox<- function(  boxHeight = 1.2,
                    boxColor="green3",
                    boxBorder="black",
                    lwd=1.0
                    )
{box
#box height is in units of the margin spacing.
# compute the spacing of a line in user coordinates
  usrCoord<- par("usr")
# height of one margin line, in inches
  lineInches<- par("mai")[3] / par("mar")[3]                   
 # in user y-units (bottom/top margins)
  lineY  <- lineInches * (usrCoord[4]- usrCoord[3]) / par("pin")[2] 
#  lineX <-  lineInches * (usrCoord[2]- usrCoord[1]) / par("pin")[1]
   newCoord<- c(
                usrCoord[1] ,
                usrCoord[2] ,
                usrCoord[4] ,
                usrCoord[4] + boxHeight*lineY
                ) 
  # allow plotting in usr coordinates outside of the 
  # plotting rectangle -- otherwise this would be clipped (omitted)
  par( xpd=NA) # NA allows plotting everywhere on the device.
  #box with color
  rect(newCoord[1],newCoord[3],newCoord[2],newCoord[4],
       col=boxColor, border=boxBorder, lwd=lwd)
  # now reset the clipping so plotting outside is not done. 
  par( xpd=FALSE)
  
}
