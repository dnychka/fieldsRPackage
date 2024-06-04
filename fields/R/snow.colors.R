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
#
# from development -- the script to tease out colors from 
# piecewise linear color scale. 
#library( fields)
#look<- read.csv( "cmap_sca_256.csv")
#look<- as.matrix( look)
#rawRGBTable<- look/255
#x<- 1:255
# y1<- splint( x,look[,1], x,
#              df=15, derivative=2 )
# y2<- splint( x,look[,2], x,
#             df=15, derivative=2 )
# y3<- splint( x,look[,2], x,
#              df=15, derivative=2 )
# f<- y1^2+y2^2+y3^2
# indm<- 1:253
# ind0<- 2:254
# indp<- 3:255
# maxInd<-  (f[indp] < f[ind0]) & (f[indm]< f[ind0])
# temp<- t(col2rgb( colors()))
# colorsFound<-  rbind( look[1,],
#                  look[maxInd,],
# look[255,] )/255
# HexColors<-  rgb( colorsFound[,1],
#                  colorsFound[,2],
#                 colorsFound[,3])
# rawTable<- rgb( look[,1], look[,2],look[,3])
# snowColorsMap<- list( x= c(0,x[maxInd], 255)/255,
#                      col= HexColors,
#                      rawRGBTable=rawRGBTable)
#image(matrix(1:length(HexColors)), col=snow.colors() )
 
snow.colors<- function( n=256, alpha=1){
  x<- c(
    0,   7, 31,  63,  94, 127, 159,
    190, 222, 243, 248, 255)/255
 
  col<- c("#091E5A", "#0E2266", "#243493",
          "#225DA8", "#1D8FC0", "#41B6C5", 
          "#7FCEBC", "#C5E9B5", "#EDF8B2",
          "#F9FDE3", "#FCFEEF", "#FFFFFE") 
  colorTable<- designer.colors(n-1,
                col=col,
                   x=x ,
                   alpha = alpha)
#grey added to represent minimum value  
  withGrey<- c( "#808080", colorTable)
  return( withGrey )
 } 
 
 
 
 
 
 
