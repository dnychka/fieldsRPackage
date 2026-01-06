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
"tim.colors" <- function(n = 64, alpha = 1.0) {
    # tims original 64 color definition ala MATLAB:
    # orig <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", 
    #     "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", 
    #     "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", 
    #     "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", 
    #     "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", 
    #     "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", 
    #     "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", 
    #     "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", 
    #     "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", 
    #     "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", 
    #     "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", 
    #     "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", 
    #     "#AF0000", "#9F0000", "#8F0000", "#800000")
  # switching to the turbo viridisLite choices
  # turbo(64)
  orig <- 
    c( 
   "#30123BFF", "#351E58FF", "#392A74FF", "#3D358CFF", "#4040A3FF",
   "#424CB6FF", "#4457C8FF", "#4662D7FF", "#466CE4FF", "#4777EFFF",
   "#4681F7FF", "#458BFDFF", "#4195FFFF", "#3C9FFDFF", "#36AAF9FF",
   "#2EB3F3FF", "#27BEE9FF", "#20C7E0FF", "#1BD0D5FF", "#18D7CAFF",
   "#18DEC0FF", "#1AE4B6FF", "#20EAACFF", "#2AEFA1FF", "#35F393FF",
   "#44F786FF", "#53FA79FF", "#62FC6BFF", "#72FE5EFF", "#81FF52FF",
   "#90FF48FF", "#9EFD40FF", "#A8FB39FF", "#B3F836FF", "#BDF434FF",
   "#C7EF34FF", "#D2E935FF", "#DBE336FF", "#E3DB38FF", "#EBD339FF",
   "#F1CB3AFF", "#F6C33AFF", "#FABA39FF", "#FCB136FF", "#FEA732FF",
   "#FE9B2DFF", "#FE8F29FF", "#FC8323FF", "#F9771EFF", "#F66B19FF",
   "#F25F14FF", "#EC5410FF", "#E84A0CFF", "#E1420AFF", "#DB3A07FF",
   "#D33205FF", "#CB2A04FF", "#C22402FF", "#B81D02FF", "#AD1701FF",
   "#A11201FF", "#950D01FF", "#880802FF", "#7A0403FF" 
   )
  
  designer.colors( n, col=orig, alpha=alpha)
}
