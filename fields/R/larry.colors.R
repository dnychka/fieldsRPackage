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
larry.colors<- function(){
  ctemp<- matrix(  c(
  182,  106,   40,
  205,  133,   63,
  225,  165,  100,
  245,  205,  132,
  245,  224,  158,
  255,  245,  186,
  255,  255,  255,
  205,  255,  205,
  153,  240,  178,
   83,  189,  159,
  110,  170,  200,
    5,  112,  176,
    2,   56,   88 ),
  ncol=3, byrow=TRUE)
  ctemp<- ctemp/255
  rgb(ctemp[,1], ctemp[,2], ctemp[,3])
}
