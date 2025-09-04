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
"Mines.colors" <- function(color=1:14) {
  # tims original 64 color definition definition:
  orig <- c(
    "#09396C",
    "#CC4628",
    "#21314d",
    "#75757D",
    "#0272DE",
    "#F0F600",
    "#80C342",
    "#F1B91A",
    "#879EC3",
    "#AEB3B8",
    "#57A2BD",
    "#CFDCE9",
    "#B42024",
    "#81848A")
    nm <- c(
      "blasterBlue",
      "coloradoRed",
      "darkBlue",
      "darkGray" ,
      "earthBlue",
      "energyYellow",
      "envrGreen"  ,
      "goldenTech"  ,
      "lightBlue"  ,
      "lightGray"  ,
      "mutedBlue" ,
      "paleBlue" ,
      "redFlannel",
      "silver"
    )
    names(orig) <- nm
    orig[color]
}
