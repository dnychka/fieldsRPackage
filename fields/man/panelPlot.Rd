%#
%# fields  is a package for analysis of spatial data written for
%# the R software environment.
%# Copyright (C) 2024 Colorado School of Mines
%# 1500 Illinois St., Golden, CO 80401
%# Contact: Douglas Nychka,  douglasnychka@gmail.edu,
%#
%# This program is free software; you can redistribute it and/or modify
%# it under the terms of the GNU General Public License as published by
%# the Free Software Foundation; either version 2 of the License, or
%# (at your option) any later version.
%# This program is distributed in the hope that it will be useful,
%# but WITHOUT ANY WARRANTY; without even the implied warranty of
%# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%# GNU General Public License for more details.
%#
%# You should have received a copy of the GNU General Public License
%# along with the R software environment if not, write to the Free Software
%# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%# or see http://www.r-project.org/Licenses/GPL-2
%##END HEADER
%##END HEADER

\name{panelPlot}
\alias{panelPlotExample}
\title{
Tips for a panel plot using base R graphics. 
}
\description{ This help section gives some tips to create panels of plots in a row/column
layout that can have common axes and a single legend or color bar. It also shows how to 
customize a single plot in the panel to add additional features. All the code  is in base 
R graphics and so is useful as a graphics tutorial and for the user comfortable with base R
coding. 
}

\usage{
panelPlotExample()
}

\details{
The function is really a fairly elaborate "script" that creates a  
 3X2 panel of plots. Usually there will not be as much variation among the individual
 plots. That is the whole point of organizing them in a matrix format!
However, the intent is for the user to copy this code and modify it according.
The commenting should help
in naviagting the different parts. We have resisted making this into a more polished function
for two reasons. 
First, panel plots typically need customization and so handling all kinds of
flexibilty would make
a general function complex.  There are useful frameworks such as \code{ggplot2} 
and \code{lattice}
that provide this kind of function but in doing so require some effort  to 
deviate from the default settings. The second reason is our realization that 
that basic elements of
grouping plots together and adding  common axes on the outsides is really not
too complicated to code up. 
We implement this as a \code{for} loop over the panel plots
with simple R coding and base graphics. Making this explicit then helps
adding particular features to 
specific plots in the panel -- an option that we found awkward in higher level packages for 
creating panel plots. See the comments in our example. For a data science student 
it is useful 
to become familar with building a  complex figure from  a sequwence of 
simpler graphical elements and 
to be able to do this using base R graphics.  

The most exotic part of this example is the use of the the 
\code{eval} function to 
evaluate the plot functions for each separate panel. The plotting code is setup as 
components of
a list, \code{plotCode}, at the beginning of the script. Within
the \code{for} loop over the plots these
graphics commands  are just evaluated. This helps to
avoid modifying the code within the \code{for} loop keeping  it 
 unchanged for different figures.
The comments in the function indicate the parts to change for each plot and those 
that usually are not touched. 

Here we just generate trivial vectors to plot although the reader is probably aware that any data object 
refered to in this function will be looked for in the work space. Or to be more rigorous add some
figure specific arguments to the function and just pass the object as an argument. This would help to 
isolate drawing the figure from stray versions of the data objects in the workspace with the same name. 

To run the example  use \code{panelPlot()} from this package. Unfortunately the default package 
loading is to strip comments from functions and so \code{panelPlot} has
its comments removed. You can find the commented version in 
\url{https://github.com/dnychka/fieldsRPackage/tree/main/fields/R}
in the file \code{panelPlotExample.R} Use this code to make a copy of the function,  
rename it, then modify as needed. Enjoy!

}

\examples{
# run as a function
panelPlotExample()
# see 
# https://github.com/dnychka/fieldsRPackage/tree/main/fields/R/panelPlotExample.R
# for the function with comments 
#
# We usually would dump finsihed figures out to a pdf file using  the pdf output function
# to control the size and reproducibilty.
# E.g.
\dontrun{
   pdf("panelPlotExample.pdf", width=8, height=8 )
   panelPlotExample()
   dev.off()
}

}

\keyword{datasets}
