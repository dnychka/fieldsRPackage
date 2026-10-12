%#
%# fields  is a package for analysis of spatial data written for
%# the R software environment.
%# Copyright (C) 2026 Colorado School of Mines
%# 1500 Illinois St., Golden, CO 80401
%# Contact: Douglas Nychka,  douglasnychka@gmail.com,
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
\alias{panelPlotCommonAxes}
\alias{addTitleBox}
\title{
Creates a panel of plots with common axes using base R graphics.
}
\description{ 
This function creates panels of plots in a row/column
layout that can have common axes and a single legend or color bar.
Its main benefit is that one can 
customize a single plot in the panel to add additional features. 
All the code is in base 
R graphics, so it is useful as a graphics 
tutorial or as a starting point for modifying the code, 
and for the user comfortable with base R coding. 
This function is similar to the lattice package, but the intent is 
that it be easier to use and modify.
}

\usage{
panelPlotCommonAxes(m, n, xr, yr,
                               plotCode,
                               xlabels= "",
                               ylabels="",
                               outerXLabel=NULL,
                               outerYLabel=NULL,
                               PANESPACE=0,
                               LEGENDSPACE=c(0,0,0,0),
                               axes=TRUE,
                               addAxis=c(TRUE, TRUE, TRUE, TRUE),
                               alternateAxes=TRUE,
                               axisMargin=c( 3.1,3.1,3.1,3.1),
                               axisList=list( line=0),
                               axisLabelList=list(line=2),
                               outerLabelMargin= c(3.1,3.1,0,0),
                               outerLabelList=list(line=2),
                               colBg="grey90",
                               colGrid="grey98",
                               boxLwd= 1,
                               titleBoxHeight= 0,
                               titleBoxList= list( 
                               boxColor = "lightgreen",
                               boxBorder= "black",
                               lwd=boxLwd),
                               verbose=FALSE
                               )
                               
      addTitleBox(  boxHeight = 1.2,
                    boxColor="green3",
                    boxBorder="black",
                    lwd=1.0
                    )
      
}
\arguments{
  \item{m}{ Number of panel rows.

}
  \item{n}{ Number of panel columns.

}
  \item{xr}{ Range of the horizontal axis/axes. If a two-column matrix, row J 
  of the matrix gives the axis range for column J of the panel. 

}
  \item{yr}{
   Range of the vertical axis/axes. If a two-column matrix, row I 
  of the matrix gives the axis range for row I of the panel.
}
  \item{plotCode}{
  A list of expressions for plotting in each panel. See Examples and Details.

}
  \item{xlabels}{ Labels for the horizontal axes.

}
  \item{ylabels}{Labels for the vertical axes.

}
  \item{outerXLabel}{ Outside horizontal label. 

}
  \item{outerYLabel}{Outside vertical label.

}
  \item{PANESPACE}{  Spacing between and around panes in units of a character width.
  

}
  \item{LEGENDSPACE}{ A vector of four line spaces in the same order as the
  \code{mar} and \code{oma} margins. 
  This is the spacing added around the panel in units of a character width.

}
  \item{axes}{ If TRUE, add axes to the panel.

}

  \item{addAxis}{ A logical vector of four values specifying which sides should 
  get axes.

}
\item{alternateAxes}{If TRUE, only every other axis is labeled with numerical
values. 
This 
makes for less busy panel margins. }

  \item{axisMargin}{Space added for tick marks and the axis labeling.}

  \item{axisList}{ A list of other axis arguments. See \code{\link{mtext}} for 
  the
  possible choices. The default is a list with just a
   line component, \code{list(line=0)}, but the user can add to this for
   more modifications.
}
  \item{axisLabelList}{ Same as axisList but for the outer labels. 

}
  \item{outerLabelMargin}{ Same as axisMargin but for the outer label.

}
  \item{outerLabelList}{ Same as \code{axisLabelList} but for the outer label.

}
  \item{colBg}{ Background color for the plots. }

\item{colGrid}{ Color of grid lines}

\item{titleBoxHeight}{Height of the title box in units of "lines".
     This argument is passed to \code{addTitleBox} }
     
\item{ titleBoxList}{A list with other arguments to send to the 
\code{addTitleBox} function.}

\item{boxHeight}{Height of the title box in units of "lines"
     (the same units used to position 
  text and axes in the margins). }
  \item{verbose}{ 
  If TRUE, print some information about the layout. This is for debugging. 
  }
  \item{boxColor}{Fill color of the box.}
  \item{boxBorder}{Color of the box border.}
  \item{boxLwd}{Line width for the box around each plot and for the title box.}
  \item{lwd}{Thickness of the border lines.}
}



\details{
This function was created to make plots in an 
mXn panel with consistent spacing and to add 
common axes on the sides. 
The idea is that if two plots side by side have a 
common vertical axis, 
then 
it should 
not be replotted in between the two plots. Instead, the second axis
should be added on the right side -- or omitted.
We have resisted making this into a
more polished function
for two reasons. 
First, panel plots typically need customization, and 
so handling all kinds of
flexibility would make
a general function complex and its code more difficult to follow. 
Even so, we are not totally happy with how this version reads!
There are useful frameworks such as \code{ggplot2} 
and \code{lattice}
that provide this kind of added functionality, but in doing so they require 
some effort to 
deviate from the default settings. 

The second reason is our 
realization that
the basic elements of
grouping plots together and adding common axes on 
the outsides are really not
too complicated to code up. 
The goal of this function 
was to make the whole process accessible  
using just base R graphics. 
In this way, it is easier to modify if one wants 
some additional 
customization. 

\strong{Difference from the lattice package } The \code{lattice} package is an
elegant
framework following Cleveland's work on statistical graphics. The main 
idea is that the same graphics function is applied to subsets of a data set,
for example, making separate plots for different categories. In this
way, the "for" loop over the separate plots is implicit -- determined by the 
number of different categories. 
The challenge is that customizing a specific plot in this group, 
e.g., adding an arrow or an additional line for explanation, is hard. And, of 
course, while the default spacing and layout are thoughtful, inevitably they
will 
need to be tweaked for the final version in a publication or for a talk. 
This function takes the opposite point of view in that you specify each plot in
the panel separately to have full control. This can be tedious for a panel 
with many plots and is a disadvantage. However, see the note at the end of this
section on how to mimic lattice and draw multiple plots with a single graphics 
function. 

\strong{Space}
The many arguments for spacing are included because, inevitably,
to get a final
figure, one needs to tweak 
the spaces around a plot. Sometimes this happens simply 
when converting the
plot viewed in an interactive session to a pdf for a paper or talk.

\strong{Basic idea of the function}
 The figure panel is mXn -- m rows and n columns -- with each plot referred to
 as a pane. 
 The main part of the 
function 
loops through plotting these \emph{row by row}. 
So if the panel is 3X2, then 
it plots in the order (1,1), (1,2), (2,1), ..., (3,2).  
Each pane of the panel is initially set up blank with a generic call
to the \code{plot} function using the limits in 
\code{xr} and \code{yr}, without any axes or plotting. The expressions in 
\code{plotCode} then plot into this pane. Just use R graphics functions 
that \emph{add} to a plot, e.g., \code{points, lines, matpoints}. Functions 
such as \code{image}, \code{imagePlot}, \code{contour}, and \code{map}
with the argument \code{add=TRUE} also avoid creating a new plot. 
The real estate is 
already set up for you!   

\strong{Specifying the plotting code}
So what about the plotting code? This is, at first glance, a bit strange (maybe
at second glance too).
It is a list with each component being an \emph{expression}.
Each expression can have
multiple plotting/R code lines separated by commas. So, for example, if one 
planned to
have 4 scatterplots with common axes, based on \code{xr} and 
\code{yr}, and using \code{points}:
\preformatted{
plotCode<- list( 
                 expression( points(x1,y1)),
                 expression( points(x2,y2)),
                 expression( points(x3,y3)),
                 expression( points(x4,y4),
                       title("fourth plot"))
                 )
}

Here we added another plot function to the fourth plot -- adding a title
to show 
how this is done. If a 
component of the list is \code{NULL}, then that plot is skipped.

\strong{Exporting the for loop indices} The function has the variables \code{I}
(row index), \code{J}
(column index), and \code{counter} (plots in order) to keep track of the 
current plot. These can be used to avoid creating separate plotting code 
for each plot or to customize a particular case. 
For example, 
a single expression can be used with these variables to create
the different plots. This is more in the style of the \code{lattice} package;
see the example below for its use. 

\strong{Adding legends, color bars and other annotation.}
Finally, we note that often one wants to draw a panel and then also add a 
common legend, color bar, or both. Again, we want to do this using base
R graphics
in a simple way that is easy to customize.

We suggest the following 
steps: 1) Use LEGENDSPACE to reserve space for the legends, typically to
the right or at the bottom. 2) Draw your 
plot. 3) Use \code{addLegend()} to reset your figure region to a
single [0,1]X[0,1] plotting region.
 
4) Add your extra legend or color bar using the [0,1] ranges to place
these elements. 
The reason for the figure reset is 
to make it
easier to think about plotting coordinates for locating the legend. 
See the examples below. 
The advantage of step 3 is that both the color bar  
and legend placement work nicely with 
a unit square figure region.  

\strong{Adding a strip of color to the top}
Lastly, we have included a handy function, \code{addTitleBox}, intended
for the specific purpose of adding a strip of a different color at the top of 
a plot to hold a title. This is the style of the \code{lattice} package, and
it would be annoying to keep coding this if the function were not around.
See the examples below. BTW, for a humble statistician this was very tricky to 
code up! In this function, note how the line spacing in the margin is translated
to user coordinates. 

This 
function is also a good example of the care needed using base R graphics. 
The function turns off the clipping that prevents plotting outside of the
plot region. Then it adds a colored rectangle, still using the 
plot coordinate system. After drawing the rectangle, the clipping is 
turned back on. In other circumstances, this strategy of plotting outside 
the normal plot rectangle might 
easily go awry and give strange results. 

}

\author{
Douglas Nychka

}

\seealso{
\code{\link{imagePlot}},\code{\link{addLegend}}, \code{\link{addColorBar}}

}

\examples{
# A 2x2 panel to show sampling variability of the N(0,1)
# scatter plots  from 4 samples 

set.seed(123)
x<- matrix( rnorm(50*6), ncol=6,nrow=50)
y<- matrix( rnorm(50*6), ncol=6, nrow=50)

# begin list of expressions
plotCode<-  list(
expression( points(x[,1], y[,1],pch=1)
),
expression( 
points(x[,2], y[,2],pch=2)
),
expression( 
points(x[,3], y[,3],pch=3)
),
expression( 
points(x[,4], y[,4], pch=4)
)
)# end list 

#
xr<- range( c(x))
yr<- range( c(y))

############################################################
# a basic plot relying on all the default choices 
############################################################
panelPlotCommonAxes(
m=2,n=2,
xr,yr,
plotCode
)
############################################################
# same plot but using the indexing within the for loop
############################################################
plotCode1<-  list(
expression( points(x[,counter], y[,counter],pch=counter, col=counter)
)
)
# plotting code reused for each pane of the panel
# because this list has length 1
fields.style() # handy reset of some  default graphics sizes and colors
panelPlotCommonAxes(
m=2,n=2,
xr,yr,
plotCode1
)

############################################################
# add some titles in a color strip
############################################################

plotCode2<-  list(
expression( points(x[,1], y[,1],pch=1),
title("Case 1", adj=.2, line=.2)
),
expression( 
points(x[,2], y[,2],pch=2),
title("Case 2", adj=.2, line=.2)
),
expression( 
points(x[,3], y[,3],pch=3),
title("Case 3", adj=.2, line=.2)
),
expression( 
points(x[,4], y[,4], pch=4),
title("Case 4", adj=.2, line=.2)
)
)# end list 

panelPlotCommonAxes(
m=2,n=2,
xr,yr,
plotCode2,
titleBoxHeight=1.2
)


############################################################
# adding a legend, turning off some axes, and adding labels 
# also mess with the x ranges so there are different axes
############################################################

# different ranges for x -- something harder to do in the lattice package
xr<- rbind( c( -3,3),
            c(-4,4))
#  

panelPlotCommonAxes(m=2,n=2,
# ranges 
xr,yr,
# the code for the four plots 
plotCode2,
# outer labels
outerXLabel ="X",
outerYLabel= "Y",
# picky about which axes are plotted 
addAxis=c(TRUE, TRUE, FALSE, FALSE),
alternateAxes=FALSE,
# room for a title
titleBoxHeight=1.2,
# amount of space (in lines) for the axes on the sides of the panel
axisMargin=c( 4.1,4.1,0,0),
# extra room for outer labels 
# tight margin
outerLabelMargin= c( 1,1,0,0), 
# 5 extra spaces on the right reserved for the legend
LEGENDSPACE=c(0,0,0,5),
boxLwd=1
)
#
# Now the base R magic of adding the legend 
# reset the figure region to (0,1) X (0,1) coordinates 
addLegend()
cXY<- par()$cxy # character size for these dimensions
# this is useful for thinking in units of "lines"
# now call the usual legend function ...
legend(
"topright",
# set the number of lines in from the top-right corner
# adjust as needed ...
inset=c(cXY[1]*1.5,cXY[2]*3.5), #  about 1.5 and 3.5 lines from corner
legend=c("do","re","mi", "fa"), 
pch= 1:4, cex=1.2)

###

\dontrun{
dataLattice<- data.frame( X=c(x), Y=c(y), cases=  
as.factor(rep(1:6, rep(nrow(x),6))) )

library( lattice)

# clean call to lattice using all its defaults 
# but plot in reverse order  ( a table)

xyplot( Y~X|cases, data=dataLattice,
main= "lattice panel plot ", as.table=TRUE )

# use the model of lattice/trellis: the same graphics applied to each subset
plotCode3<-  list(
expression( points(x[,counter], y[,counter],
col="#0080ff", #medium sky blue default for lattice
cex=1.5),
title(format(counter),  line=.2))
)

xr<- range( c( x))
yr<- range( c( y))

# quick with defaults 
 panelPlotCommonAxes(m=2,n=3,
                     xr,yr,
                     plotCode3,
                     outerXLabel ="X",
                     outerYLabel= "Y",
                     titleBoxHeight=1.2
 )
 
# matching lattice 
# note some futzing around with spacing and color
# grid and background are turned off -- not our first choice
# 
 panelPlotCommonAxes(m=2,n=3,
                     xr,yr,
                     plotCode3,
                     outerXLabel ="X",
                     outerYLabel= "Y",
                     axisMargin=c(1,1,1,1) ,
                     outerLabelMargin= c(4,4,4,4), #tight margin
                     titleBoxHeight=1.2,
                     outerLabelList=list(line=3),
                     titleBoxList= list( 
                       boxColor = "grey90",
                       boxBorder= "black",
                       lwd=1), 
                       colBg="white",
                       colGrid="white",
                     axisLabelList=list(line=3),
 )
 
 # this is complicated but most shows how to control
 # exactly where this will go
 mtext( text= "fields panel plot ", outer=TRUE, 
 side =3, line = 3, 
 font = 2, cex = 1.2)
# BTW typically in a publication the caption subsumes any title in the figure

 
# and doing this in ggplot2 
# not all the feature are easily reproduciable 
#-- we  suspect the ggplot folks would prefer the plot the way it is!

library(ggplot2)
ggplot(dataLattice, aes(x = Y, y = X)) +      
  geom_point(shape = 1, colour = "#0080ff" ) +
  facet_wrap(~ cases, nrow = 3, ncol = 2) +      
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"),      
        strip.background = element_rect(fill = "grey80"))
}
}

\keyword{hplot}
