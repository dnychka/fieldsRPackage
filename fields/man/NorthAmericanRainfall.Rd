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
\name{NorthAmericanRainfall}
\alias{NorthAmericanRainfall}
\alias{NorthAmericanRainfall2}
\docType{data}
\title{
Observed North American summer precipitation from the
historical climate network.
}
\description{
Average rainfall in tenths of millimeters for the months
of June, July and August (JJA) for the period 1950-2010 (old version, 1700+ locations)
or 1971-2023 (new version, 4800+ locations) and for North America. The old version uses monthly data from the Global Historical Climate Network (GHCN) version 2 and the new version uses Version 4. 
}
%\usage{data(NorthAmericanRainfall)}
\format{
 For both versions  the format is a list with components:
"longitude" "latitude"  "precip"    "elevation" "precipSE" "trend" "trendSE" "type" "x.s" "sProjection"
with elevation in meters, longitude as (-180,180), latitude as (-90, 90) and precipitaion in 1/10 mm 
( precip/254  converts to inches of rainfall). Note that the usual convention is to report a total amount of rainfall but here we give the mean over the three summer months. Multiple by 3 to get total precipitation over June through August.

\code{precip} is the intercept for 1980.5 (1997 for new version) when a straight line least squares regression is fit to 
each station's record. SE is the companion standard error from the least squares fit.
If the station is complete, then \code{precip} and \code{precipSE} will just be the mean and standard deviation adjusted for a linear trend. The estimated trend  \code{trend} and and its standard error \code{trendSE} are also included. 

 Also due to the centering, for complete data the intercept and trend estimate will be uncorrelated. The component \code{type} indicates whether the station has been "adjusted" (see below) or is still in
  "unadjusted" form. 
  
\code{x.s} is a useful transformation of locations into stereographic coordinates that reduces the
inflation of North Canada due to the usual lon/lat coordinates. Specifically it is found by: 
\preformatted{
  library(mapproj)
  xStereo<- mapproject( NorthAmericanRainfall$lon,NorthAmericanRainfall$lat, projection="stereographic")
  NorthAmericanRainfall$x.s<- cbind( xStereo$x, xStereo$y)
  hold<- .Last.projection()
   NorthAmericanRainfall$orientation<- hold$orientation
   NorthAmericanRainfall$projection <- hold$projection
}

}

\source{

The monthly data used to construct 
\code{NorthAmericanRainfall} this summary was generously provided by Xuebin Zhang, however,
the orignal source is freely available as the Global Historical Climate Network Version 2 Precipitation
quality controlled, curated and served by the US National Center for Environmental Information. 
The adjusted data from this archive has been modified from its raw form to make the record more homogenous. Heterogenities can come from a variety of sources such as a moving the station a short distance or changes in instruments. See the  National Centers for Envrionmental Information then Access Data  and then GHCN. 

The new version \code{NorthAmericanRainfall2} uses the GHCN "Version 4" data product. There are substantially more stations reported in this data set.
To reduce the size the stations for this R data set stations were selected if they had more than 80 percent of the months reported over the period 1971-2023. This subset amounts to  4893 locations. A simple linear regression was fit to each station record and the OLS trend, trend SE, trend line at 1997 is evaluated to get a "mean" rainfall estimate for this location. 

See the fields github repository for the source code and source data to create this data set. 

}
\examples{
data(NorthAmericanRainfall)
x<- cbind(NorthAmericanRainfall$longitude,  NorthAmericanRainfall$latitude)
y<- NorthAmericanRainfall$precip
quilt.plot( x,y)
world( add=TRUE)

Zstat<- NorthAmericanRainfall$trend / NorthAmericanRainfall$trendSE
quilt.plot( x, Zstat)

# with a better projection
library( maps)
bubblePlot( NorthAmericanRainfall$x.s, 
NorthAmericanRainfall$precip,
size=.3, col=rev(tim.colors()) )
map("world", 
projection=NorthAmericanRainfall$projection,
orientation= NorthAmericanRainfall$orientation,
add=TRUE, col="grey30")

# note station sparsity outside CONUS and lower CA.

data(NorthAmericanRainfall)
bubblePlot( NorthAmericanRainfall$lon, NorthAmericanRainfall$lat,
NorthAmericanRainfall$precip,
size=.3, col=rev(tim.colors()) )
world( add=TRUE, col="grey30")

data(NorthAmericanRainfall2)
bubblePlot( NorthAmericanRainfall2$lon,NorthAmericanRainfall2$lat,
NorthAmericanRainfall2$precip,
size=.3, col=rev(tim.colors()), xlim= c(-130,-50))
world( add=TRUE, col="grey30")

}
\keyword{datasets}
