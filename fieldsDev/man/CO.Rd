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

\name{Colorado Monthly Meteorological Data}
\alias{COmonthlyMet}
\alias{CO.elev}
\alias{CO.id}
\alias{CO.loc}
\alias{CO.names}
\alias{CO.ppt}
\alias{CO.ppt.MAM}
\alias{CO.tmax}
\alias{CO.tmax.MAM}
\alias{CO.tmin}
\alias{CO.tmin.MAM}
\alias{CO.years}
\alias{CO.ppt.MAM.climate}
\alias{CO.tmax.MAM.climate}
\alias{CO.tmean.MAM.climate}
\alias{CO.tmin.MAM.climate}
\alias{CO.elevGrid}
\alias{CO.Grid}

\title{Monthly surface meteorology for Colorado 1895-1997}

\description{

Source:
 This is a group of R data sets for monthly min/max temperatures and
precipitation over the period 1895-1997. It is a subset extracted from
the more extensive US data record.
Temperature is in degrees C and
precipitation is total monthly accumulation in millimeters. Note that
minimum (maximum) monthly temperature is the mean of the daily minimum
(maximum) temperatures.

Data domain:

A rectangular lon/lat region [-109.5,-101]x [36.5,41.5] larger than the
boundary of Colorado comprises approximately 400 stations. Although
there are additional stations reported in this domain, stations that
only report precipitation or only report temperatures have been
excluded. In addition, stations that have mismatches between locations
and elevations from the two metadata files have also been excluded. The
net result is 367 stations that have colocated temperatures and
precipitation.

}

\format{
This group of data sets is organized with the following  objects:
\describe{
\item{CO.info}{A data frame with columns: station id, elev, lon, lat, station name}


\item{CO.elev}{elevation in meters}

\item{CO.elevGrid}{An image object being elevation in meters on a 4 km grid covering Colorado. }

\item{CO.id}{ alphanumeric station id codes}

\item{CO.loc}{locations in lon/lat}

\item{CO.Grid}{Just the grid.list used in the CO.elevGrid.}

\item{CO.ppt CO.tmax CO.tmin}{Monthly means as three dimensional arrays ( Year, Month, Station).
Temperature is in degrees C and precipitation in total monthly
accumulation in millimeters.}

\item{CO.ppt.MAM CO.tmax.MAM CO.tmin.MAM}{Spring seasonal means
(March, April, May) as two dimensional arrays
   (Year, Station).}

 \item{CO.MAM.ppt.climate CO.MAM.tmax.climate CO.MAM.tmin.climate}{Spring seasonal means
(March, April, May) by station for the period 1960-1989. If fewer than 15 years are present over this period, an NA is recorded.
No detrending or other adjustments have been made for these mean estimates.
}
}


}

\examples{

data(COmonthlyMet)

#Spatial plot of 1997 Spring average daily maximum temps
 quilt.plot( CO.loc,CO.tmax.MAM[103,]  )
 US( add=TRUE)
 title( "Recorded MAM max temperatures (1997)")

# min and max temperatures against elevation

matplot( CO.elev, cbind( CO.tmax.MAM[103,], CO.tmin.MAM[103,]),
  pch="o", type="p",
  col=c("red", "blue"), xlab="Elevation (m)", ylab="Temperature (C)")
title("Recorded MAM max (red) and min (blue) temperatures 1997")

#Fitting a spatial model:
obj<- Tps(CO.loc,CO.tmax.MAM.climate,  XMat= CO.elev )
\dontrun{
out<- spatialProcess(CO.loc,CO.tmax.MAM.climate,
          smoothness=1.0, XMat= CO.elev)
surface( out)
	  }


}

\keyword{datasets}
% docclass is data
% Converted by Sd2Rd version 1.21.
