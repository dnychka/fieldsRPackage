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
\name{AudiA4}
\alias{AudiA4}
\docType{data}
\title{Used car prices for AudiA4 from cars.com}
\description{
These data are a small regression data set that was collected in the process
of trying to find a "deal" on a used Audi A4. (Note that a good deal could be
taken as a car where the regression residual is very negative.)
}
 \usage{ data(WorldBankCO2)}
\format{
  This is a 1565X4 matrix with the rows being individual cars listed for sale on
  cars.com and the columns being the covariables:
 
\code{ "year"     "price"    "mileage"  "distance"}
\itemize{
\item year: Model year 
\item price: Price in dollars
\item mileage: Mileage. 
\item distance: Distance of car from Boulder, Colorado
}

}
\examples{
data(AudiA4)
fields.style()
plot( AudiA4$mileage/1000,AudiA4$price/1000, pch=16, 
xlab="Mileage (1000s)", 
ylab="Price (1000s)", col="grey" )
fit<- Tps(AudiA4$mileage/1000,AudiA4$price/1000 )
mileageGrid<-seq(0, 1700 )
priceHat<- predict( fit, mileageGrid)
lines(mileageGrid,priceHat, col="orange3"  )
}

\keyword{datasets}
