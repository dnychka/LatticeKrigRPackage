%%BEGIN HEADER
%
% LatticeKrig is a package for analysis of spatial data written for
% the R software environment.
% Copyright (C) 2026 Colorado School of Mines
% 1500 Illinois St., Golden, CO 80401
% Contact: Douglas Nychka,  douglasnychka@gmail.com,
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% A copy of the GNU General Public License is included
% along with the R software environment if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
% or refer to  http://www.r-project.org/Licenses/GPL-2
%
%%END HEADER
\name{KrigingExampleData}
\alias{KrigingExampleData}
\docType{data}
\title{Synthetic data for kriging examples}
\description{
This is a sample of 50 (x,y) ordered pairs: the x values are in sorted order,
drawn from the interval [0,1]. The smallest x value is exactly 0, and the largest is exactly 1.
The y values are the function y=9*x*(1-x)^3 evaluated at each x value, with a small measurement error (normal with mean 0, standard deviation 0.01) added on.
}
\usage{
 data(KrigingExampleData)
 }
\details{
The following code was used to generate the data:

\preformatted{
set.seed(123)
x <- c(0,sort(runif(48)),1 )  
y <-  9*x*(1-x)^3  + rnorm(50, sd=0.01)
KrigingExampleData<- list( x=x, y=y)
save( KrigingExampleData, file="KrigingExampleData.rda")
}

}
\format{The data in \code{KrigingExampleData} is listed in two vectors, \code{x} and \code{y}.}
\examples{
\dontrun{
data(KrigingExampleData)
x<- KrigingExampleData$x
y<- KrigingExampleData$y
plot(x,y)
kFit <- LatticeKrig(x, y)
xGrid <- seq(0,1,0.001)
lines(xGrid, predict(kFit, xGrid), col='red')
}
}
