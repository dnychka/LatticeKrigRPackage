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

\name{LKinfoCheck}
\alias{LKinfoCheck}
\alias{LKinfoCheck.default}
\alias{LKinfoCheck.LKRectangle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Check the LKinfo object
%%  ~~function to do ... ~~
}
\description{
This method performs some simple checks on the LKinfo object before it is 
returned by LKrigSetup. One benefit is that it checks for the essential
components in each of the parts of the LKinfo object. 
Currently there is only a default method supplied. 

Here is an example of how the checks work
\preformatted{
	LKinfo<- LKrigSetup( cbind(c(0,1)), LKGeometry="LKInterval",
    a.wght=5, nlevel=2, nu=1.5, NC=4)
    LKinfo$alpha
# corrupt it  
    LKinfo$alpha <- list( 1,.5,.25)
try( LKinfoCheck( LKinfo) )
}

}
\usage{
LKinfoCheck(object, ...)
\method{LKinfoCheck}{default}(object, ...)
\method{LKinfoCheck}{LKRectangle}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  	An LKinfo object to be checked.
}
  \item{\dots}{
  	Additional arguments to method. 
}
}

\value{
There is no returned value. The benefit is the side effect of
an error message  and a "stop if the are problems in the LKinfo object. 
}

\author{
	Doug Nychka
}

\seealso{
\code{\link{LKrigSetup}}
}
\keyword{ spatial }

