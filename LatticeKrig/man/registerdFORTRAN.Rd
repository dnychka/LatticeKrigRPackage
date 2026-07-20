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

\name{registeredFORTRAN}
\alias{findnorm}
\alias{lkdist}
\alias{lkdistcomp}
\alias{lkdistgrid}
\alias{lkdistgridcomp}
\alias{lkdiag}
\docType{data}
\title{
Internal FORTRAN routines for working with grids and finding distances. 
}
\description{
These are objects of class  \code{FortranRoutine} and also 
\code{\link{NativeSymbolInfo}}
They provide information for compiled functions called with \code{.Call},
or \code{.Fortran}.
Ordinarily one would not need to consult these and they are used to make
the search among dynamically loaded libraries ( in particular the fields library
) have less ambiguity and
also be faster. These are created when the package/library is loaded 
are have their definitions from the compilation of  \code{init.c} in the
package source (src) directory. 
}
\format{
  The format is a list with components:
  \describe{
  \item{name}{The (registration ?) name of the C function.}
  \item{address}{See \link{NativeSymbolInfo}. }
  \item{dll}{Dynamically linked library information.}
  \item{numParameters}{Number of calling arguments in function.}
  }
  }

\details{
Registered routines are
\describe{
 \item{findnorm}{Finds marginal variance quickly for the case of a rectangular
 lattice,
 a specific Wendland basis, and constant a.wght at a given level. Called from 
 \code{LKRectangleFastNormalization}}
 \item{lkdist}{Euclidean distance between two sets of coordinates but 
 restricted to a
 maximum distance.Called from \code{LKDist}.}
 \item{lkdistcomp}{Same as lkdist but returns the component distances for each
 coordinate.Called by \code{LKDistComponents}}
 \item{lkdistgrid}{Finds all distances in x1 that are within delta of a set 
 of grid
 points. Results are returned in the row, column, value sparse matrix 
 format called from \code{LKDistGrid}.}
 \item{lkdistgridcomp}{Same as lkdistgrid but the component distances are found
 for
 each coordinate.\code{LKDistGridComponents}}
  \item{lkdiag}{Fills in diagonal and off diagonal elements
  in sparse matrix format.}
 }
 See \code{package_native_routine_registration_skeleton} for the 
utility used to create these data objects.
%%  ~~ If necessary, more details than the __description__ above ~~
}

\references{
For background on registering C, C++ and Fortran functions see 5.4 of
Writing R Extensions. 
}
\examples{
print(lkdistgridcomp)

}
\keyword{datasets}
