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

\name{LatticeKrig internal}

\alias{create_awght_obj}
\alias{generate_param_field}
\alias{plot_cov_surface}

\title{
  LatticeKrig internal and secondary functions 
}
\description{
Listed below are supporting functions for some of the help file examples. 
}
\usage{
plot_cov_surface(LKinfo, 
    real_field,
    loc_x, 
    loc_y, 
    plotting = TRUE, 
    cov_contour = TRUE,
    cov_lwd = 1.3,
    nrows_advanced,
    sGrid_advanced
)

create_awght_obj(
    sidelen, 
    awght_config,
    awght_override = FALSE,
    awght_override_val,
    theta_config, 
    theta_override = FALSE,
    theta_override_val,
    rho_config,
    rho_override = FALSE,
    rho_override_val,
    plotting = FALSE,
    data_rows,
     awghts,
    low_awghts,
    high_awghts
)

generate_param_field(
    paramtype,
    config,
    sGrid,
    rows,
    columns,
    n,
     awghts,
    low_awghts,
    high_awghts)
}
\keyword{internal}
% docclass is function
% Converted by Sd2Rd version 1.21.
