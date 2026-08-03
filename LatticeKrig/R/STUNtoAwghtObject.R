##BEGIN HEADER
#
# LatticeKrig is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2026 Colorado School of Mines
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
# A copy of the GNU General Public License is included
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or refer to  http://www.r-project.org/Licenses/GPL-2
#
##END HEADER

STUNtoAwghtObject <- function(
    gridlist,
    # if one already has the parameters in R
    kappa2 = NULL,
    theta = NULL,
    rho = NULL,
    # if the parameters are in an h5 file from STUN
    file_path = NULL,
    dataset_name = NULL,
    index = NULL,
    # optional visualization of parameter fields 
    sanity_plotting = FALSE
){
  
  # By default we expect kappa2/theta/rho to be supplied directly
  # If any are missing, fall back to reading and decoding from an HDF5 file
  if (is.null(kappa2) || is.null(theta) || is.null(rho)) {
    
    if (is.null(file_path) || is.null(dataset_name) || is.null(index)) {
      stop(
        "Supply kappa2, theta, and rho directly, or provide file_path, ",
        "dataset_name, and index to read them from an HDF5 file.",
        call. = FALSE
      )
    }
    
    # Optional: rhdf5 is a Bioconductor package needed for loading h5 files
    # Check if installed at runtime, suggest install if not
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
      stop(
        "Reading from an HDF5 file requires the 'rhdf5' package.\n",
        "Install it with:\n",
        "  install.packages(\"BiocManager\")\n",
        "  BiocManager::install(\"rhdf5\")\n",
        "Alternatively, supply kappa2, theta, and rho directly."
      )
    }
    
    params <- rhdf5::h5read(file_path, dataset_name)
    rhdf5::H5close()
    
    # Flipping axis when moving from python to R
    params <- params[, dim(params)[2]:1, , index]
    
    # Recover params
    kappa2 <- exp(params[,,1])    # kappa needs to be scaled
    theta  <- params[,,2] + pi/2  # theta needs to be rotated
    rho    <- params[,,3]         # rho is fine as is
  }
  
  if (sanity_plotting){
    par(mfrow = c(1,3))
    imagePlot(as.surface(gridlist, kappa2),       main = "kappa2",     col = viridis(256))
    imagePlot(as.surface(gridlist, theta - pi/2), main = "theta(adj)", col = viridis(256))
    imagePlot(as.surface(gridlist, rho),          main = "rho",        col = viridis(256))
    par(mfrow = c(1,1))
  }
  
  # need these for encoding into LK
  rhox <- sqrt(rho)
  rhoy <- 1 / rhox
  
  # create H tensor out of params
  H11 <- (rhox^2 * (cos(theta))^2) + (rhoy^2 * (sin(theta))^2)
  H12 <- (rhoy^2 - rhox^2) * (sin(theta) * cos(theta))
  H21 <- H12 
  H22 <- (rhox^2 * (sin(theta))^2) + (rhoy^2 * (cos(theta))^2)
  
  rows <- length(gridlist$x)
  
  # fill the high dimensional stencil (9 fields)
  stencil_tensor <- array(NA, c(rows, rows, 9))
  stencil_tensor[,,1] <- 0.5 * H12
  stencil_tensor[,,2] <- -H22
  stencil_tensor[,,3] <- -0.5 * H12
  stencil_tensor[,,4] <- -H11
  stencil_tensor[,,5] <- kappa2 + 2 * H11 + 2 * H22
  stencil_tensor[,,6] <- -H11
  stencil_tensor[,,7] <- -0.5 * H12
  stencil_tensor[,,8] <- -H22
  stencil_tensor[,,9] <- 0.5 * H12
  
  # put everything into awght obj of a particular class
  awght_obj <- list(x = gridlist$x, y = gridlist$y, z = stencil_tensor)
  class(awght_obj) <- "multivariateSurfaceGrid"
  
  return(awght_obj)
}
