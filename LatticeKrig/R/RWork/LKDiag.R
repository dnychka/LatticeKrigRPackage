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

LKDiag <- function(entries, nrow, diags = NULL, ncol = nrow, full = FALSE) {
  entries = as.double(entries)
  nEntries = as.integer(length(entries))
  nrow = as.integer(nrow)
  ncol = as.integer(ncol)
  mat <- as.double(matrix(0, nrow=nrow, ncol=ncol))
  
  if (is.null(diags)) {
    if (nEntries %% 2 == 0) {
      diags = c((-nEntries/2) : -1, 1 : (nEntries/2))
    } else {
      diags = (-(nEntries-1)/2) : ((nEntries-1)/2)
    }
  }
  
  diags = as.integer(diags)
  if(max(diags) >= ncol || min(diags) <= -nrow) {
    warning("One of the given diagonals is outside the matrix, will be ignored")
    diags = diags[(diags < ncol) && (diags > -nrow)]
  }
  if (length(entries) == 1) {
    entries = as.double(rep(entries[1], length(diags)))
    nEntries = as.integer(length(entries))
  }
  if (length(entries) != length(diags)) {
    stop("The length of entries and length of diagonals don't match")
  }
  if(full) {
    out <- .Fortran("LKDiag", entries = entries, nEntries = nEntries, diags = diags,
                    nRow = nrow, nCol = ncol, matrix = mat, PACKAGE="LatticeKrig")
    return(matrix(out$mat, nrow=nrow, ncol=ncol))
  } else {
    ind <- NULL
    ra <- NULL
    da <- c(nrow, ncol)
    for(idx in 1:length(diags)) {
      diag <- diags[idx];
      entry <- entries[idx];
      if (diag < 0) {
        diagLength <- min(ncol, nrow + diag)
        startRow <- 1 - diag
      } else {
        diagLength <- min(nrow, ncol - diag)
        startRow <- 1
      }
      rowInd <- (1:diagLength) - 1 + startRow
      colInd <- rowInd + diag
      diagInd <- cbind(rowInd, colInd)
      ind <- rbind(ind, diagInd)
      ra <- c(ra, rep(entry, diagLength))
    }
    ord <- order(ind[,1], ind[,2])
    ind <- ind[ord,]
    ra <- ra[ord]
    return(spind2spam(list(ind = ind, da = da, ra = ra)))
  }
}
