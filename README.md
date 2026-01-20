# LatticeKrigRPackage
LatticeKrig R Package: spatial interpolation for large data sets.

LatticeKrig is frozen at version 8.0 in the NCAR repository.

This repository will now be used for current source and development for the 
R package.

# LatticeKrig: Multiresolution Kriging Based on Markov Random Fields

 
This repository contains some supplemental material in the  top level and
 see the subdirectory **LatticeKrig** for the "head" of the standard R package. 
The most current package on CRAN  is listed here as
*LatticeKrig_VERSION.tar.gz* .

To create a possibly new version from this repository download the **LatticeKrig** subdirectory and in UNIX
```
 R CMD build --force LatticeKrig
```
To install this version for your system try
``` 
R CMD INSTALL LatticeKrig
```

or from the tar ball
```
R CMD INSTALL LatticeKrig_VERSION.tar.gz
```
where VERSION are the correct version numbers in this file.

LatticeKrig uses some fortran that will have to be compiled see as an
example for macs

https://cran.r-project.org/bin/macosx/toolsfor


