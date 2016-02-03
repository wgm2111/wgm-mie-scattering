
# Summary of wgm-mie-scattering
This Python package is allows you to compute the scattering properties
of airborne spherical particles.  This is useful for remote sensing of
clouds (and aerosols) in the atmosphere.  Mie calculations are in FORTRAN
codes written by Michael Mishchenko.  This package allows you to run
Michael's calculations from Python.  I have also provided some documentation
and argument checking for these routines.

The output is given as an array of generalized spherical function
coefficients (scaled by the extinction), but the output arrays have a method
called `angle_eval` which returns the scattering matrix at specified
angles. 


## Quick try out

(optional) Install the "Anaconda python distribution" by continuum analytics.

Check out this repository:
```
$ git clone https://github.com/wgm2111/wgm-mie-scattering
```

Change to the repository directory and build FORTRAN source code
```
$ python setup.py distribute
```

Run the example ipython notebook
```
$ ipython notebook 
```


## Install to import from any directory 

For a system wide install use:
```
  $ python setup.py install
```


## Requirements
Python with numpy, scipy, f2py, and some other common packages.  The
Anaconda python distribution is good, and will have all the necessary 
third party packages. 

## References
Fortran source code thanks to Michael Mishchenko: 
http://www.giss.nasa.gov/staff/mmishchenko/

Python wrapping thanks to Will Martin:
https://github.com/wgm2111/wgm-mie-scattering

