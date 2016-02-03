
# Summary of wgm-mie-scattering
This Python package is meant to be a simple utility for computing
the scattering properties of airborne spherical particles.  This is
useful for remote sensing of clouds (and airborne particles) in the 
atmosphere.

This code wraps mie scattering codes written by Michael Mishchenko,
and provides documented routines for gamma, and lognormal size averages. 
Also, argument checking is provided.

The output is given as an array of generalized spherical function
coeficients (scaled by the extinction).  This array has a method
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


## You can also install to import `mie_scattering` from any directory. 

For a system install use:
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

