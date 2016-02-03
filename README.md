
# Summary of wgm-mie-scattering
This Python package is meant to be a simple utility for computing
the scattering properties of airborne spherical particles.  This is
useful for remote sensing of clouds (and airborne particles) in the 
atmosphere. 

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


## Install (Not working yet)

For a system install use:
```
  $ python setup.py install
```

For developing in place:
```
  $ python setup.py develop
```

For building in the local directory:
```
  $ python setup.py build
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

