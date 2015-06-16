
# 
# author: William G.K. Martin 
# email: willful_liam@gmail.com 
# copyright (c) 2011
# liscence: BSD style
# 

"""
# Summary of wgm-mie-scattering
This Python package is meant to be a simple utility for computing
the scattering properties of airborne spherical particles.  This is
useful for remote sensing of clouds (and airborne particles) in the 
atmosphere. 

## Install (Not working yet)
For a system install use:
  $ python setup.py install
For building in the local directory 
  $ python setup.py build

## Requirements
Python with numpy, scipy, f2py, and some other common packages.  The
Anaconda python distribution is good, and will have all the necessary 
third party packages. 

## Usage
A non-exhaustive list of ways to do single-scattering calculations 
are as follows.  See function documentation for details.

  >>> from miescattering import mono, mieave, lognorm
  >>> single_scattering_array = mono(lam, reff, veff, mr, mi)
  >>> single_scattering_array = lognorm(lam, reff, veff, mr, mi)
  >>> single_scattering_array = mieave(
        pdf_of_radius, lam, veff, mr, mi)

The 'single_scattering_array' objects allow for calculation of 
scattering at arbitrary scattering angles: 
  >>> single_scattering_array.at_angles(theta). 

See also, ssa.atensor(theta).  The single_scattering_array is a 
subclass of "numpy.ndarray".  Use its extra methods to transform 
from GSF coefficients to specific scattering-angles.

## References
Fortran source code thanks to Michael Mishchenko: 
http://www.giss.nasa.gov/staff/mmishchenko/

Python wrapping thanks to Will Martin:
https://github.com/wgm2111/wgm-mie-scattering

"""


import size_average

# Routines from Michael's Fortran code
import spher_interface
from spher_interface import mono, lognorm # mono and lognormal averages
import spher_interface_test               # some test routines

# Averaging routines in Python
import size_average

# Extra functions
import z_scale_mie




# from size_average import mieave
# from spher_interface import lognorm, mono
