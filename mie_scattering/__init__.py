

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

# Routines from Michael's Fortran code
from spher_interface import lognorm

import size_distribution 



# from size_average import mieave
# from spher_interface import lognorm, mono



# import special, wig

# from single_scattering_array import (
#     SingleScatteringArray, ss_array)

# from size_parameter import (
#     ModeBasis, GammaAerosol, LognormAerosol, MultiLognormAerosol)

# from size_parameter_log import (
#     LogModeBasis, LogModeVec)

# from refractive_index import mr_water, mi_water		
# from refractive_index_daimon import mr_water_daimon_21p5
