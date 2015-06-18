#!/usr/bin/env python

# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2011
# liscence: BSD style


# future
from __future__ import division, absolute_import, print_function


# Import the numpy version of distutils (which has support for f2py)
from mie_scattering.spher_interface import lognorm
 
# Parameters for this example
WAVELENGTH = .532               # [microns]
REFF = 0.1                      # effective radius [microns]
VEFF = 0.05                     # a little be of variance
MR = 1.45                       # real refractive index
MI = 0.001                      # imaginary ref. index

# Call lognorm to do the calculation 
ssa_lognorm = lognorm(WAVELENGTH, MR, MI, REFF, VEFF) 


