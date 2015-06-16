#!/usr/bin/env python

# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2011
# liscence: BSD style


# future
from __future__ import division, absolute_import, print_function


# Import the numpy version of distutils (which has support for f2py)
import setuptools
from numpy.distutils.core import Extension
from os.path import join



# Global names
# ------------------
NAME = "miescattering"
VERSION = '0.0.1'
SOURCE_DIR = 'src'#join(PACKAGE_DIR, 'src')
TARGET_SO_LIB = ''
PACKAGES = ['']
PACKAGE_DIR = {'':'src'}



def package_path(fname, package_dir=PACKAGE_DIR):
    "Join filename with the PACKAGE_DIR"
    return join(package_dir, fname)
def source_path(fname, package_dir=PACKAGE_DIR):
    "Join filename with the PACKAGE_DIR"
    return join(package_dir, fname)


# Packages

# PACKAGES = setuptools.find_packages()#['src']

# Files
ext_modules = [Extension(TARGET_SO_LIB+'.'+'wig', 
                         sources = [join(SOURCE_DIR, 'wig.f95'), 
                                    join(SOURCE_DIR, 'fact.f95')]),
               Extension(TARGET_SO_LIB+'.'+'fact', 
                         sources = [join(SOURCE_DIR, 'fact.f95')]),
               Extension(TARGET_SO_LIB+'.'+'lacis', 
                         sources = [join(SOURCE_DIR, 'lacis.f')]), 
               Extension(TARGET_SO_LIB+'.'+'matr', 
                         sources = [join(SOURCE_DIR, 'matr_.f')]), 
               Extension(TARGET_SO_LIB+'.'+'spher', 
                         sources = [join(SOURCE_DIR, 'spher_.f')])]



# Define the configuration function 
def configuration(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(NAME, parent_package, top_path)
    config.make_config_py()
    return config


# MIE_SOURCE_FILES = ['spher_.f', 'matr_.f']
# LACIS_SOURCE_FILES = [wig.f95]


# # Extension to call 
# ext_spher = Extension(
#     name = 'spher',
#     sources = [os.path.join(SOURCE_DIR, fname) 
#                for fname in MIE_SOURCE_FILES])





if __name__ == "__main__":
    
    # import setup
    from numpy.distutils.core import setup
    

    # ------------------------------------------------------------------------------
    setup(name = NAME,
          version = VERSION,
          description = "Single scattering computations of airborne particles.",
          author = "William GK Martin",
          author_email = "willful_liam@gmail.com",
          packages=PACKAGES, 
          package_dir=PACKAGE_DIR,
          configurateion=configuration(),
          license = 'GPL2', 
          ext_modules=ext_modules)


