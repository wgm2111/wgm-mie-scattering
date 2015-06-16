
# 
# author: William G.K. Martin 
# email: willful_liam@gmail.com 
# copyright (c) 2011
# liscence: BSD style
# 

"""
Module for setting up a database computation for interpolation of 
scattering properties based on a scaled z parameter. 

  z = 4*pi*r(m-1) / lam

This parameter has the convenience of aligning the lowest frequency 
oscilations of cross-sectional scattering behavior. 

"""



# standard library 
from __future__ import print_function

# scipy as friends
import scipy as sp
# from scipy import signal
# import scipy.interpolate as interp

# my imports
import single_scattering.mie.core.size_average as size_average
from  spher_interface import mono 
_mono = mono





# define a tent function constructor 
# --
def get_tent(z0, z1):
    " Return a tent function with integral equal to 1 "
    def out(z, z0=z0, z1=z1):
        "Tent function"
        f = 2.0 / (z1-z0) * (1 - abs(2.0 / (z1 - z0) * (z - (z0 + z1)/2.0)))
        return sp.maximum(f, 0.0)
    return out
def get_rect(z0, z1):
    " Return a tent function with integral equal to 1 "
    def out(z, z0=z0, z1=z1):
        "Tent function"
        v = abs(z-(z0+z1)/2.0)
        f = sp.zeros_like(z)
        f[v<((z1-z0)/2.0)] = 2.0 / (z1-z0)
        return f
    return out


def mono(r, lam, a, b, lmax=None):
    " Return the SSA for mono-disperse spheres "
    return _mono(lam, a, b, r, lmax)

def qmono(r, lam, a, b, lmax=None):
    " Return the SSA for mono-disperse spheres "
    out = mono(r, lam, a, b, lmax) / (sp.pi * r**2)
    return out


# Define a class for managing single-refractive Mie scattering calculations
# --
class ZScaleMie(object):
    """
    Class for managing grid construction for a single refractive index 


    Construct a scaling and z grid for the supplied parameters.

    members
    ========
    a       - mrr
    b       - mri
    rmin    - min radius needed
    rmax    - max radius needed
    dr      - approximate resolutin in radius
    wlmin   - min wavelength
    wlmax   - max wavelength
    dwlrel  - approx. resolution in wavelength

    zmin    - min z parameter value
    zmax    - max z parameter value
    dz      - upper bound on absolute z resolution 
    dzrel   - upper bound on relative z resolution 

    """


    # Choose the maximum number of coefficients 
    # --
    LMAX = 500


    def __init__(self, a, b, rmin, rmax, dr, wlmin, wlmax, dwlrel):
        """
        argin
        ========
        a       - mrr
        b       - mri
        rmin    - min radius needed
        rmax    - max radius needed
        dr      - approximate resolutin in radius
        wlmin   - min wavelength
        wlmax   - max wavelength
        dwlrel  - approx. resolution in wavelength

        """

        # Argin
        # --
        self.a, self.b  = a, b
        self.rmin, self.rmax, self.dr = rmin, rmax, dr
        self.wlmin, self.wlmax, self.dwlrel = wlmin, wlmax, dwlrel

        # Prepare a set of znodes
        # --
        self.zmin, self.zmax, self.dz, self.dzrel = self.get_resolution()
        self.zmid, self.Nlo, self.Nhi = self.get_znodes_parameters()
        self.N = self.Nlo + self.Nhi
        self.znodes = self.get_znodes()

    # Size parameter mappings
    # --
    def z(self, r, wavelength=1.0):
        "Return the size parameter"
        out = 4*sp.pi*r*(self.a-1.) / wavelength
        return out 
    def r(self, z, wavelength=1.0):
        "return the radius assuming lam=1"
        out = z * wavelength / (4 * sp.pi * (self.a-1.0))
        return out

    def mono_z(self, z, lmax=None):
        "The scattering of size parameter z refractive index a+ib"
        lmax = self.LMAX if (lmax is None) else lmax
        _r = self.r(z)
        return mono(_r, 1.0, self.a, self.b, lmax)
    def qmono_z(self, z, lmax=None):
        "The scattering efficiency of size parameter z refractive index a+ib"
        lmax = self.LMAX if (lmax is None) else lmax
        _r = self.r(z)
        return qmono(_r, 1.0, self.a, self.b, lmax)


        # methods for computing mono scattering
    # --
    def zq_rect(self):
        " Compute the rect scattering average "
        zstarts = self.znodes[:-1]
        zstops = self.znodes[1:]
        zmids = .5 * (zstarts + zstops)
        out = sp.array([
                1.0 / (z1 - z0) * size_average.vromb(
                    self.qmono_z, z0, z1, 
                    basemesh=3, countmax=9, relerr=1e-7, verbosity=1)[0]
                for z0, z1 in zip(zstarts, zstops)]).transpose(1,2,0)
        return zmids, out
    
    # Resolution calculators
    # --
    def get_resolution(self):
        " Compute the resolution (zmin, zmax, dz, dzrel)."
        zmin = 4 * sp.pi * self.rmin * (self.a-1) / self.wlmax
        zmax = 4 * sp.pi * self.rmax * (self.a-1) / self.wlmin
        dz = 4 * sp.pi * (self.a-1) * self.dr / self.wlmax
        dzrel = abs(self.dwlrel)
        return zmin, zmax, dz, dzrel
    def get_znodes_parameters(self):
        " Compute the zmid and Nlo, and Nhi needed constructing nodes. "
        zmid = self.dz / self.dzrel
        print("zmid = {0}".format(zmid))
        if zmid <= self.zmax:
            Nlo = int(round((sp.log(zmid) - sp.log(self.zmin)) / self.dzrel))
            Nhi = int(round((self.zmax - zmid) / self.dz))
        else:
            zmid = self.zmax
            Nlo = int(round((sp.log(zmid) - sp.log(self.zmin)) / self.dzrel))
            Nhi = 0#int(round((self.zmax - zmid) / self.dz))
        return zmid, Nlo, Nhi


            
    def get_znodes(self):
        " Compute a nodes for a log-lower and linear-upper grid. "
        zlower = sp.exp(sp.linspace(sp.log(self.zmin), sp.log(self.zmid), self.Nlo))
        zupper = sp.linspace(self.zmid, self.zmax, self.Nhi)
        znodes = sp.concatenate([zlower, 
                                 zupper[1:]])
        return znodes




    # overloading
    # --
    def __repr__(self):
        """
        Overloading for helpful print message.
        """
        _repr = []

        out = (
            "\nThe refractive indices are MRR={0} and MRI={1}".format(self.a, self.b) +
            "\n====================================================="+
            "\n\nResolution parameters are set as follows:"
            "\n\tRMIN={0}, RMAX={1}, DR={2}".format(self.rmin, self.rmax, self.dr)+
            "\n\tWLMIN={0}, WLMAX={1}, DWL={2}".format(self.wlmin, self.wlmax, self.dwlrel)+
            "\n\nThese give resolution constraints for the size parameter:"+
            "\n\tzmin={0}, zmax={1}, dz={2}, dzrel={3}".format(
                self.zmin, self.zmax, self.dz, self.dzrel)+
            "\n\nA hybrid log-to-linear set of nodes will require:"+
            "\n\tzmid={0}, Nlo_log={1}, Nhi_lin={2}, Ntotal={3}".format(
                self.zmid, self.Nlo, self.Nhi, self.N))
                
        return out




# Testing script
# --
if __name__ == "__main__":

    from time import time
    

    A = 2.5#1.28
    B = 1e-10
    RMIN = .001                 
    RMAX = 20#50.0
    DR = .010
    WLMIN = .28
    WLMAX = 3.9
    DWLREL = .001#5
    
    zscale = ZScaleMie(A, B, RMIN, RMAX, DR, WLMIN, WLMAX, DWLREL)
    print(zscale)

    # start = time()
    # zmids, ssa = zscale.zq_rect()
    # time_local_means = time() - start
    # print("Computeing local means took {0} seconds".format(round(time_local_means, 2)))
    

