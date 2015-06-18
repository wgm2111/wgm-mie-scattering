
#
# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2010
# liscence: BSD style
#
# ================================================================================


# standard python imports
from __future__ import print_function

# third party imports 
import scipy as sp
import nose.tools as nt
import numpy.testing as npt

# my imports 

from special import wignerd, gensph

# Global variables
# ================================================================================

lmax = 20

mns = [(0,0), (2,2), (2,-2), (0,2)]
mns += [(0,1), (3,0), (3,5), (7,1)]
theta = sp.pi * sp.rand(2,3)
ae = lambda x,y: npt.assert_array_almost_equal(x,y, decimal=14)


d00 = sp.array([sp.ones(theta.shape), 
                sp.cos(theta), 
                .5 * (3 * sp.cos(theta)**2-1.0),
                -.5 * sp.cos(theta)*(3.0-5.0*sp.cos(theta)**2), 
                1.0/8.0 * (3.0 - 30*sp.cos(theta)**2 + 35*sp.cos(theta)**4)])
d22 = sp.array([sp.zeros(theta.shape), 
                sp.zeros(theta.shape), 
                .25 *  (1.0+sp.cos(theta))**2,
                -.25 * (1.0+sp.cos(theta))**2 * (2-3*sp.cos(theta)), 
                .25 *  (1.0+sp.cos(theta))**2 * (1-7*sp.cos(theta)+7*sp.cos(theta)**2)])
d2m2 = sp.array([sp.zeros(theta.shape), 
                 sp.zeros(theta.shape), 
                 .25 * (1.0-sp.cos(theta))**2,
                 .25 * (1.0-sp.cos(theta))**2 * (2+3*sp.cos(theta)),
                 .25 * (1.0-sp.cos(theta))**2 * (1+7*sp.cos(theta)+7*sp.cos(theta)**2)])
d02 =  sp.array([sp.zeros(theta.shape), 
                 sp.zeros(theta.shape), 
                 .5 * sp.sqrt(3.0/2.0) * sp.sin(theta)**2,
                 1. * sp.sqrt(15./8.0) * sp.sin(theta)**2 * sp.cos(theta),
                 -sp.sqrt(5./32.) * sp.sin(theta)**2 * (1-7*sp.cos(theta)**2)])


# function and class definitions 
# ================================================================================

def test_wignerd_values():
    """
    use analytic expressions to test the wignerd function values
    """
    dlist = [d00, d22, d2m2, d02]
    mlist = [0, 2, 2, 0]
    nlist = [0, 2, -2, 2]
    
    for d, m, n in zip(dlist, mlist, nlist):
        dpy = wignerd(theta, m, n, 5)
        ae(dpy, d.transpose(1,2,0))
        

def test_wignerd_symetry():
    """
    test that wignerd has the right symetry relations 
    """
    
    for m, n in mns:
        
        # first set
        dMN = wignerd(theta, m, n, lmax)
        dNM = wignerd(theta, n, m, lmax)
        dmn = wignerd(theta, -m, -n, lmax)
        dnm = wignerd(theta, -n, -m, lmax)
        ae(dMN, (-1)**(m-n)*dmn)
        ae(dMN, (-1)**(m-n)*dNM)
        ae(dMN, dnm)        

        # second set
        dmN = wignerd(theta, -m, n, lmax)
        dMn = wignerd(theta, m, -n, lmax)        
        _dMN = wignerd(sp.pi - theta, m, n, lmax)
        
        ls = sp.arange(lmax)
        print(ls)
        aux = (-1)**(ls-n)*dmN +_dMN
        print(abs(aux))
        ae(_dMN, (-1)**(ls-n)*dmN)# (-1)**(ls-n)*dmN)
        ae(_dMN, (-1)**(ls+m)*dMn)

def test_gensph():
    """
    Check that wingerd and gensph agree
    """
    for m, n in mns:
        dMN = wignerd(theta, m, n, lmax)
        pMN = gensph(sp.cos(theta), m, n, lmax)
        ae(pMN, 1j**(m-n)*dMN)
        






# ====
# Example script
# ================================================================================

if __name__=="__main__":

    # imports

    #script
    pass

