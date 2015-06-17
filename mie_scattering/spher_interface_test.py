

#
# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2010
# liscence: BSD style
#
# ==============================================================================

"""
This example script plots the small-scale mie resonances that will complicate 
the integration over size bins. 

"""


# standard python imports
from __future__ import print_function

# third party imports 
import scipy as sp
import numpy.testing as nptest

# my imports 
import spher_interface as sph
mono, mono_dense_array, lognorm = (sph.mono, 
                                   sph.mono_dense_array,
                                   sph.lognorm)



# Global variables
# ==============================================================================
lam = 1.6
mrr = 1.5
mri = 8.0e-5
rfun = lambda x: lam * x /  (2 * sp.pi)


# angles sampling
angles = sp.linspace(0,sp.pi, 90)

def mono_test():
    "Call mono and check for sanity."
    lam = .41
    mr = 1.33
    mi = .0005
    radius = .15
    out = mono(lam, mr, mi, radius)
    assert out.shape[0] == 6
    
def mono_dense_array_test():
    "Call mono and check for sanity."
    lmax = 50
    lam = .41
    mr = 1.33
    mi = .0005
    radius = .1 * sp.rand(5)
    out = mono_dense_array(lam, mr, mi, radius, lmax)
    assert out.shape == (6, lmax, 5)
    # assert out.shape[2:] == radius.shape


def lognorm_test():
    "Call lognorm both ways and check for sanity."
    lams = 1.0 + sp.rand(3) 
    mr = 1.33 + sp.rand(3)
    mi = .0005 + sp.rand(3)
    reff = .01 
    veff = .005
    out = lognorm(lams[0], mr, mi, reff, veff)
    assert out.shape[0] == 6
    out = lognorm(lams, mr, mi, reff, veff)
    assert out.shape[2:] == (3,)


# # make calculations
# lmax = 35
# lams = sp.array([.41, .55])
# mrs = sp.array(2*[1.33, 1.55, 1.58])
# mis = sp.array(2*[0.0008, 0.015, 0.03, 0.06])
# rs = sp.array(2*[1,1.2,1.3,1.4,1.5])
# veff = .001 
# shape = (6, lmax, lams.size, mrs.size, mis.size, rs.size)
# # out = sph.SingleScatteringDense(sp.zeros(shape))




# for i, lam in enumerate(lams):
#     for j, mr in enumerate(mrs):
#         for k, mi in enumerate(mis):
#             for l, r in enumerate(rs):
#                 o = sph.mono_dense(lam, mr, mi, r).reshape_lmax(lmax)
#                 out[:, :, i, j, k, l] = o

# # Testing the shape of the new array 
# def test_SSD_scattering_tensor():
#     "SSD.scattering_tenor new==old"
#     anew = out.scattering_tensor(angles)
#     aold = out.scattering_tensor_old(angles)
#     aold = aold.transpose(0,1,3,4,5,2)
#     nptest.assert_array_almost_equal(anew, aold)
    

# def test_SSD_scattering_tensor_shape():
#     "SSD.scattering tensor reshapes dim 1 only"
#     assert out.shape[2:] == out.scattering_tensor(angles).shape[2:]
    

# out = sp.array([
#         mono(lam,mr,mi,r) 
#         for lam in lams
#         for mr in mrs
#         for mi in mis
#         for r in rs])
# print(shape)
# print(out.shape)

# function and class definitions 
# ===============================================================================



# ====
# Example script
# ===============================================================================

