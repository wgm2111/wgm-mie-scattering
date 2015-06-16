
#
# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2010
# liscence: BSD style
#
# ================================================================================

"""

This file is used to interface with a modified version of Michael Mishchenko's 
MIE scattering code.  A copy of the original code is available on Michael's web 
page, ad in the file spher.f 

This module imports from spher.so, which is a compiled python module made from 
spher_.f (an edited version of Michael's code).

** ToDo **
    Reuse repeated computations of "mono" when advantageous.  




# AA=0.6D0
# BB=0.2 D0
# AA1=0.17D0
# AA2=3.44D0
# BB1=DLOG(1.96D0)*DLOG(1.96D0)
# BB2=DLOG(2.37D0)*DLOG(2.37D0)
# GAM=1D0                                                               
# LAM=0.63D0                                            
# MRR=1.53 D0                                               
# MRI=0.008 D0                                            
# NDISTR=3
# NK=100
# N=100
# NP=4
# R1=0 D0
# R2=35D0       
# DDELT=1D-7


  
"""

# standard python imports
from __future__ import print_function
from copy import copy
# from my_tools.parallel import parmap

# third party imports 
import scipy as sp
from scipy import sparse
import scipy.linalg as la

# my imports 
import spher
from single_scattering.core import ss_array, SingleScatteringArray

# Global variables
# ==============================================================================



# maximum number of GSF coefficients (set in spher_.f)
GSFMAX = SingleScatteringArray.GSFMAX                 

# default variables when not supplied explicitly
default_kwds = dict(
    aa = 0.6, bb = 0.2, aa1 = 0.17, aa2 = 3.44, 
    bb1 = sp.log(1.96)**2, bb2 = sp.log(2.37)**2, 
    gam = 1.0,
    area = 0.0,
    volume = 0.0,
    rvw = 0.0,
    lam = 0.63,
    mrr = 1.53,
    mri = .008,
    ndistr = 3, nk = 10, n = 2000, np = 4, r1 = 0.0, r2 = 35.0, ddelt=1e-10)


def _add_default_kwds_doc(fun):
    """
    Append documentation giving the default list of kwd values for spher.
    """
    newdoc = "\nThe default values for sphere are as follows"
    for key, val in default_kwds.iteritems():
        newdoc += "{0} = {1}\n".format(key, val) 
    fun.__doc__ += newdoc
    return fun


# helper functions to define an interface
@_add_default_kwds_doc
def build_args(verbose=False, **kwds):
    """
    Build input args for fortran subroutine "spher" in the file "spher_.f".

    input args for spher routine
    =============================
    ndistr  -  type of distribution to average over (1,2,3,4,5)
    aa      -  first parameter for distribution averaging
    bb      -  second parameter for distribution averaging 
    gam     -  third parameter for the distribution (if needed)
    r1      -  minimum radius
    r2      -  maximum radius
    mrr     -  real refractive index
    mri     -  imaginary refractive index
    lam     -  wavelength of light (in backround media) 

    n       -  number of subintervals for the integration 
    nk      -  gaussian quadrature points used during averaging
    
    """
    # define the list of looked up parameters for spher
    argins = ('aa', 'bb', 'gam', 'lam', 'mrr', 'mri', 'r1', 'r2', 'n', 'np', 
              'ndistr', 'nk', 'area', 'volume', 'rvw', 'rmean', 'aa1', 'bb1', 
              'aa2', 'bb2', 'ddelt')
    
    # copy standard update
    look = copy(default_kwds)
    look.update(kwds) 
    if verbose: print(look)
    
    out = tuple([look.get(argin, 0.0) for argin in argins])
    return out



# ============================================================= #
# Post processing of spher output into dense array structure    #
# ============================================================= #

def output_wrapper(lmax_, a1, a2, a3, a4, b1, b2, cext, csca, lmax=None):
    """
    Wrap the GSF coefficients into a sparse matrix (after un-normalizing) 
    and return lmax, cext, and gsf.
    """
    # optional lmax
    lmax = lmax if lmax else lmax_

    # Check that cext and csca are sane
    if round(cext,10) >= round(csca,10): pass
    else: 
        msg = "expected cext >= csca, but cext={0} and csca={1}"
        raise ValueError(msg.format(cext, csca))
    assert csca >= 0.0
    
    # make a sparse array in lil format
    gsf_data = sp.zeros((6, lmax))
    gsf_data[0,:lmax] = csca * a1[:lmax]
    gsf_data[1,:lmax] = csca * a2[:lmax]
    gsf_data[2,:lmax] = csca * a3[:lmax]
    gsf_data[3,:lmax] = csca * a4[:lmax]
    gsf_data[4,:lmax] = csca * b1[:lmax]
    gsf_data[5,:lmax] = csca * b2[:lmax]

    return ss_array(cext, gsf_data)# SingleScatteringDense(gsf_data)





# function and class definitions 
# ==============================================================================

def mono_dense(lam, mrr, mri, r, **kwds):
    """
    Calculate the scattering matrix for mono-disperse spheres of radius r.
    """
    # eps = 1e-13
    scale_r1 = 0.9999999
    scale_r2 = 1.0000001
    args = build_args(aa=r, bb=1e-1, lam=lam, mrr=mrr, mri=mri, ndistr=4, nk=1, n=1, 
                      r1=r*scale_r1, r2=r*scale_r2)
    # with silence: out = spher.spher(*args)
    out = spher.spher(*args)

    return output_wrapper(*out)


def mono_dense_array(lam, mrr, mri, r_array, lmax, **kwds):
    """
    Calculate the scattering gsf coeficients for mono-disperse spheres
    or radius r.
    """
    assert r_array.ndim == 1
    out_shape = (6, lmax, r_array.size)
    out = sp.zeros(out_shape)

    
    # call to define blocks in parallel
    # if False:                    # parallel

    #     # define a function to pass to parmap
    #     def f(r):
    #         # initialize output
    #         out = sp.zeros(out_shape[:-1])
    #         # call mono and 
    #         _out = mono_dense(lam, mrr, mri, r,**kwds)
    #         _lmax = _out.shape[1]
    #         l = min(_lmax, lmax)
    #         out[:,:l] = _out[:,:l]
    #         return out
    #     # call parmap
    #     out_list = parmap(f, list(r_array))
    #     # store in NOT parallel
    #     for i, _out in enumerate(out_list):
    #         out[:,:,i] = _out[:,:]
    # else:
    
    for i, r in enumerate(r_array):
        _out = mono(lam, mrr, mri, r, lmax=lmax, **kwds)
        _lmax = _out.shape[1]
        l = min(_lmax, lmax)
        out[:,:l,i] = _out[:,:l]

    # md_array = sp.array([mono(lam, mrr, mri, r, lmax=lmax, **kwds)
    #                      for r in r_array]).transpose([1,2,0])
    return SingleScatteringArray(out)

def mono(lam, mrr, mri, r, lmax=None, **kwds):
    """ 
    Define a general purpose mono-disperse scattering function 
    """
    # call mono_dense if r is a scalar
    if sp.isscalar(r):
        out = mono_dense(lam, mrr, mri, r)
        if ~(lmax is None):
            out.reshape_lmax(lmax)
        return out
    # assume r is a vector and call mono_dense_array
    else:
        # Check that lmax is specified for passing a vector of radii
        if lmax is None:
            raise ValueError("Specify lmax when passing a vector for 'r'")
        return mono_dense_array(lam, mrr, mri, r, lmax, **kwds)



# Define an interface to lognormal averaging [spher.f, Mishchenko]
# -----------------------------------------------------------------

def add_lognorm_doc(fun):
    doc = """
Note that the effective parameters reff and veff are equivalent to 
those used by [Mishchenko].

log(veff+1) = log(sigma_g)**2
reff = r_g * exp(5/2 * log(sigma_g)**2)

The use of effective radius gives a better notion of the particle radii
contributing to the scattering.
"""
    fun.__doc__+= doc
    return fun


def lognorm_dense(lam, mrr, mri, reff, veff, rmin=None, rmax=None,verbose=False, **kwds):
    """
    Calculate the scattering matrix for spheres using a lognormal pdf for 
    the size dependence.
    """
    # set the integration bounds
    if rmin: 
        kwds['r1'] = rmin
    else:
        if 'r1' not in kwds: kwds['r1'] = 1e-20
    if rmax: 
        kwds['r2'] = rmax
    else:
        if 'r2' not in kwds: kwds['r2'] = reff + 5 * veff        

    # define incoming args
    bb = sp.log(1+veff)
    aa = reff * sp.exp(-5.0 / 2.0 * bb)
    args = build_args(verbose=verbose, 
                      aa = aa, bb = bb,
                      lam = lam, mrr = mrr, mri = mri, ndistr=2, **kwds)
    
    # with silence: out = spher.spher(*args)
    out = spher.spher(*args)
    lmax = kwds['lmax'] if kwds.has_key('lmax') else None
    return output_wrapper(*out, lmax=lmax)




def lognorm_multilam(lams, mrrs, mris, reff, veff, 
                     rmin=None, rmax=None, lmax=None, 
                     verbose=False, **kwds):
    """
    Calculate the scattering matrix for spheres using a lognormal pdf for 
    the size dependence.    
    """
    assert (lams.size == mrrs.size) and (lams.size == mris.size)

    if not lmax:
        lam = lams[0]
        mrr = mrrs[0]
        mri = mris[0]
    
        test = lognorm_dense(lam, mrr, mri, reff, veff, rmin, rmax, 
                             verbose=verbose, **kwds)
        lmax = test.shape[1]

    # compute the output
    # -------------------
    out = [lognorm_dense(lam, mrr, mri, reff, veff, rmin, rmax, lmax=lmax)
           for lam, mrr, mri in zip(lams, mrrs, mris)]

    out = sp.array(out).transpose([1,2,0])
    
    return SingleScatteringArray(out)


@add_lognorm_doc
def lognorm(lam, mrr, mri, reff, veff, 
            rmin=None, rmax=None, lmax=None, verbose=False, **kwds):
    """
    Calculate the scattering tensor for lognormally distributed sphereical 
    particles.  

    ** input **
        lam       - number or array (lam, mrr, mri same shape)
        mrr       - number or array > 1
        mri       - number or array > 0
        reff      - number > 0
        veff      - number > 0
        rmin      - number > 0
        rmax      - number > rmin

    ** output **
        out       - the scattering tensor (6, lmax, len(lam))

    ** other **
       This routines always returns a DenseSingleScattering object
    
    """

    if sp.isscalar(lam):
        out = lognorm_dense(
            lam, mrr, mri, reff, veff, 
            rmin=rmin, rmax=rmax, 
            verbose=verbose, **kwds)
    elif hasattr(lam, 'size'):
        if lam.size > 1:
            out = lognorm_multilam(
                lam, mrr, mri, reff, veff, 
                rmin=rmin, rmax=rmax, lmax=lmax, 
                verbose=verbose,**kwds)
        else:
            out = lognorm_dense(
                lam[0], mrr, mri, reff, veff, 
                rmin=rmin, rmax=rmax,
                verbose=verbose, **kwds)
    else: 
        raise ValueError('unexpected value passed as "lam"')
    
    return SingleScatteringArray(out)





# ====
# Example script
# ==============================================================================

if __name__=="__main__":

    # imports

    #script
    pass

