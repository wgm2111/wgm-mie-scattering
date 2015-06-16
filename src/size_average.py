
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
from scipy.integrate import quad
from scipy.optimize import brentq
import tables 

# my imports 
import spher_interface as sph

# Global variables
# ================================================================================



# function and class definitions 
# ================================================================================



# subdivide a interval by equal sub integrals
# -----------------------------------------------------------------------------
def subdivide(g, rmin, rmax, m):
    """
    Subdivide the given interval by equal integrals of g.

    **example**
     >>> g = lambda r : r**2
     >>> rmin = 0.0
     >>> rmax = 1.0
     >>> m = 2
     >>> r = subdivide(g, rmin, rmax, m)
     >>> abs(r - array([0.0, .5**(1./3.), 1.0])).sum() <= 1e-4
     True
    """
    r = sp.linspace(rmin, rmax, m+1)

    # total integral value
    tot = quad(g, rmin, rmax)[0]
    cuts = sp.linspace(0.0, tot, m+1)[1:-1]
    
    # define intervals r_n to r_n+1 to give equal area under g(r)
    r[1:-1] = [brentq((lambda r_: quad(g, rmin, r_)[0]-cut), rmin, rmax,
                      xtol=1.0e-4)
                for cut in cuts]
    # return the sequence of subinterval boundaries
    return r



def vromb(fun, xmin, xmax, 
          relerr=1e-6, abserr=None, basemesh=6, countmax=8,
          verbosity=0):
    """ 
    Integrate the vector valued function 'f' with Romberg adaptive integration 
    on a uniform mesh.  Stop when there is convergence to the specified number 
    of places.
    """
    # error convenience functions
    def L1_relerr(a, b):
        return abs(b-a).sum() / abs(a).sum()
    def L1_abserr(a, b):
        return abs(b-a).mean()

    # initialize told and tnew
    tnew = [None]
    
    # starting grid interval count
    n_new = 2**basemesh

    # calculate trapazoidal integral on base mesh
    x = sp.linspace(xmin, xmax, n_new+1)
    fvals = fun(x)
    tnew[0] = (x[1]-x[0])/2.0 * (fvals[...,:-1]+fvals[...,1:]).sum(-1) # do trap
    # print("x.shape i {0} and fvals.shape is {1}".format(x.shape, fvals.shape))

    # loop until convergence or above count limit
    count = 0
    converged = False

    for count in range(1, countmax+1):
        
        # initialize new data
        n_new, n_old = 2**(basemesh+count), n_new
        told, tnew = tnew, [None]*(count+1)

        # Refined trapz using midpoint =>  TRAP_n+1 = .5 * (TRAP_n + MID_n)
        x = sp.linspace(xmin, xmax, n_old+1)
        xmid = 1.0 / 2.0 * (x[:-1] + x[1:])
        wmid = (xmid[1]-xmid[0]) 
        fsum = fun(xmid).sum(-1)
        # print("told.shape is {0}, wmid.shape is {1}, and fsum.shape is {2}".
              # format(told[0].shape, wmid.shape, fsum.shape))
        tnew[0] = 1./2. * (told[0] + wmid * fsum)

        # extrapolate 
        for k in range(count):
            c = 4**(k+1)
            tnew[k+1] = 1.0 / (c - 1. ) * (c * tnew[k] - told[k])

        # estimate errors
        L1relerr = L1_relerr(tnew[-1], told[-1]) 
        L1abserr = L1_abserr(tnew[-1], told[-1]) 
        
        # check for convergence
        relok = L1relerr < (relerr)
        absok = (L1abserr < (abserr)) if abserr else True
        converged = relok and absok
        
        if converged: 
            break

        # print a progres report for curent subinterval
        # ----------------------------------------------
        if verbosity >=2: 
            if (count in [1, 3, 5, 7]):
                template = """
--
 count = {0},   mesh_size = {1},   dr = {2}

 L1relerr = {3}
 L1abserr = {4}
"""
                dr = (xmax - xmin) / n_new
                args = (count, n_new, dr, L1relerr, L1abserr)
                print(template.format(*args))
            

    out = tnew[-1]

    # print final results for current interval
    # ----------------------------------------------
    if verbosity >= 1:
        template = """

The final error on the interval rmin={0}, rmax={1} is 
 
 relerr    : {2}

 abserr    : {3}

 mesh_size : {4}   with dr={5}

"""
        relerr_out = L1relerr
        abserr_out = L1abserr
        n_out = n_new
        dr_out = (xmax - xmin) / float(n_out)
        args = (xmin, xmax, relerr_out, abserr_out, n_out, dr_out)
        print(template.format(*args))

    return out, sp.array([L1abserr])



# ====
# averaging routines for single scattering
# ================================================================================


# average scattering due to a rv
# -----------------------------------------------------------------------------
def mieave(lam, mrr, mri, rv, lmax=None):
    " average scattering for a given pdf "

    # Setup the integration       
    rmin, rmax = rv.rmin, rv.rmax
    ss = sph.mono_dense(lam, mrr, mri, rmax)
    lmax = lmax if lmax else ss.lmax()

    # define the function to integrate
    def fun(r, weights=None):
        " Compute an array of numbers related to the SingleScattering "
        out = sph.mono_dense_array(lam, mrr, mri, r, lmax=lmax)
        return out * rv.ndf(r)

    # subdivide the integral
    dx_total = 2*sp.pi*(rmax-rmin) / lam
    nint = max(1, dx_total // 10)
    rstart = subdivide(lambda r : sp.pi * r**2 * rv.ndf(r), rmin, rmax, nint)
    rend = rstart[1:]
    
    # loop over intervals to compute the total
    # -----------------------------------------
    out = vromb(fun, rmin, rmax)
    tot = out[0]
    err = out[1]
    count = 0.0
    for a, b in zip(rstart, rend):          
        if count == 0.0:
            out = vromb(fun, a, b)
            tot = out[0]
            err = out[1]
        else:
            out = vromb(fun, a, b)
            tot += out[0]
            err += out[1]
        count+=1
        
    # make a SingleScattering object out of the array 
    # ss = ss._undensify(tot)
    return tot



# ====
# Example script
# ================================================================================

if __name__=="__main__":
    pass
    # # imports
    # import size_parameter as sip 
    # from scatmo.src_fortran import plots
    # lam = .410
    # mrr = 1.5
    # mri = .008
    # reff = 1.5
    # veff = .01
    # rmin = 1e-14
    # rmax = 4.0

    # ss_mike = sph.lognorm_dense(
    #     lam, mrr, mri, reff, veff, rmin=rmin, rmax=rmax, ddelt=1e-10)    
    # la = sip.LognormAerosol(reff, veff, rmin, rmax)
    # ss = mieave(lam, mrr, mri, la, lmax = ss_mike.lmax())
    # # ss_trunk = ss.sparsify(lmax=ss_mike.lmax)
    # # ss = ss.sparsify()

    # gsf = ss.gsf(ss_mike.lmax())
    # gsfmike = ss_mike.gsf()

    
    # print("l1 relative Error is {0}".format(
    #         abs(gsfmike-gsf).sum() / abs(gsfmike).sum()))

    # plots.plot_compare(ss, ss_mike)
    # # plots.plot_single_scattering(ss, holdon=False, label='my calc')
    # # # plots.plot_single_scattering(ss_trunk, label="my truncated")
    # # plots.plot_single_scattering(ss_mike, label="mike's calc")
    
    # err = ss - ss_mike

    # plots.plot_single_scattering(err, label="err", fnumber=1)
    # #script
    # pass

