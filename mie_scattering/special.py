
# 
# author: William G.K. Martin 
# email: willful_liam@gmail.com 
# copyright (c) 2011
# liscence: BSD style
# 
"""
Special provides routines for computing the wigner-d function and general
spherical function values.  


  theta                   # arbirary shape (or scalar)
  x = cos(theta)         

  # FORTRAN:
  gensph(    x, m, n, lmax)  returns  P_mn^l(x)      for  0 <= l < lmax

  # PYTHON:
  wigner(theta, m, n, lmax)  returns  d_mn^l(theta)  for  0 <= l < lmax

The arrays returns have shape,  

  theta.shape + (lmax,)

** resolved issues **
  The recursions are calculated with a python loop  Could probably be 
  sped up using FORTRAN.  Evaluation took about 15 times longer 
  than pure array multiplication.  Order(10 MFlops).
  
"""


# standard python imports
from __future__ import print_function

# third party imports 
import scipy as sp

# my imports 
from lib import wig

# Global variables
# =============================================================================
wigfact = wig.fact


# function and class definitions 
# =============================================================================

def gensph(x, m, n, lmax):
    """
    call the fortran wigner d routine 

    ** to do **
    
      (1) Add a layer of logic around the FORTRAN call to store previously 
          computed in memory values for the last used set of parameters.   

    """
    if len(x.shape) == 2:
        out = wig.wignerd_cos(x, m, n, lmax)
        coef = (-1)**((m-n)/2) if (m-n)%2 == 0 else 1.0j**(m-n)
        return coef*out
    elif len(x.shape) == 1:
        out = wig.wignerd_cos(x[sp.newaxis], m, n, lmax).squeeze()
        coef = (-1)**((m-n)/2) if (m-n)%2 == 0 else 1.0j**(m-n)
        return coef*out
    else:
        return gensph_python(x, m, n, lmax)



def wignerd(theta, m, n, lmax):
    """ 
    Compute and return the first lmax generalized spherical functions
    """
    #Compute the cosine of the angles
    x = sp.cos(theta)

    # starting l
    lmin = max(abs(m), abs(n))
    xi = 1 if n>=m else (-1)**(m-n)
    
    # initialize wigner-d function value array
    if hasattr(x,'shape'):
        d = sp.zeros(x.shape + (lmax,))
    else:
        d = sp.zeros((lmax,))
    
    # Itteration if m=n=l=0
    if lmin == 0:
        d[..., 0] = 1.0
        d[..., 1] = x
    # loop up to lmax
        for l in sp.arange(lmin+1, lmax-1):
            d[..., l+1] = (
                ((2*l+1.0)*(l*(l+1)*x - m*n) * d[..., l] - 
                 (l+1.0)*sp.sqrt(l**2-m**2)*sp.sqrt(l**2-n**2)*d[..., l-1])
                /(l * sp.sqrt((l+1)**2-m**2) * sp.sqrt((l+1)**2-n**2)))
    # Itteration for all other m, n, l
    else:
    
        d[..., lmin] = xi * 2**(-lmin) * (
            sp.sqrt(wigfact(2*lmin) / 
                    (wigfact(abs(m-n)) * wigfact(abs(m+n))))
            )*(1.0-x)**(abs(m-n)/2.0) * (1.0+x)**(abs(m+n)/2.0)
    
    # loop up to lmax
        for l in sp.arange(lmin, lmax-1):
            d[..., l+1] = (
                ((2*l+1.0)*(l*(l+1)*x - m*n) * d[..., l] - 
                 (l+1.0)*sp.sqrt(l**2-m**2)*sp.sqrt(l**2-n**2)*d[..., l-1])
                /(l * sp.sqrt((l+1)**2-m**2) * sp.sqrt((l+1)**2-n**2)))
    # normalize
    # coef = sp.sqrt(sp.arange(d.shape[-1]) + 1.0/2.0)
    # return d*coef
    return d

def gensph_python(x, m, n, lmax):
    """ 
    compute and return the first lmax generalized spherical functions
    """
    # starting l
    lmin = max(abs(m), abs(n))
    xi = 1 if n>=m else (-1)**(m-n)
    
    # initialize wigner-d function value array
    if hasattr(x,'shape'):
        d = sp.zeros(x.shape + (lmax,))
    else:
        d = sp.zeros((lmax,))
    
    # Itteration if m=n=l=0
    if lmin == 0:
        d[..., 0] = 1.0
        d[..., 1] = x
    # loop up to lmax
        for l in sp.arange(lmin+1, lmax-1):
            d[..., l+1] = (
                ((2*l+1.0)*(l*(l+1)*x - m*n) * d[..., l] - 
                 (l+1.0)*sp.sqrt(l**2-m**2)*sp.sqrt(l**2-n**2)*d[..., l-1])
                /(l * sp.sqrt((l+1)**2-m**2) * sp.sqrt((l+1)**2-n**2)))
    # Itteration for all other m, n, l
    else:
    
        d[..., lmin] = xi * 2**(-lmin) * (
            sp.sqrt(wigfact(2*lmin) / 
                    (wigfact(abs(m-n)) * wigfact(abs(m+n))))
            )*(1.0-x)**(abs(m-n)/2.0) * (1.0+x)**(abs(m+n)/2.0)
    
    # loop up to lmax
        for l in sp.arange(lmin, lmax-1):
            d[..., l+1] = (
                ((2*l+1.0)*(l*(l+1)*x - m*n) * d[..., l] - 
                 (l+1.0)*sp.sqrt(l**2-m**2)*sp.sqrt(l**2-n**2)*d[..., l-1])
                /(l * sp.sqrt((l+1)**2-m**2) * sp.sqrt((l+1)**2-n**2)))
    
    # scale to make a generalized spherical function
    coef = (-1)**((m-n)/2) if (m-n)%2 == 0 else 1.0j**(m-n)
    # coef = coef * sp.sqrt(sp.arange(d.shape[-1]) + 1.0/2.0)
    return coef * d



def atensor(theta, lmax=500):
    """ 
    Compute the tensor for given angles and lmax. 

    ** usage **
      To construct an array for sampling angles to the following.
      
      A = atensor(theta, lmax)
      ssa = SingleScatteringArray(defining_scattering_data)
      [f_i(theta_j)]_ij... = sp.einsum('ijkl, kl...->ij...', A, ssa)
       
    """
    x = sp.cos(theta)
    out = sp.zeros((6, theta.size, 6, lmax))
    out[(0,3),:,(0,3),:] = gensph(x, 0, 0, lmax)
    out[(4,5),:,(4,5),:] = gensph(x, 0, 2, lmax)
    out[1,:,1,:] = gensph(x,2,2,lmax)
    out[2,:,2,:] = gensph(x,2,-2,lmax)

    # pre and post arrays
    a = sp.eye(6)
    a[1:3, 1:3] = .5 * sp.array([[1,1],[1,-1]])
    b = sp.eye(6)
    b[1:3, 1:3] = sp.array([[1,1],[1,-1]])
    
    # make the final array
    # out = sp.einsum('ij,jklm->iklm', a, out) 

    # print(dict(a=a.shape,b=b.shape,out=out.shape))

    out[:] = sp.tensordot(a,out,axes=[1,0])
    


    # next step
    # out = sp.einsum('iklm,lo->ikom', out, b) 
    I = sp.arange(out.ndim)
    out[...] = sp.dot(
        out.transpose((I-1)%I.size), b).transpose((I+1)%I.size)
    return out


# ===
# Example script
# =============================================================================

if __name__=="__main__":

    #script
    theta = sp.linspace(0,sp.pi,300)
    A = atensor(theta)



