
# 
# author: William G.K. Martin 
# email: willful_liam@gmail.com 
# copyright (c) 2011
# liscence: BSD style
# 



"""

This module defines a numpy.ndarray subclass called SingleScatteringArray.

The point of this array is to store arbitrary single scattering functions 
as generalized spherical functions.  

** examples **

# make a SingleScatteringArray using the constructor ss_array
>>> from scipy import linspace, rand
>>> ext = 1.0
>>> lmax = 10
>>> gsf = rand(6,lmax, 2, 3, 4)
>>> gsf[(1,2,4,5),0:2,...] = 0.0     # these coefs are always zero
>>> ssa = ss_array(ext, gsf)         # makes a SingleScatteringArray 
>>> print(ssa.shape)
(6, 10, 2, 3, 4)

# Check that sampling at angles gives the correct shape 
>>> angles = linspace(0., 1., 13)
>>> ssa_angles = ssa.scattering_tensor(angles)
>>> print(ssa_angles.shape)
(6, 13, 2, 3, 4)

"""

# standard python imports
from __future__ import print_function
from copy import copy

# third party imports 
import scipy as sp
import scipy.linalg as la

# my imports 
from special import gensph, atensor 




# Class definition
# ==============================================================================


# maximum number of GSF coefficients (set in spher_.f)
GSFMAX = 20000                 



class SingleScatteringArray(sp.ndarray):
    """
    A numpy.ndarray subclass for storing single scattering.  For storing the 
    generalized spherical function coefficients used to compute block diagonal
    scattering matrices.

    ** Note **

      These scattering arrays are NOT nessesarily normalized.  Mostly, they
      store quantities with units of [cross-section / volume], but also
      work for [cross-section / particle].  Addition is well defined for the 
      [area/volume] scattering arrays, but not for the [area/particle] arrays. 

      For example, the output from mono-disperse calculations (see also) has 
      units of [area/particle] so to average, we integrate over a number 
      concentration density function with units of [particles/volume].      

      Also, the extinction coefficient is stored in gsf_data[6,0], and may 
      or may not be the extinction crossection. But the single scattering 
      albedo is 
      
        gsf_data[0,0]/gsf_data[6,0].
        

    ** See Also **
      src_fortran.special - for angle calculations.
      single_scattering.mie_mono and single_scattering.mie_lognorm - for 
        mono disperse and lognormally averaged scattering due to spheres
        
    """
    GSFMAX = GSFMAX
    @staticmethod
    def __intest(arr):
        assert (arr.shape[0]==6) #and (len(arr.shape)==2)        
        
    # ----
    # Array subclass method definitions
    # -----------------------------------
    def __new__(cls, gsf_data):
        """
        Construct a single scattering representation.
        """
        # chech that gsf_csr is dense and the correct shape
        cls.__intest(gsf_data)
        out = sp.asarray(gsf_data).view(type=cls)
        return out
    def __array_finalize__(self, obj):
        " finalize array "
        if obj is None: return 
    def __array_wrap__(self, arr_out, context=None):
        "wrap array"
        arr_out = arr_out.view(type=self.__class__)

        
    # def __add__(self, a, b):
    #     " Add SingleScatteringArrays with possible resizing of axis=1 "
    #     almax = a.lmax()
    #     blmax = b.lmax()
        
    #     if almax>=blmax:
    #         out = a.copy()
    #         sp.ndarray.__iadd__(out[:,:blmax,...], b) 
    #         return out
    #     elif almax<blmax:
    #         out = b.copy()
    #         sp.ndarray.__iadd__(out[:,:almax,...], a)
    #         return out

    def _iadd(self, b):
        " Add the smaller SingleScatteringArray to the larger "
        a = self
        almax, blmax = a.lmax(), b.lmax()
        if almax>=blmax:
            sp.ndarray.__iadd__(a[:,:blmax,...], b)
            return a
        elif almax<blmax:
            sp.ndarray.__iadd__(b[:,:almax,...], a)
            return b
       
    # ----
    # Data retrieving methods
    # -------------------------
    
    def lmax(self):
        " Return the maximum non-zero coefficient "
        return self.shape[1]

    def ext(self):
        """ 
        Return the extinction coefficient (not nessesarily the extinction 
        cross-section). 
        """
        return self[5,0].copy()

    def gsf(self, lmax=None):
        " Return a gsf array (6,lmax) "
        out = (self.view(type=sp.ndarray).copy()[:,:lmax]
               if lmax else self.view(type=sp.ndarray).copy())
        out[5,0] = 0.0
        return out

    def sca(self):
        """
        Return the scattering crossection of particle population. 
        """
        return self[0,0].copy()


    # ----
    # Data conversion methods 
    # -------------------------
    def scattering_tensor_depritiated(self, angles):
        """
        Evaluate for the scattering tensor at specifiec angles.

        ** note **
        10/13/2011 : Depriciated.  Use angle_eval
          
        """
        # get the dense generalized spherical function coefficients
        # ------------------------------------------------------------
        gsf = self.gsf()
        gsf[5,0,...] = 0.0

        # calculate the cosines of the angles
        # --------------------------------------
        # print(sp)
        cosines = sp.cos(angles)

        # evaluate the generalized spherical functions using  
        # -----------------------------------------------------
        p00 = gensph(cosines, 0, 0, gsf.shape[1])
        p22 = gensph(cosines, 2, 2, gsf.shape[1])
        p2m2 = gensph(cosines, 2, -2, gsf.shape[1])
        p20 = gensph(cosines, 2, 0, gsf.shape[1])

        # compute the scattering matrix values at each angle
        # -----------------------------------------------------        

        sh= gsf.shape
        out_shape = (6, angles.size) + sh[2:]

        # initialize the output.  This forces the correct shape =)
        out = sp.zeros(out_shape)
        out[0] = sp.dot(gsf[0].T, p00.T).T
        out[1] = sp.dot((gsf[1]+gsf[2]).T, p22.T).T
        out[2] = sp.dot((gsf[1]-gsf[2]).T, p2m2.T).T
        out[1] = .5 * (out[1]+out[2])
        out[2] = out[1] - out[2]
        out[3] = sp.dot(gsf[3].T, p00.T).T
        out[4] = sp.dot(gsf[4].T, p20.T).T
        out[5] = sp.dot(gsf[5].T, p20.T).T
        return out

    
    def atensor(self, theta):
        """
        Use special.atensor and the local lmax information to compute 
        the angular tensor assosiated with the current instance.

        ** Use **
        
          The following two arrays are equivilent:
          (1) self.scattering_tensor(theta) 
          (2) tensordot(self.atensor(theta), self, axes=[(2,3),(0,1)])
          
        """
        return atensor(theta, lmax=self.lmax())


    def angle_eval(self, theta):
        """
        Use atensor to evaluate the scattering matrix elements at vector
        of angles. 
        * in *
          theta - angles to sample scattering matrix elements  [radians]
        """

        # OLD WAYS to compute this
        # 
        # # this is a big temporary array 
        # big = self.atensor(theta) * self
        # return big.sum(3).sum(2)
        # return sp.einsum("ijkl...,kl...->ij...", self.atensor(theta), self)
        # 
        # These were raising SEG FAULTS!!! WHY???
        # 
        
        # Scattering matrix 
        aten = self.atensor(theta)
        P = sp.tensordot(aten, self, axes=[(-2, -1), (0, 1)])
        return P



    def scattering_tensor(self, theta):
        " Old name for routine.  Use the more clear name 'angle_eval'. "
        return self.angle_eval(theta)

    def scattering_tensor_really_old(self, angles):
        """
        Evaluate for the scattering tensor at specifiec angles.
        """
        # get the dense generalized spherical function coefficients
        # ------------------------------------------------------------
        gsf = self.gsf()
        gsf[5,0,...] = 0.0

        # calculate the cosines of the angles
        # --------------------------------------
        cosines = sp.cos(angles)

        # evaluate the generalized spherical functions using  
        # -----------------------------------------------------
        p00 = gensph(cosines, 0, 0, gsf.shape[1])
        p22 = gensph(cosines, 2, 2, gsf.shape[1])
        p2m2 = gensph(cosines, 2, -2, gsf.shape[1])
        p20 = gensph(cosines, 2, 0, gsf.shape[1])

        # compute the scattering matrix values at each angle
        # -----------------------------------------------------    
        sh= gsf.shape
        I = sp.arange(len(sh))
        I[1], I[-1] = I[-1], I[1]
        gsf = gsf.transpose(I)

        f1 = sp.inner(p00, gsf[0,...]).squeeze()

        f2 = sp.inner(p22, gsf[1,...]+gsf[2,...]).squeeze()
        f3 = sp.inner(p2m2, gsf[1,...]-gsf[2,...]).squeeze()
        f2[...] = .5 * (f2 + f3)
        f3[...] = f2 - f3
        f4 = sp.inner(p00, gsf[3,...]).squeeze()
        f5 = sp.inner(p20, gsf[4,...]).squeeze()
        f6 = sp.inner(p20, gsf[5,...]).squeeze()

        # return the scattering tensor 
        # ------------------------------------------
        return sarray(sp.array([f1, f2, f3, f4, f5, f6]))

    def scatmats(self, angles):
        """
        Evaluate the scattering tensor at specific angles and store in 
        an array of shpae (4,4,#angles, ...) 
        """
        # calculate the angle sampling 
        st = self.scattering_tensor(angles)
        # populate a (4,4,...) array with the returned values
        out_shape = (4,4)+st.shape[1:]
        out = sp.zeros(out_shape)
        out[0,0,...] = st[0]
        out[1,1,...] = st[1]
        out[2,2,...] = st[2]
        out[3,3,...] = st[3]
        out[0,1,...] = st[4]
        out[1,0,...] = st[4]
        out[2,3,...] = st[5]
        out[3,2,...] = -st[5]
        return out

    def reshape_lmax(self, lmax):
        " Change the length of the second dimension to new lmax value "
        oldlmax = self.shape[1]
        if lmax <= oldlmax:
            return self[:,:lmax,...]
        else:
            new_shape = sp.array(self.shape)
            new_shape[1] = lmax
            gsf_data = sp.zeros(new_shape)
            gsf_data[:, 0:oldlmax,...] = self
            return SingleScatteringArray(gsf_data)
        
    def sparsify(self, lmax=None):
        " convert a gsf array to a gsf array to a SingleScattering "
        # get the data to make a SingleScattering
        # -----------------------------------------
        lmax = lmax if lmax else self.lmax()
        oldlmax = self.lmax()
        ext = self.ext()
        gsf = self.gsf()

        # make sparse array
        # -------------------
        gsf_lil = sparse.lil_matrix((6, self.GSFMAX))
        l = min(lmax, oldlmax)
        gsf_lil[:,:l] = gsf[:,:l]
        gsf_csr = gsf_lil.tocsr()
        return SingleScattering(lmax, ext, gsf_csr)

    def norm(self, angles=None, ord=2):
        """ 
        Compute one number to represent the norm.
        the ord is passed to 
          scipy.linalg.norm(array, ord=ord)
        To compute a particular norm, ord = 2 for largest singular
        value norm by default.
        """        
        if angles is None: 
            Na = min(20+self.lmax(), 60)
            angles = sp.linspace(0,sp.pi, Na)
        else: 
            assert angles.ndim == 1
            Na = angles.size

        # Calculate angular sampling 
        scat = self.scatmats(angles)

        # Calculate the norm        
        tran_vec = range(scat.ndim)[2:]+[0,1]
        out_shape = sp.array(scat.shape)[range(scat.ndim)[2:]]
        norms = sp.array([la.norm(a, ord=ord)
                         for a in scat.transpose(tran_vec).reshape(-1,4,4)])
        norms = norms.reshape(out_shape).max(0)
        return norms
        



# Packing extinction into the same array as the gsf coefficients
def _pack(ext, gsf, copy=True):
    """
    Store the gsf data in an appropiatly sized array.
    """
    gsfdata = gsf.copy() if copy else gsf # copy array
    gsfdata[5, 0] = ext
    return gsfdata

def _unpack(gsfdata, copy=True):
    """
    unpack an array 
    """
    gsf = gsfdata.copy() if copy else gsfdata # copy the array 
    ext = gsf[5, 0]
    gsf[5, 0] = 0.0
    return ext, gsf

def ss_array(ext, gsf):
    """
    Make a SingleScatteringArray object out of (ext, gsf)
    """
    return SingleScatteringArray(_pack(ext, gsf, copy=False))
    



# script
if __name__=="__main__":

    # With doctests
    import doctest
    doctest.testmod()
            
    # sanity testg 
