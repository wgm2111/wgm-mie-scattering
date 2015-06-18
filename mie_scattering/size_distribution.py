
#
# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2010
# liscence: BSD style
#
# ===============================================================================

"""
Module for constructing particle size distributions (actually probability
density functions of radius).  Support for lognormal and gama distributions 
is provided, and there is also support for piecewise smooth polynomials 
known as b-splines.  

Each distribution can be constructed with its particular class definition: 
(1) ModeBasis - for bsplines
(2) LognormAerosol - for lognomal size distribution 
(3) Gamma Aerosol - for a gamma size distribution 

For (1) : 

    ModeBasis and ModeVec are objects used to make a b-spline basis to 
    represent Cinf compact support functions on a given interval.
 
    # Construct a ModeBasis with a interval, (a,b), and a number of bins, n. 
    mb = ModeBasis(a,b,n)

    # Construct a ModeVec by passing mb an array of coefficients or an rv_aerosol 
    mv = mb(*args)                 


For (2  or 3) : 

    LognormAerosol and GammaAerosol subclass rv_continuous object in the 
    scipy.stats.distributions module.

    # Construct LognormAerosol 
    la = LognormAerosol(reff, veff, rmin=None, rmax=None)

    # Construct GammaAerosol 
    ga = GammaAerosol(reff, veff, rmin=None, rmax=None)


* Notes *
For lognormal and gamma distributions rmin and rmax are optional, and they 
help routines that average mie scattering properties.  

"""

# standard python imports
from __future__ import print_function
import os
import shelve
import copy

# python imports
from copy import deepcopy

# scipy/numpy imports
from scipy import linspace, array, ones, zeros, eye, arange, dot, sqrt, log
import scipy as sp
from scipy.interpolate import splrep, splev, splint
from scipy.integrate import quad, romberg
import scipy.sparse as sparse
from numpy import ndarray, unique
import numpy.testing as ntest

from scipy.stats.distributions import gamma_gen, lognorm_gen, rv_continuous 


# my imports 


# Routines for evaluating |vec| and grad|vec| given a norm_array
def evaluate_norm(norm_array, vec):
    """
    Evaluate the 'functional' norm of vec corresponding to norm_array
    """
    return sqrt(dot(dot(norm_array, vec), vec))

def evaluate_dnorm(norm_array, vec):
    """
    Evaluate the derivitive of the 'functional' norm of vec with respect
    to scale parameters stored in vec. 
    """
    return 1.0 / evaluate_norm(norm_array, vec) * dot(vec, norm_array)




# function and class definitions 
class ModeBasis(list):
    """ 
    Itterator, manager, and constructor for a b-spline basis for the number
    concentration pdf for an aerosol population's size parameter. 

    ** Supported functionality **

      polynomial basis for C^k compact support 
      on interval (a, b) with nbins basis elements 
      polynomial order k=1, 3, and 5 for uniform nodes
      l2norm and dl2norm calculations (for k=1,3,5)
      h2semi and dh2semi calculations (for k=5 only)
      persistant l2norm and h2semi matrix calculation
     
    ** ToDo **
      h2semi norm calculations for lower order polynomials k=1,3. 

    """
    
    # Open a shelve and clear if it is old
    # ----------------------------------------
    fname_root = ".size_parameter"
    fname = (os.path.abspath(__file__).split(__file__.split('/')[-1])[0] 
             + fname_root)

    database = shelve.open(fname) # open database file
    path = os.path.abspath(__file__)

    # Wipe the database clean if the module is older than the db files
    try:
        module_mtime = os.path.getmtime(path)
        db_files = [f for f in os.listdir('.') if (fname_root in f)]
        if len(db_files) > 0:
            database_ctime = max(
                *[os.path.getctime(f) for f in db_files])
            if module_mtime > database_ctime: 
                database.clear()
    except OSError as err:
        print("OSError when getting size_parameter database ctime:")
        print(err)
    except Exception as ex:
        raise ex

    # Template for printing
    _print_msg = "ModeBasis (list of basis funcitons):"\
                 "\t {0} B-spline's of order {1} on interval ({2},{3})."


    def __init__(self, a, b, nbins, k=5):
        " construct a b-spline basis manager for nodes x "

        # Check if case is precomputed
        # -------------------------------
        # database = shelve.open(ModeBasis.fname) # open database file
        _a = repr(a)[:-2]
        _b = repr(b)[:-2]
        self.input_key = _a + _b + repr(nbins) + repr(k)

        # Check incoming arguments for sanity
        # ---------------------------------------
        assert (k==1 or k==3 or k==5    # Only supported k is k=3
                ), "k={0} is not supported, use odd values".format(k)
        assert 0.0 < a < b          # need a<b both positive 

        # determine knots
        # --------------------
        nbins = int(nbins)
        pad = k+1                            # pad points to close ends        
        x_open = linspace(a, b, nbins+2)[1:-1] # open interval
        open_slice = slice(pad, pad+nbins, 1)
        x_left = linspace(a, x_open[0], pad+1)[:-1]
        x_right = linspace(x_open[-1], b, pad+1)[1:]

        # define the padded x (for splrep)
        # ------------------------------------
        x_pad = zeros(2*pad + nbins)    
        x_pad[0:pad] = x_left
        x_pad[pad:-pad] = x_open
        x_pad[-pad:] = x_right

        # big and element-wise supports
        # ------------------------------------
        rmin = a
        rmax = b
        steps = 1 + (k-1)/2
        rmins = [x_pad[pad-steps+i] for i in range(nbins)]
        rmaxs = [x_pad[pad+steps+i] for i in range(nbins)]

        # define t, c, k and store the normalized values in c
        # -----------------------------------------------------
        t, c, k = splrep(x_pad, zeros(x_pad.size), k=k)
        def tck_maker(i):
            " return the tck rep for the ith basis function "
            out = 0 * c
            out_open = out[open_slice]
            out_open[i] = 1.0
            return (t, out, k)
        cnorm = 0.0 * c
        cnorm[open_slice] = array([1.0 /
                                   splint(rmins[i], rmaxs[i], tck_maker(i))
                                   for i in range(nbins)])
        tck_norm = (t, cnorm, k)

        # make self a list of ModeVec objects
        # --------------------------------------
        eyedentity = eye(nbins)
        argslist = zip([v for v in eyedentity], rmins, rmaxs)
        keys = ['tck_norm', 'open_slice', 'mode_basis']
        vals = [tck_norm, open_slice, self]
        kwds = dict(zip(keys, vals))
        list.__init__(self, (ModeVec(*args, **kwds) for args in argslist))

        # hang attributes
        # -----------------
        self.nbins = nbins
        self.k = k
        self.rmin, self.rmax = rmin, rmax
        self.x_open, self.x_pad = x_open, x_pad
        self.cnorm = cnorm
        self.tck_norm = tck_norm
        self.open_slice = open_slice
        self.argslist = argslist
        self.kwds = kwds
        (self.l2norm_array, 
         self.h1semi_array) = self.get_norm_array()

        (self.l2norm_array_volume, 
         self.h1semi_array_volume) = self.get_norm_array_volume()

    def __call__(self, argin):
        """
        construct a ModeVec object with either a vector of number
        concentrations or a rv_aerosol object 
        """
        # assemble arguemnets for ModeVec from a ndarray object
        #----------------------------------------------------------
        if hasattr(argin, 'size'):
            try: 
                assert argin.size == self.nbins
            except: 
                temp = 'expected input array of size {0}, got {1}'
                msg = temp.format(self.nbins, argin.size)
                raise ValueError(msg)
            nconc = argin
            args = (nconc, self.rmin, self.rmax)
            kwds = self.kwds

        # assemble arguments for a analytic aerosol object
        #----------------------------------------------------------
        elif hasattr(argin, 'pdf'):
            # shortname variables
            x_pad = self.x_pad
            # expand pdf in bsplines
            t, c, k = splrep(x_pad, argin.pdf(x_pad), k=self.k)
            nconc = c[self.open_slice] / self.cnorm[self.open_slice]
            args = (nconc, self.rmin, self.rmax)
            kwds = self.kwds

        else: 
            raise ValueError("ModeBasis.__call__ input error")

        # construct and return the ModeVec
        return ModeVec(*args, **kwds)
 
    def __repr__(self):
        " Define for reasonable looking print output."
        out = self._print_msg.format(self.nbins, self.k, self.rmin, self.rmax)
        return out
# msg = "ModeBasis(list): "\
#               "\n\t min_radius = {0}, \tmax_radius={1}" \
#               "\n\t Nbins = {2}, \t\tpolynomial order={3}."
#         out = msg.format(self.rmin, self.rmax, self.nbins, self.k)
#         return out
        
    

    def has_uniform_nodes(self):
        " Return True if the nodes are uniform False if they are not "
        k = self.k
        x = self.x_open[k+1:-(k+1)-1]
        dx = (x[1:] - x[:-1]).round(9)
        # print('dx = {0}'.format(dx))
        if unique(dx).size == 1:
            return True
        else:
            return False

    def norm(self, vec, kind='l2norm'):
        " Evaluate a norm of specified kind."
        # argument checking
        assert vec.size == self[0].size
        assert kind in ['l2norm', 'h1semi',
                        'l2norm_volume', 'h1semi_volume']
        # compute the norm/seminorm
        if kind == 'l2norm':
            return evaluate_norm(self.l2norm_array, vec)
        if kind == 'h1semi':
            return evaluate_norm(self.h1semi_array, vec)
        # if kind == 'h2semi':
        #     return evaluate_norm(self.h2semi_array, vec)
        if kind == 'l2norm_volume':
            return evaluate_norm(self.l2norm_array_volume, vec)
        if kind == 'h1semi_volume':
            return evaluate_norm(self.h1semi_array_volume, vec)

    def l2norm(self, vec): 
        return self.norm(vec, 'l2norm')

    def h1semi(self, vec): 
        return self.norm(vec, 'h1semi')

    def l2norm_volume(self, vec): 
        return self.norm(vec, 'l2norm_volume')
    def h1semi_volume(self, vec): 
        return self.norm(vec, 'h1semi_volume')
        

    def dnorm(self, vec, kind='l2norm'):
        """
        Evaluate a the derivitive of the specified kind.  This is taken with 
        respect to the mode vector parameters
        """
        # argument checking
        assert vec.size == self[0].size
        assert kind in ['l2norm', 'h1semi', 
                        'l2norm_volume', 'h1semi_volume']
        # compute the norm/seminorm
        if kind == 'l2norm':
            return evaluate_dnorm(self.l2norm_array, vec)
        if kind == 'h1semi':
            return evaluate_dnorm(self.h1semi_array, vec)
        if kind == 'l2norm_volume':
            return evaluate_dnorm(self.l2norm_array_volume, vec)
        if kind == 'h1semi_volume':
            return evaluate_dnorm(self.h1semi_array_volume, vec)


    def dl2norm(self, vec): 
        return self.dnorm(vec, 'l2norm')
    def dh1semi(self, vec): 
        return self.dnorm(vec, 'h1semi')
    def dl2norm_volume(self, vec): 
        return self.dnorm(vec, 'l2norm_volume')
    def dh1semi_volume(self, vec): 
        return self.dnorm(vec, 'h1semi_volume')


    def get_norm_array(self):
        """
        Calculate the matricies needed to compute the l2norm and h2seminorm
        of superpositions of bin-functions (b-splines).  The routine is faster
        for uniformly spaced bin functions. 
        
        The norm matrix is the matrix M such that
        
        ||x_function||  =  sqrt(dot(x, dot(M, x))
        D_x ||x_function||  = dot(x, M)) / ||x_function||
        
        ** input **
        mode_basis_instance - an instance of the class ModeBasis 
        
        ** output **
        (l2norm_array, h1semi_array) - two norm arrays
        
        """

        # try looking up in the shelve
        # -----------------------------
        ikey = self.input_key
        database = self.database
        if ikey in database:
            l2norm, h1semi = database[ikey]
            return l2norm, h1semi
        
        
        # initialize the norm_arrays
        # ----------------------------
        _size = self[0].size
        l2norm = zeros((_size, _size)) 
        h1semi = zeros((_size, _size)) # initialized

        # Slower version for non-uniform nodes
        # ---------------------------------------
        
        # Loop over diagonal offset to compute values
        for k in range(self.k+1):
            # Integrate to get current (off)diagonal of the norm matricies
            out = array(
                [array([quad(lambda r: f.__call__(r)*g.__call__(r), 
                             max(f.rmin, g.rmin), min(f.rmax, g.rmax))[0],
                        quad(lambda r: f.__call__(r,der=1)*g.__call__(r,der=1),
                             max(f.rmin, g.rmin), min(f.rmax, g.rmax))[0]])
                 for f, g in zip(self, self[k:])])
            I = arange(0, _size-k)
            J = arange(k, _size)
            l2norm[I,J] = out[:,0]
            h1semi[I,J] = out[:,1]
            if k!=0: 
                l2norm[J,I] = out[:,0]
                h1semi[J,I] = out[:,1]

        # store the new calculation and return 
        # -------------------------------------
        database[ikey] = (l2norm, h1semi)
        return l2norm, h1semi


    # --------------------------------------------------------
    # Get arrays for computing the volume norm of the mode
    # --------------------------------------------------------
    def get_norm_array_volume(self):
        """
        Calculate the matricies needed to compute the l2norm and h2seminorm
        of superpositions of bin-functions (b-splines).  The routine is faster
        for uniformly spaced bin functions. 
        
        The norm matrix is the matrix M such that
        
        ||x_function||  =  sqrt(dot(x, dot(M, x))
        D_x ||x_function||  = dot(x, M)) / ||x_function||
        
        ** input **
        mode_basis_instance - an instance of the class ModeBasis 
        
        ** output **
        (l2norm_array, h1semi_array) - two norm arrays
        
        """
        # try looking up in the shelve
        # -----------------------------
        ikey = 'vol' + self.input_key
        database = self.database
        if ikey in database:
            l2norm = database[ikey]
            return l2norm

        # initialize the norm_arrays
        # ----------------------------
        _size = self[0].size
        l2norm = zeros((_size, _size)) 
        h1semi = zeros((_size, _size)) 
        
        # Slower for non-uniform nodes
        # ----------------------------
        def get_h1_integrand(f, g):
            def integrand(r):
                fvals = f.__call__(r)
                dfvals = f.__call__(r, der=1)
                gvals = g.__call__(r)
                dgvals = g.__call__(r, der=1)
                vals = (4*sp.pi / 3.0)**2 * ((3 * r**2 * fvals + r**3 * dfvals)*
                                             (3 * r**2 * gvals + r**3 * dgvals))
                return vals
            return integrand
        
        # Loop over diagonal offset to compute values
        for k in range(self.k+1):
            # Integrate to get current (off)diagonal of the norm matricies
            
            # for l2norm
            out = array(
                [(4*sp.pi/3.0)**2 * quad(lambda r: f.__call__(r)*g.__call__(r)*r**6, 
                                         max(f.rmin, g.rmin), min(f.rmax, g.rmax))[0]
                 for f, g in zip(self, self[k:])])
            I = arange(0, _size-k)
            J = arange(k, _size)
            l2norm[I,J] = out[:]
            if k!=0: 
                l2norm[J,I] = out[:]
                
            # for h1semi
            out = array(
                [quad(get_h1_integrand(f,g), max(f.rmin, g.rmin), min(f.rmax, g.rmax))[0]
                 for f, g in zip(self, self[k:])])
            # I = arange(0, _size-k)
            # J = arange(k, _size)
            h1semi[I,J] = out[:]
            if k!=0: 
                h1semi[J,I] = out[:]
                
        # store the new calculation and return 
        # -------------------------------------
        database[ikey] = l2norm, h1semi
        return l2norm, h1semi




# Subclass ndarray to store a basis representation of a aerosol mode 
class ModeVec(ndarray):
    " vector storing the nconc values for a b-spline representation "

    # added because multiprocessing was looking for these
    rmin = None
    rmax = None

    def __new__(cls, nconc, rmin, rmax, mrr=None, mri=None,
                tck_norm=None, open_slice=None, mode_basis=None):
        " backend construction of a ModeVec object "
        # cast as ModeVec class array
        out = nconc.view(cls)
        # hang data 
        out.tck_norm = deepcopy(tck_norm)
        out.k = tck_norm[2]
        # out.mrr, out.mri = mrr, mri
        out.rmin, out.rmax = rmin, rmax
        out.open_slice = open_slice
        out.mode_basis = mode_basis
        return out

    def __array_finalize__(self, obj):
        " special finalize to hang attributes "
        if obj is None: return
        self.tck_norm = getattr(obj, 'tck_norm', None)
        # self.mrr = getattr(obj, 'mrr', None)
        # self.mri = getattr(obj, 'mri', None)
        self.rmin = getattr(obj, 'rmin', None)
        self.rmax = getattr(obj, 'rmax', None)
        self.open_slice = getattr(obj, 'open_slice', None)
        self.mode_basis = getattr(obj, 'mode_basis', None)

    def __array_wrap__(self, array_out, contest=None):
        " wrap data on array again after ufunc "
        array_out = array_out.view(self.__class__)
        ModeVec.__array_finalize__(array_out, self)
        return array_out

    def __call__(self, x, der=0):
        " evaluate spline for given nconc values "
        t, c, k = deepcopy(self.tck_norm)
        c[self.open_slice] = c[self.open_slice] * self
        return splev(x, (t, c, k), der=der)    

    def ndf(self, r):
        return self.__call__(r)

    def vdf(self, r):
        " Evaluate the spline and scale by volume "
        volume = 4/3. * sp.pi * r**3
        return volume * self.__call__(r, der=0)

    # def pdf(self, r):
    #     """ 
    #     evaluate the spline at radii r but normalize to make into a pdf 
    #     """
    #     # since splines are already normalized this job is made simple
    #     return self(r)/(self.sum())

    # def dVdlnr(self, r):
    #     "Evaluate for comparison with Oleg "
    #     return self.vdf(r) / r

    def norm(self, kind='l2norm'):
        " Evaluate the norm of self"
        return self.mode_basis.norm(self, kind=kind)
    def dnorm(self, kind='l2norm'):
        " Evaluate the norm of self "
        return self.mode_basis.dnorm(self, kind=kind)


class mode_basis(ModeBasis):
    """
    Mode basis with support for computing functional norms/seminorms
    """
    pass


# ===========================================================================
# Analytically defined aerosol distributions 
# ===========================================================================
rmaxfactor = 20                 # rmax = reff + rmaxfactor * sqrt(veff)


# define rv_aerosol constructors 
# ===========================================================================

class AnalyticAerosol(object):
    """
    Class for handling the size parameter pdf poly-disperse aerosols.  

    ** members ** 
    reff, veff
    rmin, rmax

    ** methods **
    see scipy.stats.rv_continuous
    

    ** examples **
    >>> reff, veff = .5, .2
    >>> ga = GammaAerosol(reff, veff, mrr, mri)
    >>> la = LognormAerosol(reff, veff, mrr, mri)

    """
    nconc = 1.

    def pdf(self, x):
        " rescale the pdf of rv and return "
        I = (x >= self.rmin) * (x < self.rmax)
        return 1.0 / self.integral * sp.where(I, self.rv.pdf(x), 0.0)

    def ndf(self, x):
        " Number concentration "
        return self.nconc * self.pdf(x)

    def vdf(self, x):
        " Volume concentration "
        return 4.0 / 3.0 * sp.pi * x**3 * self.ndf(x)

    def dVdlnr(self, x):
        " Size distribution for comparing with Oleg Dubovic's dV/dln(r)"
        return self.vdf(x)/x

class GammaAerosol(AnalyticAerosol):
    """
    define a gamma random variable with effective parameters.
    
    """
    __doc__ = __doc__ + AnalyticAerosol.__doc__

    def __init__(self, reff, veff, rmin=None, rmax=None):
        """
        construct a gamma aerosol mode.
        """
        # initialize the stats.rv_continuous with gamma_gen
        # ---------------------------------------------------
        rv = gamma_gen(name="gamma_aerosol", extradoc=self.__doc__)  
        args = ((1.0-2.0*veff) / veff, ) 
        kwargs = dict(scale=reff*veff)#, 
                      # name= 'gamma_aerosol') # prepare to freeze rv
        rv = rv.freeze(*args, **kwargs)      # freeze
        
        # set bounds and get the integral inbounds 
        # ------------------------------------------
        if rmin is None: rmin = max(reff - rmaxfactor*sqrt(veff), 0.0)
        if rmax is None: rmax = reff + rmaxfactor * sqrt(veff)    
        integral_inbounds = quad(rv.pdf, rmin, rmax)[0]

        # make the namespace
        # ---------------------
        self.rv = rv
        self.name = "gamma_aerosol"
        self.reff, self.veff = reff, veff
        self.rmin, self.rmax = rmin, rmax
        self.integral = integral_inbounds
 
class LognormAerosol(AnalyticAerosol):
    """
    Define a lognorm random variable with effective parameters.
    
    """
    __doc__ = __doc__ + AnalyticAerosol.__doc__
    
    def __init__(self, reff, veff, rmin=None, rmax=None):
        """
        Construct a lognormal random variable
        """
        # initialize the stats.rv_continuous with lognorm_gen
        # ------------------------------------------------------
        rv = lognorm_gen(name="lognorm_aerosol", extradoc=self.__doc__) 
        args = (sqrt(log(veff + 1.)), )
        kwargs = dict(scale=reff / (veff+1.)**(5./2))#,
                      # name= 'lognorm_aerosol') # prepare to freeze
        rv = rv.freeze(*args, **kwargs)        # freeze
        
        # set bounds and get the integral inbounds 
        # ------------------------------------------
        if rmin is None: rmin = max(reff - rmaxfactor*sqrt(veff), 0.0)
        if rmax is None: rmax = reff + rmaxfactor * sqrt(veff)
        integral_inbounds = quad(rv.pdf, rmin, rmax)[0]

        # make the namespace
        # ---------------------
        self.rv = rv
        self.name = "lognorm_aerosol"
        self.reff, self.veff = reff, veff
        self.rmin, self.rmax = rmin, rmax
        self.integral = integral_inbounds


class MultiLognormAerosol(object):
    """
    Define a lognorm random variable with effective parameters.
    
    """
    __doc__ = __doc__ + AnalyticAerosol.__doc__
    
    def __init__(self, reffs, veffs, rmin, rmax, scale=1.0):
        """
        Construct a multi-mode lognormal random variable
        """
        
        # Make a list of Lognormal Aerosols
        # ------------------------------------------------------
        self.nmodes = len(reffs)
        self.scale = scale
        self.reffs = reffs
        self.veffs = veffs
        self.rmin = rmin
        self.rmax = rmax
        self.rvlist = [LognormAerosol(reff, veff, rmin=rmin, rmax=rmax)
                       for reff, veff in zip(reffs, veffs)]

    def pdf(self, r):
        " Caculate the pdf evaluated at positions r "
        rvals = sp.array([rv.pdf(r) for rv in self.rvlist])
        return rvals.mean(axis=0)

    def __call__(self, r):
        return self.pdf(r) * self.scale
    
    def ndf(self, r):
        return self.__call__(r)

    def vdf(self, r):
        volume = 4/3. * sp.pi * r**3
        return volume * self.ndf(r)




# ====
# Example script
# ===============================================================================

if __name__=="__main__":


    import matplotlib.pyplot as plt

    # make a couple of mode basis objects
    a = .01
    b = 1.25
    n = 23
    mb = ModeBasis(a,b,n)


    # make a plot showing the bin functions
    plt.close('all')
    f0 = plt.figure(0, (21,8))
    subs = [f0.add_subplot(1,3,i) for i in xrange(1,4)]

    # plot the log mode basis 
    for mv in mb:
        r = sp.linspace(mv.rmin*(1+.5e-2), 
                        mv.rmax*(1-.5e-2))
        y = log(r)
        subs[0].plot(y, mv(r)) 
        subs[1].plot(r, mv.vdf(r))
        subs[2].plot(r, mv.ndf(r))#semilogy(r, mv.ndf(r))
    subs[0].set_title('LogModeBasis volume density')
    subs[0].set_ylabel('vdf')
    subs[0].set_xlabel('y = log(r)')
    subs[1].set_title('LogModeBasis volume density')
    subs[1].set_ylabel('vdf')
    subs[1].set_xlabel('r')
    subs[2].set_title('LogModeBasis number density')
    subs[2].set_ylabel('ndf')
    subs[2].set_xlabel('r')

    # plot the standard mode basis

    f0.show()


    # make a lognormal aerosol to test out interpolation 
    reff = .24
    veff = .04
    rv = LognormAerosol(reff, veff, a, b)
    rv_mv = mb(rv)

    r = sp.linspace(a, b, 250)
    f1 = plt.figure(1, (10,10))
    s = f1.add_subplot(1,1,1)
    s.plot(r, rv.ndf(r), label="rv")
    s.plot(r, rv_mv.ndf(r), label="ModeVec")

    s.set_title('Interpolation comparison')
    s.legend()
    f1.show()


    






