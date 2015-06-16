C  The program computes far-field light scattering by polydisperse     
C  homogeneous spherical particles using the Lorenz-Mie theory.    
                                                                     
C  Last modified 08/06/2005                                            
                                                                        
C  The program is based on the book                                    
C 
C  1. Mishchenko, M. I., L. D. Travis, and A. A. Lacis (2002). 
C     Scattering, Absorption, and Emission of Light
C     by Small Particles. Cambridge University Press, Cambridge.
C
C  The book is available in the .pdf format at the following
C  web site:
C
C   http://www.giss.nasa.gov/~crmim/books.html
C                                                                      
C  One can also use the following paper as a brief user guide:               
C                                                                      
C  2. M. I. Mishchenko, J. M. Dlugach, E. G. Yanovitskij,            
C     and N. T. Zakharova, Bidirectional reflectance of flat,        
C     optically thick particulate laters:  an efficient radiative    
C     transfer solution and applications to snow and soil surfaces,  
C     J. Quant. Spectrosc. Radiat. Transfer, Vol. 63, 409-432 (1999).         
C
C  This paper is available upon request from Michael Mishchenko
C  (please send a mesage to crmim@giss.nasa.gov) and is also 
C  available in the .pdf format at the following web site:
C 
C  http://www.giss.nasa.gov/~crmim/publications/index.html 
                                                                       
C  The code has been developed by Michael Mishchenko at the NASA       
C  Goddard Institute for Space Studies, New York.  The development     
C  of the code was supported by the NASA Radiation Sciences Program.    
                                                                       
C  The code can be used without limitations in any not-for-            
C  profit scientific research.  The only request is that in any            
C  publication using the code the source of the code be acknowledged   
C  and proper references be made.                                      
                                                                       
                                                                       
C  Input parameters:                                                   
C
C     NDISTR specifies the type of particle size distribution 
C     as follows:        
C
C     NDISTR = 1 - modified gamma distribution                         
C          [Eq. (5.242) of Ref. 1]                                        
C              AA=alpha                                                
C              BB=r_c                                                  
C              GAM=gamma                                               
C     NDISTR = 2 - log normal distribution                             
C          [Eq. (5.243) of Ref. 1]                                        
C              AA=r_g                                                  
C              BB=[ln(sigma_g)]**2                                      
C     NDISTR = 3 - power law distribution                              
C          [Eq. (5.244) of Ref. 1]                                        
C               AA=r_eff (effective radius)                            
C               BB=v_eff (effective variance)                          
C               Parameters R1 and R2 (see below) are calculated        
C               automatically for given AA and BB                      
C     NDISTR = 4 - gamma distribution                                  
C          [Eq. (5.245) of Ref. 1]                                        
C               AA=a                                                   
C               BB=b                                                   
C     NDISTR = 5 - modified power law distribution                     
C          [Eq. (5.246) of Ref. 1]                                        
C               BB=alpha                                               
C     NDISTR = 6 - bimodal volume log normal distribution              
C              [Eq. (5.247) of Ref. 1]             
C              AA1=r_g1                                                
C              BB1=[ln(sigma_g1)]**2                                   
C              AA2=r_g2                                                
C              BB2=[ln(sigma_g2)]**2                                   
C              GAM=gamma                                               
C
C   - R1 and R2 - minimum and maximum                                  
C         radii in the size distribution for NDISTR=1-4 and 6.         
C         R1 and R2 are calculated automatically                       
C         for the power law distribution with given r_eff and v_eff    
C         but must be specified for other distributions.               
C         For the modified power law distribution (NDISTR=5), the      
C         minimum radius is 0, R2 is the maximum radius,               
C         and R1 is the intermediate radius at which the               
C         n(r)=const dependence is replaced by the power law           
C         dependence.                                                  
C
C   - LAM - wavelength of the incident light in the surrounding medium
C
C   - Important:  r_c, r_g, r_eff, a, LAM, R1, and R2                  
C                 must be in the same units of length (e.g., microns)  
C
C   - MRR and MRI - real and imaginary parts of the relative refractive         
C                   index (MRI must be non-negative)                        
C
C   - N - number of integration subintervals on the interval (R1, R2)   
C         of particle radii                                            
C   - NP - number of integration subintervals on the interval (0, R1)   
C         for the modified power law distribution                      
C   - NK - number of Gaussian division points on each of the           
C          integration subintervals                                    
C
C   - NPNA - number of scattering angles at which the scattering       
C        matrix is computed                                            
C        (see the PARAMETER statement in subroutine MATR).             
C        The corresponding scattering angles are given by              
C        180*(I-1)/(NPNA-1) (degrees), where I numbers the angles.     
C        This way of selecting scattering angles can be easily changed 
C        at the beginning of subroutine MATR by properly modifying     
C        the following lines:                                          
C                                                                      
C          N=NPNA                                                       
C          DN=1D0/DFLOAT(N-1)                                          
C          DA=DACOS(-1D0)*DN                                           
C          DB=180D0*DN                                                 
C          TB=-DB                                                      
C          TAA=-DA                                                     
C          DO 500 I1=1,N                                               
C             TAA=TAA+DA                                               
C             TB=TB+DB                                                 
C                                                                      
C        and leaving the rest of the subroutine intact.                
C        This flexibility is provided by the fact that                 
C        after the expansion coefficients AL1,...,BET2 (see below)     
C        are computed by subroutine SPHER, the scattering matrix       
C        can be computed for any set of scattering angles without      
C        repeating Lorenz-Mie calculations.                                   
C
C   - DDELT - desired numerical accuracy of computing the scattering   
C        matrix elements                                               
                                                                       
                                                                       
C  Output information:                                                 
C                                                                      
C   - REFF and VEFF - effective radius and effective variance          
C          of the size distribution [Eqs. (5.248)-(5.250) of Ref. 2]   
C   - CEXT and CSCA - average extinction and scattering cross sections      
c          per particle                                         
C   - <COS> - asymmetry parameter                 
C   - ALBEDO - single-scattering albedo                                
C   - <G> - average projected area per particle                        
C   - <V> - average volume per particle                                
C   - Rvw - volume-weighted average radius                             
C   - <R> - average radius                                                
C   - F11 = a_1, F33 = a_3, F12 = b_1, and F34 = b_2 - elements of the 
C          normalized scattering matrix given by Eq. (4.65) of Ref. 1.       
C   - ALPHA1, ..., BETA2 - coefficients appearing in the expansions
C          of the elements of the normalized scattering matrix in
C          generalized spherical functions [Eqs. (4.75)-(4.80) of Ref. 1].              
C
C   The first line of file 10 provides the single-scattering
C   albedo and the number of numerically significant 
C   expansion coeffieients ALPHA1, ..., BETA2 followed
C   by the values of the expansion coefficients.
C                                                                      
C   Optical efficiency factors QEXT and QSCA are computed as           
C   QEXT=CEXT/<G> and QSCA=CSCA/<G>. The absorption cross section     
C   is given by CABS=CEXT-CSCA. The absorption efficiency factor      
C   is equal to QABS=QEXT-QSCA.                                        
                                                                       
C   To calculate scattering and absorption by a monodisperse 
C   particle with radius r = AA, use the following options:          
C                                                                      
C        BB=1D-1                                                       
C        NDISTR=4                                                      
C        NK=1                                                          
C        N=1                                                           
C        R1=AA*0.9999999 D0                                            
C        R2=AA*1.0000001 D0                                            
                                                                       
C   Although the mathematical value of the parameter R2 for            
C   modified gamma, log normal, and gamma distributions is infinity,   
C   in practice a finite value must be chosen.  Two interpretations    
C   of a truncated size distributions can be given.  First, R2 can     
C   be increased until the scattering cross sections and the           
C   normalized scattering matrix elements converge within some numerical          
C   accuracy.  In this case the converged truncated                    
C   size distribution is numerically equivalent to the distribution    
C   with R2=infinity.  Second, a truncated distribution can be         
C   considered a specific distribution with optical cross sections     
C   and scattering matrix elements different from those for the        
C   distribution with R2=infinity.  In this latter case the effective  
C   radius and effective variance of the truncated distribution        
C   will also be different from those computed analytically for the    
C   respective distribution with R2=infinity.              
C   Similar considerations apply to parameter R1 whose theoretical     
C   value for modified gamma, log normal, and gamma distributions      
C   is zero, but in practice can be any number smaller than R2.        
                                                                       
C   Important:  for given size distribution parameters NDISTR, AA,     
C   BB, GAM, R1, AND R2, the                                           
C   parameters N, NP, and/or NK should be increased until convergent   
C   results are obtained for the extinction and scattering cross sections       
C   and, especially, the elements of the scattering matrix.                 
 
C   I would highly appreciate information on any problems encountered
C   with this code.  Please send your message to Michael Mishchenko
C   at crmim@giss.nasa.gov.

C   While this computer program has been tested for a variety of cases,
C   it is not inconceivable that it contains undetected errors. Also,
C   input parameters can be used which are outside the envelope of
C   values for which results are computed accurately. For this reason,
C   the author and his organization disclaim all liability for
C   any damages that may result from the use of the program. 
