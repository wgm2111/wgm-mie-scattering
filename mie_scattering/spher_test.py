
from scipy import linspace, array

from spher import spher_, spher_sum

(cext, 
 csca, 
 a1, a2, a3, a4, 
 b1, b2, 
 lmax)  = spher_(.41, .015, 1.53, .008, 1.0e-13)

print 'cext = ', cext
print 'csca = ', csca
print 'alp1 = ', csca*a1[0:4], csca*a1[-4:] 
print 'alp2 = ', csca*a2[0:4], csca*a2[-4:] 
print 'alp3 = ', csca*a3[0:4], csca*a3[-4:] 
print 'alp4 = ', csca*a4[0:4], csca*a4[-4:] 
print 'bet1 = ', csca*b1[0:4], csca*b1[-4:] 
print 'bet2 = ', csca*b2[0:4], csca*b2[-4:] 




rad = array([.015*(i+1) for i in range(10)])
wrad = 0.0 * rad
wrad[0] += 1.0


(cext, csca, 
 a1, a2, a3, a4, 
 b1, b2, 
 lmax) = spher_sum(.41, rad, wrad, 1.53, .008, 1.0e-13)


print 'cext = ', cext
print 'csca = ', csca
print 'alp1 = ', a1[0:4], a1[-4:] 
print 'alp2 = ', a2[0:4], a2[-4:] 
print 'alp3 = ', a3[0:4], a3[-4:] 
print 'alp4 = ', a4[0:4], a4[-4:] 
print 'bet1 = ', b1[0:4], b1[-4:] 
print 'bet2 = ', b2[0:4], b2[-4:] 
