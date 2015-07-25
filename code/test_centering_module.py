# coding=utf-8
#!/usr/bin/python
"""
This code tests the module 'centering', calling the function centering. 
The fits files are read here, and then the relevant functions called.

Andrés A. Plazas, JPL, July 2015
"""
from centering import *
#import centering
import sys
import subprocess as S
import matplotlib
matplotlib.use('Pdf')
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf  import PdfPages
import matplotlib.font_manager as fm
from scipy.stats import nanmean
from scipy.stats import nanstd
import pyfits

if (len(sys.argv) != 4 ):
    print >>sys.stderr, "Usage: centering.py <best_center_cat.fits> <second_best_center_cat.fits> <ofname.pdf> "
    print >>sys.stderr, "The catalogs should have the fields: 'rsum', 'wsum', 'dsum', 'mem_match_id', 'lambda_chisq', 'lambda_chisq_e', each one with nbins entries. There should be an entry 'p_cen' with the probabilities of the centers, too. "
    sys.exit(1)

data1  = pyfits.open (sys.argv[1])[1].data
data2 =  pyfits.open (sys.argv[2])[1].data
ofname = sys.argv[3]

pp = PdfPages(ofname)

#### Read data 
## Best center
id1=data1['mem_match_id']
lambda1=data1['lambda_chisq']
lambda1_err=data1['lambda_chisq_e']
rsum1=data1['rsum']
dsum1=data1['dsum']
wsum1=data1['wsum']
rsum1/=wsum1
mask_wsum1=[]
for x in wsum1:
    if (np.all(x) == 0.):
        mask_wsum1.append(False)
    else:
        mask_wsum1.append(True)
mask_wsum1=np.array(mask_wsum1)
p1=data1['p_cen'][:,0]



## Second best center
id2=data2['mem_match_id']
lambda2=data2['lambda_chisq']
lambda2_err=data2['lambda_chisq_e']
#p2= 1 - p1
rsum2=data2['rsum']
dsum2=data2['dsum']
wsum2=data2['wsum']
rsum2/=wsum2
mask_wsum2=[]
for x in wsum2:
    if (np.all(x) == 0.):
        mask_wsum2.append(False)
    else:
        mask_wsum2.append(True)
mask_wsum2=np.array(mask_wsum2)


## MATCH clusters by ID. Need to do two masks. The length of teh catalogs might be different. 
print "Matching the objects from the two catalogs. "

mask_id1= np.in1d (id1, id2) 

# apply mask
id1=id1[mask_id1]
lambda1=lambda1[mask_id1]
lambda1_err=lambda1_err[mask_id1]
p1=p1[mask_id1]
rsum1=rsum1[mask_id1]
dsum1=dsum1[mask_id1]
wsum1=wsum1[mask_id1]
osum1=osum1[mask_id1]

mask_id2= np.in1d (id2, id1)
id2=id2[mask_id2]
lambda2=lambda2[mask_id2]
lambda2_err=lambda2_err[mask_id2]
rsum2=rsum2[mask_id2]
dsum2=dsum2[mask_id2]
wsum2=wsum2[mask_id2]
osum2=osum2[mask_id2]

assert ( np.all(id1 == id2) ) # If not True, then they are not matched! 
p2 = 1-p1

r, X_T, X_F, cov , rcoeff = true_profile (rsum1, dsum1, dsum2, wsum1, wsum2, lambda1, lambda1_err, lambda2, lambda2_err, p1, p2)
#For now, cov is just zeroes.
errors = np.diag (cov)

#plot parameters
loc_label='lower left'
prop = fm.FontProperties(size=9)
marker_size=8

## Plot correlation coefficient first
fig = plt.figure()
ax = fig.add_subplot(111)
plt.errorbar (r, rcoeff, yerr=None, fmt='g.-', label='correlation', markersize=marker_size)
plt.xlabel ('R (h$^{-1}$ Mpc)')
plt.ylabel(r'Correlation Coefficient')
plt.xscale('log')
ax.legend(loc=loc_label , fancybox=True, ncol=1, numpoints=1, prop = prop)
pp.savefig()


## Now plot the Delta Sigma profiles about the true and the false centers
fig = plt.figure()
ax = fig.add_subplot(111)
plt.errorbar (r, X_T, yerr=None, fmt='b.-', label='Profile about T center', markersize=marker_size)
plt.errorbar (r, X_F, yerr=None, fmt='r.-', label='Profile about F center', markersize=marker_size)
plt.xlabel ('R (h$^{-1}$ Mpc)')
plt.ylabel(r'$\Delta \Sigma$ (M$_{\odot}$ h/pc$^2$)')
plt.xscale('log')
plt.yscale('log')
ax.legend(loc=loc_label , fancybox=True, ncol=1, numpoints=1, prop = prop)
pp.savefig()
pp.close()



