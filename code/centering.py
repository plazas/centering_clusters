#!/usr/bin/python
import numpy as np
from scipy.special import erf
from scipy.stats import nanmean
from scipy.stats import nanstd
import math
from numpy import nansum


"""
This module contains functions to implement the equations in the writeup 
'How to use the centering probabilities of galaxy clusters for cosmological investigations' by 
E. Rozo and M. Becker, which use the cenetring probabilities given by RedMapper to calculate the 
profile of an observable X (e.g, \Delta \Sigma) around the true center of a cluster from 
the most and second most likely measured centers.


input files: 1. lens catalog of \Delta \Sigma around most probale center
             2. lens catalog of \Delta \Sigma around second most probable center
     
     
# Assumes that p1 + p2 = 1
     
A. Plazas (JPL, July-2015)
"""


def compute_binning_probability ( lambda_best, lambda_err_best, lambda_2best, lambda_err_2best, p1):
    """
    This function computes the binning probability (Eqs. 18 and 19?)
    Computes the binning probability q_a^\alpha for each cluster alpha. One  number per cluster.
    """ 
    # Inputs are single numbers, not arrays. 
    # Computes the binning probability q_a^\alpha for each cluster alpha. One  number per cluster.
    ## Change of variable for the limits?
    lambda_min, lambda_max = min (lambda_best), max (lambda_best)
    umax= (lambda_max - lambda_best)/math.sqrt(2.0)/lambda_err_best
    umin= (lambda_min - lambda_best)/math.sqrt(2.0)/lambda_err_best
    binprob1=0.5*erf(umax) - 0.5*erf(umin)
    umax= (lambda_max - lambda_2best)/math.sqrt(2.0)/lambda_err_2best     #;; NOTE: catalogs should include lambda_chisq error for alternative centers. ;; using error from original center for now.
    umin= (lambda_min - lambda_2best)/math.sqrt(2.0)/lambda_err_2best
    binprob2=0.5*erf(umax) - 0.5*erf(umin)
    binprob = p1*binprob1 + (1.-p1)*binprob2
    return binprob


def weighted_mean (dsum_vec, wsum_vec):
    """
    m (x; w_x) = \sum (w_x*x) / \sum (w_x)
    """
    assert(len(dsum_vec)==len(wsum_vec))
    num = nansum (dsum_vec, axis=0)
    den = nansum (wsum_vec, axis=0)
    if den.all() == 0.:
        print "Error in function 'wieghted_mean'. den==0. Quitting."
        sys.exit(1)
    return num/den

def weighted_cov ( dsum_x, dsum_y, wsum_x, wsum_y):
    """
    cov (x,y; w_x, w_y) =  m(x*y; w_x*w_y) - m(x,w_x)*m(y, w_y) = \sum(w_x*w_y*x*y)/\sum (w_x*w_y) - m(x, w_x)*m(y,w_y)
    """
    assert(len(dsum_x)==len(dsum_y))
    assert(len(wsum_x)==len(wsum_y))
    return weighted_mean (dsum_x*dsum_y, wsum_x* wsum_y) - weighted_mean(dsum_x, wsum_x)*weighted_mean(dsum_y, wsum_y)

def compute_rcoeff ( dsum_best_vec , dsum_2best_vec, wsum_best_vec, wsum_2best_vec ):
    """
    Function to compute the correlation coefficient based on the covariance matrix.
    r=cov (x,y, w_x, w_y) /  sqrt (cov (x,x, w_x, w_x) *cov(y,y, w_y, w_y))
    """
    dsum_best_vec , dsum_2best_vec = np.array(dsum_best_vec), np.array(dsum_2best_vec)
    wsum_best_vec, wsum_2best_vec = np.array (wsum_best_vec), np.array(wsum_2best_vec)
    assert (len (dsum_best_vec) == len (dsum_2best_vec))
    assert (len (wsum_best_vec) == len (wsum_2best_vec))
    var1 =  weighted_cov (dsum_best_vec, dsum_best_vec, wsum_best_vec, wsum_best_vec)
    rms1 = np.sqrt (var1)
    var2 =  weighted_cov (dsum_2best_vec, dsum_2best_vec, wsum_2best_vec, wsum_2best_vec)
    rms2 = np.sqrt (var2)
    rcoeff = weighted_cov (dsum_best_vec, dsum_2best_vec, wsum_best_vec, wsum_2best_vec)   / (rms1*rms2)
    return rcoeff

def true_profile (rsum1, dsum1, dsum2, wsum1, wsum2, lambda1, lambda1_err, lambda2, lambda2_err, p1):
    """
    rsum1: 
    dsum1,2:
    wsum1,2:
    lambda1,2:
    lambda1,2_err:
    p1:
    """

    ### Derived quantities
    ## binning probability
    binning_prob =  compute_binning_probability (lambda1, lambda1_err, lambda2, lambda2_err, p1)
    delta_sigma1 = dsum1/wsum1   ### ?? Should be only delta_sigma, not dsum 
    delta_sigma2 = dsum2/wsum2  
    #print delta_sigma1.shape
    #print delta_sigma2.shape
    ### Define vectors and matrices
    # reshape probability from nclusters to (nclusters, nradial bins) 
    p=[]
    nradial_bins=len(dsum1[0])
    for x in p1:
        p.append (np.repeat (x, nradial_bins))
    p1=np.array(p)
    X = np.array ([ delta_sigma1, delta_sigma2]) 
    P = np.array ([p1,1-p1])
    #print X.shape
    #print P.shape
    #print "wsum1.shape", wsum1.shape
    rcoeff = compute_rcoeff (dsum1, dsum2, wsum1, wsum2)
    #print "rcoeff: "
    #print rcoeff.shape
    #print rcoeff
    #fig=plt.figure()
    #plt.plot(rcoeff, 'r')
    #pp.savefig()
    c11, c22 =1./wsum1, 1./wsum2
    c12=rcoeff*np.sqrt (c11*c22) 
    #C=np.array ( [[c11,c12], [c12,c22]] )
    S= np.array ( [[c22,-c12],[-c12, c11]]  ) / (c11*c22 - c12**2)   # inverse of C
    #print "S.shape", S.shape
    #print "X.shape", X.shape
    #print "P.shape", P.shape
    Y_T = S[0][0]*P[0]*X[0]  + S[0][1]*P[1]*X[0]  + S[1][0]*P[0]*X[1]  + S[1][1]*P[1]*X[1]
    Y_F = S[0][0]*(1-P[0])*X[0]  + S[0][1]*(1-P[1])*X[0]  + S[1][0]*(1-P[0])*X[1]  + S[1][1]*(1-P[1])*X[1]
    A=S[0][0]*P[0]*P[0]  + S[0][1]*P[0]*P[1]  + S[1][0]*P[0]*P[1]  + S[1][1]*P[1]*P[1]
    B=S[0][0]*(1-P[0])*P[0]  + S[0][1]*(1-P[0])*P[1]  + S[1][0]*P[0]*(1-P[1])  + S[1][1]*(1-P[1])*P[1]
    D=S[0][0]*(1-P[0])*(1-P[0])  + S[0][1]*(1-P[0])*(1-P[1])  + S[1][0]*(1-P[0])*(1-P[1])  + S[1][1]*(1-P[1])*(1-P[1])
    #reshape binning probability
    p=[]
    for x in binning_prob:
        p.append (np.repeat (x, nradial_bins))
    binning_prob=np.array(p)
    Y_T=binning_prob*Y_T
    Y_F=binning_prob*Y_F
    A=binning_prob*A
    B=binning_prob*B
    D=binning_prob*D
    #print "len(Y_T), len(Y_F), len(A), len(B), len(D)", len(Y_T), len(Y_F), len(A), len(B), len(D)
    #print "Y_T.shape, Y_F.shape, A.shape, B.shape, D.shape", Y_T.shape, Y_F.shape, A.shape, B.shape, D.shape
    Y_T=np.nansum( Y_T, axis=0)
    Y_F=np.nansum(Y_F, axis=0)
    A=np.nansum(A,axis=0)
    B=np.nansum(B,axis=0)
    D=np.nansum(D,axis=0)
    #print "len(Y_T), len(Y_F), len(A), len(B), len(D)", len(Y_T), len(Y_F), len(A), len(B), len(D)
    print "Y_T.shape, Y_F.shape, A.shape, B.shape, D.shape", Y_T.shape, Y_F.shape, A.shape, B.shape, D.shape
    det=(1./(A*D-B*B))
    X_T=det*(D*Y_T - B*Y_F) 
    X_F=det*(-B*Y_T + A*Y_F)
    r =rsum1[0]
    print "len(X_T), len(X_F)", len(X_T), len(X_F) 
    print "rsum1[0]", rsum1[0]
    print "X_T", X_T
    print "X_F", X_F
    #### TEMPORAL> Create a covariance matrix with zeroes as a place holder. 
    cov = np.zeros( (nradial_bins, nradial_bins) )   ### How to calculate the errors?? 
    return r, X_T, X_F, cov, rcoeff

