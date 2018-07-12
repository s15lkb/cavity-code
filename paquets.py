# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 12:21:01 2018

@author: valentin.metillon
"""

from scipy.special import binom, factorial
from numpy import*
import matplotlib.pylab as plt
from matplotlib import cm

def poisson(mu,n):
    return(exp(-mu)*mu**n/factorial(n))



def P_n_at_1d(eps,mu,n):   #Probabilité d'avoir n atomes dans un paquet dans lequel on a détécté 1 atome.
    return(exp(-mu*(1-eps))*(mu*(1-eps))**(n-1)/factorial(n-1))

mu = linspace(0,2,101)
ep = linspace(0,1,100)    
n= array(range(1,10))
mus,eps, ns = meshgrid(mu,ep,n)   # [ i_ep, i_mu, i_n]

george = P_n_at_1d(eps,mus,ns)
#fig, axes = plt.subplots(2,2)

#nrm = plt.colors.Normalize(0, 1)
#cax,kw = plt.colorbar.make_axes([ax for ax in axes.flat])
#plt.colorbar(plt1, cax=cax, norm=nrm)
#
#plt1 = axes[0,0].pcolor(mu,ep,P_n_at_1d(eps,mus,ns)[:,:,0])
#plt.colorbar(plt1)
#
#plt1 = axes[0,1].pcolor(mu,ep,P_n_at_1d(eps,mus,ns)[:,:,1])
##plt.colorbar(plt1,ax=axes[0,1])
#
#plt1 = axes[1,0].pcolor(mu,ep,P_n_at_1d(eps,mus,ns)[:,:,2])
##plt.colorbar(plt1,ax=axes[1,0])
#
#plt1 = axes[1,1].pcolor(mu,ep,P_n_at_1d(eps,mus,ns)[:,:,3])
##plt.colorbar(plt1,ax=axes[1,1])
#
#plt.show()
