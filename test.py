# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 21:16:04 2018

@author: valentin.metillon
"""

import qutip as qt
from numpy import sqrt, exp, linspace
import matplotlib.pylab as plt

#psi = (qt.fock(2,0)+qt.fock(2,1))/sqrt(2)
#H = qt.sigmaz()
#args={}
#t = linspace(0,10,31)
#e_ops = qt.sigmax()
#res = qt.mesolve([[H, lambda t,args: exp(-t**2/10)]],psi,t,[],e_ops,args=args)
##res = qt.mesolve([[H, lambda t,args: 1]],psi,t,[],e_ops,args=args)
#
#plt.plot(t,res.expect[0]) 


def f_gen(i):
    return lambda x:i
    
a = []
for i in range(10):
    a+=[f_gen(i)]