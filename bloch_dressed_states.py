# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 15:04:12 2018

@author: valentin.metillon
"""


import qutip as qt
from numpy import zeros


#NH = 3 


#psi_e = qt.basis(2,0)
#vac = qt.fock(NH,0)
#
#psi0 = qt.tensor(psi_e,vac)  

#Form of the density matrix :
#       e      g             Dimensions : 2*NH
#     ______ ______
#     0 1 2  0 1 2
#  | [      |       ] 0 
# e| [      |       ] 1 
#  | [      |       ] 2 
#    ________________
#  | [      |       ] 0 
# g| [      |       ] 1
#  | [      |       ] 2

def X_dressed(NH,n):
    X_ds = zeros((2*NH,2*NH))
    X_ds[n,NH+n+1] = 1
    X_ds[NH+n+1,n] = 1

    X_ds = qt.Qobj(X_ds)
    X_ds.dims = [[2,NH],[2,NH]]
    return(X_ds)

def Y_dressed(NH,n):
    Y_ds = zeros((2*NH,2*NH),dtype='complex')
    Y_ds[n,NH+n+1] = -1j
    Y_ds[NH+n+1,n] = 1j

    Y_ds = qt.Qobj(Y_ds)
    Y_ds.dims = [[2,NH],[2,NH]]
    return(Y_ds)
    
def Z_dressed(NH,n):    
    Z_ds = zeros((2*NH,2*NH))
    Z_ds[n,n] = 1
    Z_ds[NH+n+1,NH+n+1] = -1

    Z_ds = qt.Qobj(Z_ds)
    Z_ds.dims = [[2,NH],[2,NH]]
    return(Z_ds)
    
    

#X_ds = X_dressed(NH)
#Y_ds = X_dressed(NH)
#Z_ds = X_dressed(NH)

#
####Evolution and plot
#delta = 0.
#omega = 1.
#
#
#t = linspace(0,2*pi*5,100)
#
#H_jc = omega*X_ds + delta * Z_ds
#
#e_ops = [X_ds, Y_ds, Z_ds]
#
#opts = Options(store_states=False, store_final_state=True)
#res = qt.mesolve(H_jc,psi0,t,[],e_ops,options=opts)
#
#xyz = [res.expect[0],res.expect[1],res.expect[2]]
#
#b=qt.Bloch()
#b.add_points(xyz)
#b.show()
