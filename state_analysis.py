# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 15:20:20 2018

@author: valentin.metillon
"""


import qutip as qt
from numpy import sqrt


# |N0><0N| coherences
#Form of the density matrix :
#       0      1             Dimensions : Nh1*NH2
#     ______ ______
#     0 1 2  0 1 2
#  | [      |       ] 0 
# 0| [      |%      ] 1 
#  | [      |       ] 2 
#    ________________
#  | [      |       ] 0 
# 1| [      |       ] 1
#  | [      |       ] 2


def coherence_noon(rho,NH1,NH2,N):
    sigmax_n = qt.tensor(qt.fock(NH1,N),qt.fock(NH2,0))*qt.tensor(qt.fock(NH1,0).dag(),qt.fock(NH2,N).dag())\
    +qt.tensor(qt.fock(NH1,0),qt.fock(NH2,N))*qt.tensor(qt.fock(NH1,N).dag(),qt.fock(NH2,0).dag())             
    sigmay_n = 1j*qt.tensor(qt.fock(NH1,N),qt.fock(NH2,0))*qt.tensor(qt.fock(NH1,0).dag(),qt.fock(NH2,N).dag())\
    -1j*qt.tensor(qt.fock(NH1,0),qt.fock(NH2,N))*qt.tensor(qt.fock(NH1,N).dag(),qt.fock(NH2,0).dag())    
    sigmaz_n = qt.tensor(qt.fock(NH1,0),qt.fock(NH2,N))*qt.tensor(qt.fock(NH1,0).dag(),qt.fock(NH2,N).dag())\
    -qt.tensor(qt.fock(NH1,N),qt.fock(NH2,0))*qt.tensor(qt.fock(NH1,N).dag(),qt.fock(NH2,0).dag())       
    return sqrt((sigmax_n*rho).tr()**2+(sigmay_n*rho).tr()**2+(sigmaz_n*rho).tr()**2)
    
#NH1 = 3
#NH2 = 3 
#N = 2
#noon  = (qt.tensor(qt.fock(NH1,N),qt.fock(NH2,0))+qt.tensor(qt.fock(NH1,0),qt.fock(NH2,N)))/sqrt(2)
#rho_noon = qt.ket2dm(noon) 
#rho_2=.5*rho_noon+.5*qt.tensor(qt.qeye(NH1),qt.qeye(NH2))
#print(coherence_noon(rho_noon,NH1,NH2,N)    )
#print(coherence_noon(rho_2,NH1,NH2,N))


    