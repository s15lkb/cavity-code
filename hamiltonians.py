# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 16:19:13 2018

@author: valentin.metillon
"""
import qutip as qt
from numpy import sign, exp

def hamiltonians_atoms(NHs,N):
    #Constructs the hamiltonians for the 2-level structure of each atom
    #The result is given as a list of hamiltonians in the same order as in the tensor product
    hamiltonians = []
    H_temp = qt.qeye(NHs[0])
    for NH in NHs[1:]:
        H_temp = qt.tensor(H_temp, qt.qeye(NH))
    for i in range(N):
        for j,H in enumerate(hamiltonians):
            hamiltonians[j]= qt.tensor(qt.qeye(2),hamiltonians[j])
        hamiltonians=[qt.tensor(.5*qt.sigmaz(),H_temp)]+hamiltonians
        H_temp = qt.tensor(qt.qeye(2),H_temp)
    return(hamiltonians)

def hamiltonians_n_atoms_Ci(NHs,N,i):
    #Contructs the hamiltonians for the interaction of each atom with cavity i
    #The result is given as a list of hamiltonians in the same order as in the tensor product
    hamiltonians = [] 
    H_cav = []
    for NH in NHs:
        H_cav+=[qt.qeye(NH)]
    H_cav_d = H_cav
    H_cav_d[i] = qt.destroy(NHs[i])
    Hd_temp = qt.tensor(H_cav_d)
    H_cav_c = H_cav
    H_cav_c[i] = qt.create(NHs[i])
    Hc_temp = qt.tensor(H_cav_c)

    
    for i in range(N):
        for j,H in enumerate(hamiltonians):
            hamiltonians[j]= qt.tensor(qt.qeye(2),hamiltonians[j])
        hamiltonians=[.5*(qt.tensor(qt.sigmap(),Hd_temp)+qt.tensor(qt.sigmam(),Hc_temp))]+hamiltonians
        Hd_temp = qt.tensor(qt.qeye(2),Hd_temp)
        Hc_temp = qt.tensor(qt.qeye(2),Hc_temp)
    return(hamiltonians)
    
def local_operator_atom(M,i,NHs,N_atoms):
    #builds the tensor product corresponding to a local operator M of atom i
    ops = []    
    for j in range(N_atoms):
        ops+=[qt.qeye(2)]
    for NH in NHs:
        ops+=[qt.qeye(NH)]
    ops[i] = M
    return(qt.tensor(ops))
        
def local_operator_cavity(M,i,NHs,N_atoms):
    #builds the tensor product corresponding to a local operator M of atom i
    ops = []    
    for j in range(N_atoms):
        ops+=[qt.qeye(2)]
    for NH in NHs:
        ops+=[qt.qeye(NH)]
    ops[N_atoms+i] = M
    return(qt.tensor(ops))

def marche(t,t_c,T):
    return(.5*(1+sign(T/2-abs(t-t_c))))     
    
def killer(t,args):
    #Killer ramp for rabi oscillation 
    #delta_k_i is the killer during the oscillation in cav i
    #T_i_i is the interacion time with cavity i
    #t_k_i is the time of the middle of the pulse

    delta = args['delta']
    K = delta
    for i in range(args['N_cav']):
        t_k = args['t_k'][i]
        T_i = args['T_i'][i]
        delta_k = args['delta_k'][i]
        K+=(delta_k-delta)*marche(t,t_k,T_i)
   
    return K
        
def f_i(args,i):
    #Returns a function giving the
    #strength of the interaction of the atom with cavity i, with Gaussian mode structure
    #t_c is the time at which the atom crosses the center of the cavity
    #t_mode is the duration of the interaction, given by the width of the mode and the velocity of the atom
    #omega_r is the Rabi frequency
    t_c = args['t_c'][i]
    t_mode = args['t_mode'][i]
    omega_r = args['omega_r'][i]
    return(lambda t, args: omega_r*exp(-(t_c-t)**2/t_mode**2))

def build_hamiltonians(NHs,N_atoms,killer,args):
    H_jc = []
    
    Hs_at = hamiltonians_atoms(NHs,N_atoms)
    for j in range(N_atoms):
        H_jc+= [[Hs_at[j],killer]]
    for i in range(len(NHs)):
        for j in range(N_atoms):
            Hs_int_i = hamiltonians_n_atoms_Ci(NHs,N_atoms,i)
            H_jc+= [[Hs_int_i[j],f_i(args,i)]]
    return(H_jc)