# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 16:19:13 2018

@author: valentin.metillon
"""
import qutip as qt
from numpy import sign, exp


def hamiltonians_atoms(NH1,NH2,N):
    #Constructs the hamiltonians for the 2-level structure of each atom
    #The result is given as a list of hamiltonians in the same order as in the tensor product
    hamiltonians = []
    H_temp = qt.tensor(qt.qeye(NH1),qt.qeye(NH2))
    for i in range(N):
        for j,H in enumerate(hamiltonians):
            hamiltonians[j]= qt.tensor(qt.qeye(2),hamiltonians[j])
        hamiltonians=[qt.tensor(.5*qt.sigmaz(),H_temp)]+hamiltonians
        H_temp = qt.tensor(qt.qeye(2),H_temp)
    return(hamiltonians)

def hamiltonians_n_atoms_C1(NH1,NH2,N):
    #Contructs the hamiltonians for the interaction of each atom with C1
    #The result is given as a list of hamiltonians in the same order as in the tensor product
    hamiltonians = [] 
    Hd_temp = qt.tensor(qt.destroy(NH1),qt.qeye(NH2))
    Hc_temp = qt.tensor(qt.create(NH1),qt.qeye(NH2)) 
    for i in range(N):
        for j,H in enumerate(hamiltonians):
            hamiltonians[j]= qt.tensor(qt.qeye(2),hamiltonians[j])
        hamiltonians=[.5*(qt.tensor(qt.sigmap(),Hd_temp)+qt.tensor(qt.sigmam(),Hc_temp))]+hamiltonians
        Hd_temp = qt.tensor(qt.qeye(2),Hd_temp)
        Hc_temp = qt.tensor(qt.qeye(2),Hc_temp)
    return(hamiltonians)


def hamiltonians_n_atoms_C2(NH1,NH2,N_atoms):
    #Contructs the Jaynes_cummings hamiltonians for the interaction of each atom with C2
    #The result is given as a list of hamiltonians in the same order as in the tensor product
    hamiltonians = [] 
    Hd_temp = qt.tensor(qt.qeye(NH1),qt.destroy(NH2))
    Hc_temp = qt.tensor(qt.qeye(NH1),qt.create(NH2)) 
    for i in range(N_atoms):
        for j,H in enumerate(hamiltonians):
            hamiltonians[j]= qt.tensor(qt.qeye(2),hamiltonians[j])
        hamiltonians=[.5*(qt.tensor(qt.sigmap(),Hd_temp)+qt.tensor(qt.sigmam(),Hc_temp))]+hamiltonians
        Hd_temp = qt.tensor(qt.qeye(2),Hd_temp)
        Hc_temp = qt.tensor(qt.qeye(2),Hc_temp)
    return(hamiltonians)
    
def local_operator(M,i,NH1,NH2,N_atoms):
    #builds the tensor product corresponding to a local operator M of system i
    #if i= N_atoms (resp. N_atoms+1), the operator acts on C1 (resp. C2)
    ops = []    
    for j in range(N_atoms):
        ops+=[qt.qeye(2)]
    ops+=[qt.qeye(NH1),qt.qeye(NH2)]
    ops[i] = M
    return(qt.tensor(ops))
        

def marche(t,t_c,T):
    return(.5*(1+sign(T/2-abs(t-t_c))))     
    
def killer(t,args):
    #killer ramp for rabi oscillation 
    #delta_k_i is the killer during the oscillation in cav i
    #T_i_i is the interacion time with cavity i
    #t_k_i is the time of the middle of the pulse
    t_k_1 = args['t_k_1']
    T_i_1 = args['T_i_1']
    delta_k_1 = args['delta_k_1']
    t_k_2 = args['t_k_2']
    T_i_2 = args['T_i_2']
    delta_k_2 = args['delta_k_2']
    
    
    delta = args['delta']
    
    
    return(delta+\
    (delta_k_1-delta)*marche(t,t_k_1,T_i_1)+\
    (delta_k_2-delta)*marche(t,t_k_2,T_i_2))
    
    
    
def f1(t,args):
    #Strength of the interaction of the atom with C1, with Gaussian mode structure
    #t_c_1 is the time at which the atom crosses the center of the cavity
    #t_mode_1 is the duration of the interaction, given by the width of the mode and the velocity of the atom
    t_c_1 = args['t_c_1']
    t_mode_1 = args['t_mode_1']
    omega_r_1 = args['omega_r_1']
    return(omega_r_1*exp(-(t_c_1-t)**2/t_mode_1**2))

def f2(t,args):
    #Strength of the interaction of the atom with C2, with Gaussian mode structure
    #t_c_2 is the time at which the atom crosses the center of the cavity
    #t_mode_2 is the duration of the interaction, given by the width of the mode and the velocity of the atom
    t_c_2 = args['t_c_2']
    t_mode_2 = args['t_mode_2']
    omega_r_2 = args['omega_r_2']
    return(omega_r_2*exp(-(t_c_2-t)**2/t_mode_2**2))
    

def build_hamiltonians(NH1,NH2,N_atoms,killer,f1,f2,args):
    H_jc = []
    Hs_int_1 = hamiltonians_n_atoms_C1(NH1,NH2,N_atoms)
    Hs_int_2 = hamiltonians_n_atoms_C2(NH1,NH2,N_atoms)
    Hs_at = hamiltonians_atoms(NH1,NH2,N_atoms)
    for i in range(N_atoms):
        H_jc+= [[Hs_at[i],killer], [Hs_int_1[i],f1], [Hs_int_2[i],f2]]
    return(H_jc)