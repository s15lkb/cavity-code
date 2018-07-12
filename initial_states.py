# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 17:59:28 2018

@author: valentin.metillon
"""

import qutip as qt

def initial_state(N_atoms,atomic_states,psi_c_1,psi_c_2):
    #psi_c_1 : state of cavity 1
    #psi_c_2 : state of cavity 1
    #atomic states can either be a list of 'e' and 'g' or just 'e' or 'g' if all atoms are in the same state
    if atomic_states=='g':
        psi_temp = qt.tensor(psi_c_1,psi_c_2)
        for i in range(N_atoms):
            psi_temp = qt.tensor(qt.basis(2,1),psi_temp)
    elif atomic_states=='e':
        psi_temp = qt.tensor(psi_c_1,psi_c_2)
        for i in range(N_atoms):
            psi_temp = qt.tensor(qt.basis(2,0),psi_temp)
        
    else:        
        psi_temp = qt.tensor(psi_c_1,psi_c_2)
        for i in range(N_atoms):
            if atomic_states[N_atoms-1-i]=='g':
                psi_temp = qt.tensor(qt.basis(2,1),psi_temp)
            if atomic_states[N_atoms-1-i]=='e':
                psi_temp = qt.tensor(qt.basis(2,0),psi_temp)
    #psi_temp_dm = qt.ket2dm(psi_temp)
    return psi_temp #psi_temp_dm
        