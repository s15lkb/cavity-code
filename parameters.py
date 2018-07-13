# -*- coding: utf-8 -*-
"""
PARAMETERS

Author : Gautier Depambour

You will find in this file all usefull parameters for the simulation.
"""

""" NUMBER OF CAVITIES """
N_cav = 1

""" PREPARATION """

from numpy import pi, linspace, array

# Specify the state of the prepared atoms : 'g' or 'e'
prepared_state = 'e'

# Mean number of atom(s) in each bunch
mu = 0.1

# Atom preparation error in excited state
eta_prep_e = 0.03

# Atom preparation error in ground state
eta_prep_g = 0


""" INTERACTION """

###Computation parameters
NH1 = 3  
NH2 = 3
NHs = [3,3,3]
NHs = NHs[0:N_cav]

freq_unit = 50*10**3      #frequency unit
time_unit = 1/freq_unit    #time unit

NT = 101
t_min = -3.
t_max = 6.
dt = (t_max-t_min)/(NT-1)

###Physical parameters
w = array([6*10**-3, 6*10**-3, 6*10**-3])  #mode width
v=250.        #atom velocity

delta = 2*pi*67/50

delta_k_1 = 2*pi*(-32)/50
delta_k_2 = 2*pi*(-32)/50

omega_r = 2*pi*array([1.,1.,1.])

#omega_r_1 = 2*pi*1.
#omega_r_2 = 2*pi*1.

t = linspace(t_min,t_max,NT)
t_real = time_unit*t


t_c = array([0.,5.,10.])
t_mode  = w/v/time_unit

#t_c_1 = 0   #center cav 1
#t_mode_1  = w/v/time_unit
#
#t_c_2 = 5.  #center cav 2
#t_mode_2  = w/v/time_unit


#kappa_1 = .1  #Loss rate in C1
#kappa_2 = .1  #Loss rate in C2
#
#n_th1 = 1     #Number of thermal photons in C1
#n_th2 = 1     #Number of thermal photons in C2

#Parameters for the killer
t_k_1 = t_c[0]
T_i_1 = 4.198*10**-6/time_unit
t_k_2 = t_c[1]
T_i_2 = 8.42*10**-6/time_unit


""" DETECTION """

# Detection efficiency
epsilon = 0.5

# Probability to detect 'g' whereas the atom was in state 'e'
eta_e = 0.06

# Probability to detect 'e' whereas the atom was in state 'g'
eta_g = 0.06

""" Store all parameters in a dictionary """

# Parameters
args = {'N_cav':N_cav, 'prepared_state':prepared_state, 'mu':mu, 'eta_prep_e':eta_prep_e, 'eta_prep_g':eta_prep_g, \
'NH1':NH1, 'NH2':NH2, 'freq_unit':freq_unit, 'time_unit':time_unit, 'NT':NT, \
't_min':t_min, 't_max':t_max, 'dt':dt, 'w':w, 'v':v, \
'delta':delta, 'delta_k_1':delta_k_1, 'delta_k_2':delta_k_2, 'omega_r':omega_r, \
't':t, 't_real':t_real, 't_c':t_c,'t_mode':t_mode,  't_k_1':t_k_1, 'T_i_1':T_i_1,  't_k_2':t_k_2, 'T_i_2':T_i_2, \
'epsilon':epsilon, 'eta_e':eta_e, 'eta_g':eta_g, 'NHs':NHs }

# Units
units = {'N_cav':'', 'mu':'', 'eta_prep_e':'', 'eta_prep_g':'', \
'freq_unit':'Hz', 'time_unit':'s', 'NT':'', \
't_min':t_min, 't_max':t_max, 'dt':dt, 'w':'m', 'v':'m/s', \
'delta':delta, 'delta_k_1':delta_k_1, 'delta_k_2':delta_k_2,  \
't':t, 't_real':t_real, 't_c':t_c,'t_mode':t_mode, 'kappa_1':'', 'kappa_2':'', \
'n_th1':'', 'n_th2':'', 't_k_1':t_k_1, 'T_i_1':T_i_1,  't_k_2':t_k_2, 'T_i_2':T_i_2, \
'epsilon':'', 'eta_e':'', 'eta_g':'' }















