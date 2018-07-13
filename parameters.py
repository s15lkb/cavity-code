# -*- coding: utf-8 -*-
"""
PARAMETERS

Author : Gautier Depambour

You will find in this file all usefull parameters for the simulation.
"""

""" NUMBER OF CAVITIES """
N_cav = 3

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
NHs = [3,3,3]
NHs = NHs[0:N_cav]

freq_unit = 50*10**3      #frequency unit
time_unit = 1/freq_unit    #time unit

NT = 101
t_min = -3.
t_max = 12.
dt = (t_max-t_min)/(NT-1)

###Physical parameters
t = linspace(t_min,t_max,NT)
t_real = time_unit*t

v=250.        #atom velocity

delta = 2*pi*67/50

#Cavity parameters
w = array([6*10**-3, 6*10**-3, 6*10**-3])  #mode width
omega_r = 2*pi*array([1.,1.,1.])

t_c = array([0.,5.,10.])   #Time at which the atom crosses the center of the mode
t_mode  = w/v/time_unit    #Characteristic time of the interaction between atom and cavity

#Parameters for the killer
t_k = t_c  #interaction squares are centered on the modes
T_i = array([3,3,8.42])*10**-6/time_unit  #Values in microseconds
delta_k = 2*pi/50*array([-32,-32,-32])


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
 'freq_unit':freq_unit, 'time_unit':time_unit, 'NT':NT, \
't_min':t_min, 't_max':t_max, 'dt':dt, 'w':w, 'v':v, \
'delta':delta, 'delta_k':delta_k,  'omega_r':omega_r, \
't':t, 't_real':t_real, 't_c':t_c,'t_mode':t_mode,  't_k':t_k, 'T_i':T_i, \
'epsilon':epsilon, 'eta_e':eta_e, 'eta_g':eta_g, 'NHs':NHs }

# Units
units = {'N_cav':'', 'mu':'', 'eta_prep_e':'', 'eta_prep_g':'', \
'freq_unit':'Hz', 'time_unit':'s', 'NT':'', \
't_min':t_min, 't_max':t_max, 'dt':dt, 'w':'m', 'v':'m/s', \
'delta':delta, 'delta_k':delta_k, \
't':t, 't_real':t_real, 't_c':t_c,'t_mode':t_mode, 'kappa_1':'', 'kappa_2':'', \
'n_th1':'', 'n_th2':'', 't_k':t_k, 'T_i':T_i,  \
'epsilon':'', 'eta_e':'', 'eta_g':'' }















