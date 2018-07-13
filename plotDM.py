# -*- coding: utf-8 -*-
"""
FIDELITY WITHOUT CONTRAST

Author : Gautier Depambour

This code allows to plot several density matrices corresponding to different values of one varying parameter.
Then, it calculates the best fidelity between the ideal state |10 > + exp(i*phi)|01 > and the simulated density matrix
(it finds the value of phi which maximises the fidelity), for each value of the varying parameter.
Finally, it shows the fidelity as a function of the varying parameter.


"""
import matplotlib.pyplot as plt
from matplotlib import rc
import qutip as qt
import numpy as np
import sys

from parameters import args
from parameters import units
from state_prediction import RhoFinal

""" FILL THIS PART OF THE CODE """

# Number of plots 
N_plots = 11

# Varying parameter 
varying_param = 'eta_g'

# Range of the varying parameter (list)
Range = np.linspace(0,1,N_plots)
#Range =[0.1]

# Detection sequence (post-selection)
Detection = 'g'

# Maximum number of real atoms -> must be at least equal to the number of detected atoms
Nrealmax = 2

# Method to calculate fidelity (see below) : 1 or 2
method_fidelity = 1
# If method 2 chosen, specify the the value of the varying parameter for which you want the best fidelity (this is going to define phi)
specific_value = 0.3
if not (Range[0] < specific_value < Range[-1] or Range[-1] < specific_value < Range[0]):
	print ("The value must be in the range of the varying parameter, ie between "+str(Range[0])+" and "+str(Range[-1])+" !")
	sys.exit()


""" COMPUTATION OF DENSITY MATRICES """

# Warnings
if len(Range) < N_plots:
	print ("WARNING : there must be more values for the varying parameter than the number of plots !")
	N_plots = len(Range)
	print ("By default N_plots = "+str(N_plots))
if Nrealmax < len(Detection):
	print("WARNING : it is impossible to have more detected than real atoms !")
	sys.exit()
	
# IMPORTANT WARNING : some parameters can not be varying parameters, because they are calculated (they only depend on other parameters)
warning = "The varying parameter can not be "+str(varying_param)+" because it is a calculated quantity (it depends on other parameters)"
if varying_param == 'freq_unit':
	print (warning)
	sys.exit()
if varying_param == 't_real':
	print (warning)
	sys.exit()
if varying_param == 'dt' or varying_param == 't':
	print (warning+"\nInstead of "+str(varying_param)+", choose the varying parameter among t_min, t_max or NT")
	sys.exit()
if varying_param == 't_mode_1' or varying_param == 't_mode_2':
	print (warning+"\nInstead of "+str(varying_param)+", choose the varying parameter among w or v")
	sys.exit()


# Compute density matrices
Density_matrices = []
step = int(len(Range) / N_plots)
print ("Values for "+str(varying_param))
for i in range(N_plots):
	args[varying_param] = Range[i * step]
	print (float(int(round(args[varying_param]*1000))/1000))
	# Recalculate parameters which can depend on the varying parameter
	args['freq_unit'] = 1/args['time_unit']
	args['dt'] = (args['t_max']-args['t_min'])/(args['NT']-1)
	args['t'] = np.linspace(args['t_min'],args['t_max'],args['NT'])
	args['t_real'] = args['time_unit']*args['t']
	args['t_mode_1'] = args['w']/args['v']/args['time_unit']
	args['t_mode_2'] = args['w']/args['v']/args['time_unit']
	if varying_param == 'time_unit':
		args['T_i_1'] = 4.198*10**-6/args['time_unit']
		args['T_i_2'] = 8.42*10**-6/args['time_unit']
	if varying_param == 't_c_1':
		args['t_k_1'] = args['t_c_1']
	if varying_param == 't_c_2':
		args['t_k_2'] = args['t_c_2']
	if varying_param == 't_k_1':
		args['t_c_1'] = args['t_k_1']
	if varying_param == 't_k_2':
		args['t_c_2'] = args['t_k_2']
	# Compute density matrices for each value of the varying parameter
	Density_matrices += [RhoFinal(Nrealmax, Detection, args, fig_rhofinal=False, fig_qm=False)]

""" FIDELITY """

# Calculate the fidelity between the ideal state |10 > + exp(i*phi) |01 > and the simulated matrices
# Two methods : 
#  1 - either phi is optimized in each case, and the best fidelity is stored in a dictionary
#  2 - or we fix phi to have the best fidelity for one value of the varying parameter, and we keep this phi to calculate the fidelity for the other values taken by the varying parameter

phi = np.linspace(0,2*np.pi,100)
Psi_ideal = {}
DM_ideal = {}
for i in range(len(phi)):
	Psi_ideal[i] = (qt.tensor(qt.fock(3,0),qt.fock(3,1)) + np.exp(1j * phi[i]) * qt.tensor(qt.fock(3,1),qt.fock(3,0)))/np.sqrt(2)
	DM_ideal[i] = qt.ket2dm(Psi_ideal[i])
Fidelity = {}

	# First method
if method_fidelity == 1:
	Max_fidelity = {}
	Max_index = {}
	for i in range(N_plots):
		Fidelity[i] = {}
		Max_fidelity[i] = 0
		Max_index[i] = 0
		for j in range(len(phi)):
			Fidelity[i][j] = (qt.fidelity(DM_ideal[j],Density_matrices[i]))**2   # the qt.fidelity function returns float(np.real((A * (B * A)).sqrtm().tr()))
			if Fidelity[i][j] > Max_fidelity[i]:
				Max_fidelity[i] = Fidelity[i][j]
				Max_index[i] = j
	
		print ("\nFor "+str(varying_param)+" = "+str(float(int(round(Range[i * step]*1000))/1000))+", the best fidelity is "+str(float(int(round(Max_fidelity[i]*1000))/1000))+" for phi = "+str(float(int(phi[Max_index[i]]*100)/100)))

	# Second method
elif method_fidelity == 2:
	args[varying_param] = specific_value
	print ("The specific value of the varying parameter for which the fidelity is maximized is "+str(specific_value))
	Density_matrix_specific_value = RhoFinal(Nrealmax, Detection, args, fig_rhofinal=False, fig_qm=False)
	Max_fidelity = 0
	Max_index = 0
	Fidelity_specific_value = {}
	for i in range(len(phi)):
		Fidelity_specific_value[i] = (qt.fidelity(DM_ideal[i],Density_matrix_specific_value))**2
		print (Fidelity_specific_value[i])
		if Fidelity_specific_value[i] > Max_fidelity:
			Max_fidelity = Fidelity_specific_value[i]
			Max_index = i
			print (Max_index)
	print ("\nFor "+str(varying_param)+" = "+str(float(int(round(specific_value*1000))/1000))+", the best fidelity is "+str(float(int(round(Max_fidelity*1000))/1000))+" for phi = "+str(float(int(phi[Max_index]*100)/100))+"\n")
	for i in range(N_plots):
		Fidelity[i] = (qt.fidelity(DM_ideal[Max_index],Density_matrices[i]))**2   # the qt.fidelity function returns float(np.real((A * (B * A)).sqrtm().tr()))
		print ("\nFor "+str(varying_param)+" = "+str(float(int(round(Range[i * step]*1000))/1000))+" and phi = "+str(float(int(phi[Max_index]*100)/100))+", the fidelity is "+str(float(int(round(Fidelity[i]*1000))/1000)))


else:
	print ("Please specify a correct argument for method_fidelity")
	sys.exit()


""" PLOTS """

# Store diagonal coefficients in a dictionary and print them
Coeffs = {}
Cavity_states = [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)]
for i in range(N_plots):
	print ("\nDensity matrix for "+ str(varying_param)+" = "+str(Range[i * step])+"\n")
	Coeffs[i] = []
	for j in range(0,9):
		Coeffs[i] += [Density_matrices[i][j,j]]
		print ("coeff "+str(Cavity_states[j])+" = "+str(Coeffs[i][j]))
		

# GLOBAL PLOT

# Plot all density matrices
for i in reversed(range(0,N_plots)):
	qt.hinton(Density_matrices[i])
	plt.suptitle("Density matrix with "+str(varying_param)+" = "+str(Range[i * step])+", for detection : "+str(Detection)+", without contrast")
	if args['NH1'] == 3 and args['NH2'] == 3:
		plt.text(-28.9,1.02,np.float(int(round(Density_matrices[i][0,0].real*100))/100))
		plt.text(-26.1,0.885,np.float(int(round(Density_matrices[i][1,1].real*100))/100))
		plt.text(-23.3,0.75,np.float(int(round(Density_matrices[i][2,2].real*100))/100))
		plt.text(-20.5,0.615,np.float(int(round(Density_matrices[i][3,3].real*100))/100))
		plt.text(-17.7,0.48,np.float(int(round(Density_matrices[i][4,4].real*100))/100))
		plt.text(-15,0.345,np.float(int(round(Density_matrices[i][5,5].real*100))/100))
		plt.text(-12.2,0.21,np.float(int(round(Density_matrices[i][6,6].real*100))/100))
		plt.text(-9.4,0.075,np.float(int(round(Density_matrices[i][7,7].real*100))/100))
		plt.text(-6.6,-0.06,np.float(int(round(Density_matrices[i][8,8].real*100))/100))
		plt.text(-20.9,0.885,np.float(int(round(Density_matrices[i][1,3].real*1000))/1000))
		plt.text(-26.4,0.615,np.float(int(round(Density_matrices[i][3,1].real*1000))/1000))
	if method_fidelity == 1:
		plt.text(0,1.1,r"$\Psi_{ideal} = |10> + e^{i\phi}|01>$"+"\nFidelity max = "+str(float(int(Max_fidelity[i]*100)/100))+"\nfor $\phi$ = "+str(float(int(phi[Max_index[i]]*100)/100)))
	if method_fidelity == 2:
		plt.text(0,1.1,r"$\Psi_{ideal} = |10> + e^{i\phi}|01>$"+"\nFidelity max = "+str(float(int(Max_fidelity*100)/100))+"\nfor $\phi$ = "+str(float(int(phi[Max_index]*100)/100)))
	plt.text(-30,-0.2,r"prepared state = "+str(args['prepared_state'])+", $\eta_{prep_g}$ = "+str(args['eta_prep_g'])+", $\eta_{prep_e}$ = "+str(args['eta_prep_e'])+", $\mu$ = "+str(args['mu'])+", $\epsilon$ = "+str(args['epsilon'])+", $\eta_{g}$ = "+str(args['eta_g'])+", $\eta_{e}$ = "+str(args['eta_e']))

# Plot Fidelity as a function of the varying parameter
fig = plt.figure(figsize = (10,6))
x = []
y = []
for i in range(N_plots):
	x += [Range[i * step]]
	if method_fidelity == 1:
		y += [Max_fidelity[i]]
	if method_fidelity == 2:
		y += [Fidelity[i]]
title = r"Fidelity as a function of "+str(varying_param)+", for detection : "+str(Detection)+", without contrast\nwith "
Parameters = ['prepared_state','eta_prep_g','eta_prep_e','mu','epsilon','eta_g','eta_e']
LateX_parameters = ['prepared state','$\eta_{prep_g}$','$\eta_{prep_e}$','$\mu$','$\epsilon$','$\eta_{g}$','$\eta_{e}$']
for parameter in Parameters:
	if parameter == varying_param:
		continue
	elif Parameters.index(parameter) == len(Parameters)-1:
		title += "and "+str(LateX_parameters[Parameters.index(parameter)])+" = "+str(args[parameter])
	else:
		title += str(LateX_parameters[Parameters.index(parameter)])+" = "+str(args[parameter])+", "
if method_fidelity == 1:
	title += "\nPhi is optimized to get the best fidelity for each point"
if method_fidelity == 2:
	title += "\nPhi = "+str(float(int(round(phi[Max_index]*100))/100))+" to optimize the fidelity for "+str(varying_param)+" = "+str(float(int(round(specific_value*100))/100))


if str(units[varying_param]) == '':
	plt.xlabel(str(varying_param))
else:
	plt.xlabel(str(varying_param)+ ' (' + str(units[varying_param])+')')
plt.ylabel("Fidelity")
plt.scatter(x,y)
plt.title(title)

plt.show()











"""

I tried to put every density matrices in the same figure...


N_columns = 3
N_lines = int(N_plots / N_columns) 

N_subplot = []
for i in range(1,N_columns * N_lines + 1):
	N_subplot += [100 * N_lines + 10 * N_columns + i]	

fig = plt.figure()
for i in range(len(N_subplot)):
	print (N_subplot[i])
	fig.add_subplot(N_subplot[i])
	plt.subplot(N_subplot[i]) = plots[i]
"""	






