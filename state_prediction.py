""" 
STATE PREDICTION

Author : Gautier Depambour

This code (more precisely the RhoFinal function) calculates and plots the density matrix corresponding to the state in the cavities,
as a function of all the parameters specified in the file parameter.py (called here with "args"), 
for a specific sequence of detection (post-selection), 
and with a number of real atoms going from 0 to Nrealmax 
(normally it should be from 0 to infinity, but we don't have enough time to compute it :-) ).

"""

import matplotlib.pyplot as plt
import numpy as np
import qutip as qt


from preparation import PreparedStates, ProbaPreparedStates, ProbaNumberRealAtoms
from cavities_interaction import DensityMatrixAfterInteraction
from detection import Partial_trace, ProbasPossibleDetections, SelectDetection

debug = False

def Quantum_map(Nreal, Detection, args, fig_qm=True):
	# This function returns the quantum map for a given number of real atoms Nreal
	# 'Detection' is the sequence of detected atoms (post-selection)
	
	# Return the null matrix if the number of real atoms is 0 or if there are more detected than real atoms
	Null_vector1 = qt.Qobj([[0] * args['NH1']])
	Null_vector2 = qt.Qobj([[0] * args['NH2']])
	Null_matrix = qt.tensor(qt.ket2dm(Null_vector1),qt.ket2dm(Null_vector2))
	if Nreal == 0 or len(Detection) > Nreal:
		return Null_matrix

	# Store in a list the possible real states
	Possible_real_states = PreparedStates(Nreal)
	
	# Store in a dictionary the possible real states with their respective probabilities
	Possible_real_states_with_probas = ProbaPreparedStates(Nreal, args)
	
	# Compute all density matrices after interaction for each possible real state, and store them in a dictionary
	Rho_after_interaction = {}
	for state in Possible_real_states:
		Rho_after_interaction[tuple(state)] = DensityMatrixAfterInteraction(state, args, figs_interaction = False)
		
	# Compute all possible partial traces, for each density matrix after interaction, and for each possible ideal detection (two-level dictionary)
	All_Ptraces = {}
	for state_dm in Possible_real_states:
		All_Ptraces[tuple(state_dm)] = {}
		for state_ideal_detection in Possible_real_states:
			All_Ptraces[tuple(state_dm)][tuple(state_ideal_detection)] = Partial_trace(Rho_after_interaction[tuple(state_dm)], state_ideal_detection,args) 
	
	# Multiply each Ptrace by the probability of having the considered state before interaction
	All_Ptraces_with_weights = {}
	for state_dm in Possible_real_states:
		All_Ptraces_with_weights[tuple(state_dm)] = {}
		for state_ideal_detection in Possible_real_states:
			All_Ptraces_with_weights[tuple(state_dm)][tuple(state_ideal_detection)] = Possible_real_states_with_probas[tuple(state_dm)] * All_Ptraces[tuple(state_dm)][tuple(state_ideal_detection)] 
			
	# Compute a dictionary of all possible detections with their respective probabilities
	All_Possible_detections_with_probas = {}
	for state_ideal_detection in Possible_real_states:
		All_Possible_detections_with_probas[tuple(state_ideal_detection)] = ProbasPossibleDetections(state_ideal_detection, args)
	
	# Select detections compatible with the 'Detection' sequence given as an argument
	# For example, if the post-selection is 'g' when there are 2 real atoms, the possible detections are '0g' and 'g0', 
	# where '0' means : the real atom was not detected
	Selected_detections = SelectDetection(Nreal, Detection, args)
	
	# Store every possible partial traces with weights (probabilities) corresponding to the Detection
	Ptraces = []
	for state_dm in Possible_real_states:
		for state_ideal_detection in Possible_real_states:
			for real_detection in Selected_detections:
				Ptraces += [ All_Ptraces_with_weights[tuple(state_dm)][tuple(state_ideal_detection)] * All_Possible_detections_with_probas[tuple(state_ideal_detection)][tuple(real_detection)] ]
			
	if debug:
		print ("\n NUMBER OF REAL ATOMS = "+str(Nreal)+"\n")
		print ("\n POSSIBLE REAL STATES \n")
		print (Possible_real_states)
		print ("\n POSSIBLE REAL STATES WITH PROBABILITIES \n")
		print (Possible_real_states_with_probas)
		print ("\n DENSITY MATRICES AFTER INTERACTIONS \n")
		print (Rho_after_interaction)
		print ("\n POSSIBLE DETECTIONS \n")
		print (Selected_detections)
		print ("\n PARTIAL TRACES \n")
		print (All_Ptraces)
		print ("\n PARTIAL TRACES NORMALIZED \n")
		print (All_Ptraces_with_weights)

	# Calculate the quantum map by adding all these partial traces and dividing by the global trace
	quantum_map = Ptraces[0]
	for i in range(1,len(Ptraces)):
		quantum_map += Ptraces[i]

	quantum_map = quantum_map / quantum_map.tr()
	
	# Plot
	if fig_qm:
		qt.hinton(quantum_map) 
		plt.suptitle("State of the cavities")
		plt.show()
		
	return quantum_map 
	

def RhoFinal(Nrealmax, Detection, args, fig_rhofinal=True, fig_qm=False):
	# This function returns the density matrix of the cavities, computed with a number of real atoms from 0 to Nrealmax
	# 'Detection' is the list of detected atoms (post-selection)
	
	# Initialization : if Nreal = 0, RhoFinal is the null matrix (0 photon in the cavities)
	Rho_final = Quantum_map(0, Detection, args, fig_qm) * ProbaNumberRealAtoms(0, args)
	
	# Compute the final density matrix by considering all possible numbers of real atoms
	if Nrealmax > 0:
		for Nreal in range(1,Nrealmax+1):
			Rho_final += Quantum_map(Nreal, Detection, args, fig_qm) * ProbaNumberRealAtoms(Nreal, args)
	
	# Normalization  
	N_detected = len(Detection)
	normalisation_parameter = 0
	for i in range(N_detected): # We renormalize by considering the cases when there are less real than detected atoms (absurd cases)
		normalisation_parameter += ProbaNumberRealAtoms(i, args)
	Rho_final = Rho_final / (1 - normalisation_parameter)
	
	# Check that the trace of Rho_final is equal to 1
	trace = 0
	for i in range(args['NH1'] * args['NH2']):
		trace += Rho_final[i,i]
	if 1.0 - trace > 0.01:
		print ("\n The trace of the density matrix is "+str(np.float(int(round(trace.real*100))/100)))
		print ("WARNING : the trace should be equal to 1 !! Take a bigger number of atoms (compared to mu)")
		

	# Plot
	if fig_rhofinal:
		qt.hinton(Rho_final) 
		plt.suptitle("State of the cavities for detection "+str(Detection))
		if args['NH1'] == 3 and args['NH2'] == 3:
			plt.text(-28.8,1.02,np.float(int(round(Rho_final[0,0].real*1000))/1000))
			plt.text(-26.3,0.885,np.float(int(round(Rho_final[1,1].real*1000))/1000))
			plt.text(-23.5,0.75,np.float(int(round(Rho_final[2,2].real*1000))/1000))
			plt.text(-20.7,0.615,np.float(int(round(Rho_final[3,3].real*1000))/1000))
			plt.text(-17.9,0.48,np.float(int(round(Rho_final[4,4].real*1000))/1000))
			plt.text(-15.2,0.345,np.float(int(round(Rho_final[5,5].real*1000))/1000))
			plt.text(-12.4,0.21,np.float(int(round(Rho_final[6,6].real*1000))/1000))
			plt.text(-9.6,0.075,np.float(int(round(Rho_final[7,7].real*1000))/1000))
			plt.text(-6.8,-0.06,np.float(int(round(Rho_final[8,8].real*1000))/1000))
			plt.text(-20.9,0.885,np.float(int(round(Rho_final[1,3].real*1000))/1000))
			plt.text(-26.4,0.615,np.float(int(round(Rho_final[3,1].real*1000))/1000))
		plt.text(-30,-0.2,r"prepared state = "+str(args['prepared_state'])+", $\eta_{prep_g}$ = " \
+str(args['eta_prep_g'])+", $\eta_{prep_e}$ = "+str(args['eta_prep_e'])+", $\mu$ = "+str(args['mu']) \
+", $\epsilon$ = "+str(args['epsilon'])+", $\eta_{g}$ = "+str(args['eta_g'])+", $\eta_{e}$ = "+str(args['eta_e']))
		plt.show()
		
	return Rho_final







	


