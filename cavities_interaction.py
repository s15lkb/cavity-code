""" 
ATOMS - CAVITIES INTERACTION

Authors : Gautier Depambour & Valentin Metillon

This code calculates the density matrix after interaction with both cavities.
It plots the number of photons in each cavity as well as the probabilities for the atoms to be in the excited state, 
as functions of time.

"""

#from numpy import *
import qutip as qt
import matplotlib.pylab as plt
from qutip.solver import Options
from hamiltonians import killer, f1, f2, build_hamiltonians, local_operator
from initial_states import initial_state
import sys


def DensityMatrixAfterInteraction(atomic_states, args, figs_interaction = True):
	# atomic_states must be a list of 'e' and 'g'
	#if args['N_cav'] == 0:
	
	if args['N_cav'] > 0:
		Numb1 = qt.tensor(qt.qeye(2),qt.qeye(2),qt.num(args['NH1']),qt.qeye(args['NH2']))
		Numb2 = qt.tensor(qt.qeye(2),qt.qeye(2),qt.qeye(args['NH1']),qt.num(args['NH2']))

#		a1 = qt.tensor(qt.destroy(args['NH1']),qt.qeye(args['NH2']))
#		a2 = qt.tensor(qt.qeye(args['NH1']),qt.destroy(args['NH2']))



#		Lm1 = sqrt(args['kappa_1']*(1+args['n_th1']))*a1    #Loss of a photon in C1
#		Lm2 = sqrt(args['kappa_2']*(1+args['n_th2']))*a2    #Loss of a photon in C2
#
#		Lp1 = sqrt(args['kappa_1']*(args['n_th1']))*a1.dag()    #Absorption of a photon in C1
#		Lp2 = sqrt(args['kappa_2']*(args['n_th2']))*a2.dag()    #Absorption of a photon in C2




	
		# Initial state
		N_atoms = len(atomic_states)
		psie = initial_state(N_atoms,atomic_states,qt.fock(args['NH1'],0),qt.fock(args['NH2'],0))
		H_jc = build_hamiltonians(args['NH1'],args['NH2'],N_atoms,killer,f1,f2, args)

		# Measurement operators
		Numb1 = local_operator(qt.num(args['NH1']),N_atoms,args['NH1'],args['NH2'],N_atoms)
		Numb2 = local_operator(qt.num(args['NH1']),N_atoms+1,args['NH1'],args['NH2'],N_atoms)
	
		ee = []
		for i in range(N_atoms):
			ee += [local_operator(qt.basis(2,0)*qt.basis(2,0).dag(),i,args['NH1'],args['NH2'],N_atoms) ]    # |e><e| For atom i
		e_ops = [Numb1,Numb2] + ee

		# Evolution
		opts = Options(store_states=True, store_final_state=True)
		res = qt.mesolve(H_jc,psie,args['t'],[],e_ops,args=args,options=opts)

		rho_after_interaction = res.states[-1] #.ptrace((N_atoms,N_atoms+1)) -> the partial trace is computed in detection.py
		#print (rho_after_interaction)
	
		if figs_interaction:
		# Plot the evolution of the number of photons in the cavities and the states of the atoms during interaction
			plt.plot(args['t_real'],res.expect[0],label='N1, '+str(N_atoms)+' atoms')
			plt.plot(args['t_real'],res.expect[1],label='N2, '+str(N_atoms)+' atoms')
			for i in range(len(ee)):
				plt.plot(args['t_real'],res.expect[2+i],label='Pe '+str(i))
			plt.xlabel("t_real")
			plt.title("Pour "+str(N_atoms)+" atomes")
			plt.legend()
			if N_atoms == 1:
				plt.suptitle("Etats de l'atome et des cavités en fonction du temps")
			if N_atoms > 1:
				plt.suptitle("Etats des atomes et des cavités en fonction du temps")
		
		# Plot the density matrix after interaction : state of the cavities + the atoms, since we have not yet taken the partial trace
			qt.hinton(qt.ket2dm(rho_after_interaction)) 
			plt.suptitle("State of the cavities + atoms")
			plt.show()
		
		#Return the mean number of photons in C1, the mean number of photons in C2, and the mean number of atoms is still excited
			print ("\n APRES INTERACTION AVEC LES DEUX CAVITES \n")
			print ("Le nombre moyen de photons dans la cavite 1 est "+str(res.expect[0][-1]))
			print ("Le nombre moyen de photons dans la cavite 2 est "+str(res.expect[1][-1]))
			print ("Le nombre moyen d'atomes encore excites est "+str(res.expect[2][-1]))
	
		return qt.ket2dm(rho_after_interaction)
	
	else:
		print ("Please specify a correct number of cavities")
		sys.exit()
	
	






