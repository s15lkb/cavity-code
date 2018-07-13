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
from hamiltonians import killer, build_hamiltonians, local_operator_cavity, local_operator_atom
from initial_states import initial_state
import sys


def DensityMatrixAfterInteraction(atomic_states, args, figs_interaction = True):
    # atomic_states must be a list of 'e' and 'g'
    #if args['N_cav'] == 0:
    
    if args['N_cav'] > 0:
    
        # Initial state
        N_atoms = len(atomic_states)
        NHs = args['NHs']
        cavities_states = []
        for NH in NHs:
            cavities_states+= [qt.fock(NH,0)]
        psie = initial_state(N_atoms,atomic_states,cavities_states)
        H_jc = build_hamiltonians(NHs,N_atoms,killer, args)
        
        Numb = []
        # Measurement operators
        for i in range(len(NHs)):
            Numb+= [local_operator_cavity(qt.num(NHs[i]),i,NHs,N_atoms)]

    
        ee = []
        for i in range(N_atoms):
            ee += [local_operator_atom(qt.basis(2,0)*qt.basis(2,0).dag(),i,NHs,N_atoms) ]    # |e><e| For atom i
        e_ops = Numb + ee

        # Evolution
        opts = Options(store_states=True, store_final_state=True)
        res = qt.mesolve(H_jc,psie,args['t'],[],e_ops,args=args,options=opts)

        rho_after_interaction = res.states[-1] #.ptrace((N_atoms,N_atoms+1)) -> the partial trace is computed in detection.py
        #print (rho_after_interaction)
    
        if figs_interaction:
        # Plot the evolution of the number of photons in the cavities and the states of the atoms during interaction
            plt.plot(args['t_real'],res.expect[0],label='N1, '+str(N_atoms)+' atoms')
            plt.plot(args['t_real'],res.expect[1],label='N2, '+str(N_atoms)+' atoms')
            plt.plot(args['t_real'],res.expect[2],label='N3, '+str(N_atoms)+' atoms')
            #plt.plot(args['t_real'],H_jc[],label='N2, '+str(N_atoms)+' atoms')
            for i in range(len(ee)):
                plt.plot(args['t_real'],res.expect[args['N_cav']+i],label='Pe '+str(i))
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
    
 
#
#from parameters import args
##from numpy import exp
#atomic_states = 'ee'
#a=DensityMatrixAfterInteraction(atomic_states, args, figs_interaction = True)
#
#
#
#
#N_atoms = len(atomic_states)
#NHs = args['NHs']
#cavities_states = []
#for NH in NHs:
#    cavities_states+= [qt.fock(NH,0)]
#psie = initial_state(N_atoms,atomic_states,cavities_states)
#H_jc = build_hamiltonians(NHs,N_atoms,killer, args)
#
#for i in range(len(H_jc)):
#    plt.plot(args['t_real'],H_jc[i][1](args['t'],args))

#t_c = args['t_c'][0]
#t_mode = args['t_mode'][0]
#omega_r = args['omega_r'][0]
#
#a = lambda t, args: args['omega_r'][0]*exp(-(args['t_c'][0]-t)**2/args['t_mode'][0]**2)
#
#t_c = args['t_c'][1]
#t_mode = args['t_mode'][1]
#omega_r = args['omega_r'][1]
#b = lambda t, args: args['omega_r'][1]*exp(-(args['t_c'][1]-t)**2/args['t_mode'][1]**2)