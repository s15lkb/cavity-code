""" 
PREPARATION

Author : Gautier Depambour

This code calculates the probabilities to prepare a bunch with Ng atom(s) in the fundamental state 
and Ne atoms in the excited state, as a function of the mean number of atom(s) in each bunch, mu, and the 
preparation errors : eta_prep_g and eta_prep_e


"""
from general_functions import Poisson


def ProbaNumberRealAtoms(Nreal, args):
	# This function returns the probability of having Nreal atoms in a bunch, knowing mu = the mean number of atoms per bunch
	return Poisson(args['mu'], Nreal)


def PreparedStates(Nreal):
	# This function returns a list with all possible atomic states with Nreal atoms
	
	Atomic_states = ['e','g']
	Possible_states = {} # This dictionary allows to store intermediate list of states
	if Nreal == 0:
		return []
	else:
		Possible_states[0] = [[]]
		k = 0
		for i in range(Nreal):   # iteration atom by atom
			k += 1
			Possible_states[k] = []
			for element in Possible_states[k-1]:
				for state in Atomic_states:
					Possible_states[k] += [element + [state]]    # Add one atom in a specific state to the previous possible states of the other atoms
	
		return Possible_states[k]    # keep the last list of possible states, when all real atoms have been considered

def ProbaPreparedStates(Nreal, args):
	# This function returns a dictionary with all possible prepared states of Nreal atoms with their respective probabilities
	Prepared_states = PreparedStates(Nreal)
	Prepared_states_with_probas = {}
	for state in Prepared_states:
		proba = 1
		if args['prepared_state'] == 'g':
			for i in range(len(state)):
				if state[i] == 'g':
					proba = proba * (1 - args['eta_prep_g'])
				elif state[i] == 'e':
					proba = proba * args['eta_prep_g']
		elif args['prepared_state'] == 'e':
			for i in range(len(state)):
				if state[i] == 'g':
					proba = proba * args['eta_prep_e']
				elif state[i] == 'e':
					proba = proba * (1 - args['eta_prep_e'])
		Prepared_states_with_probas[tuple(state)] = proba
		
	return Prepared_states_with_probas
	




