""" 
DETECTOR IMPERFECTIONS

Author : Gautier Depambour

This codes models detector imperfections by taking into account :
- the detector efficiency (corresponding to the number of detected atoms vs the number of real atoms) : epsilon
- the detection errors (on atom states) : eta_e & eta_g
The function Detector_imperfections is based on the table p.58 in Xingxing ZHOU's thesis.

This code also calculates possible detections (with their probabilities) given a specific number of real atoms,
and several functions usefull for the computation of quantum_maps, done in the file state_prediction.py .

"""


#from matplotlib import *
import matplotlib.pyplot as plt
import qutip as qt
import sys

def Detector_imperfections(Ideal_detection, Detection,args):
	# This function returns the probability of observing the sequence of detected atoms (Detection) knowing the sequence of real atoms (Real_atoms).
	# Real_atoms = list of 'g' or 'e' or '0' (no atom)
	# Detection = list of 'g' or 'e' or '0' (no detection)
	# Detected_atom = epsilon
	# Non_detected_atom = 1 - epsilon
	# Atom in e detected e = 1 - eta_e
	# Atom in g detected g = 1 - eta_g
	# Atom in e detected g = eta_e
	# Atom in g detected e = eta_g
	
	N_total_real = len(Ideal_detection)
	N_total_detected = len(Detection)

	#WARNING : 
	if N_total_real != N_total_detected:
		print ("WARNING : the lists of real and detected atoms must have the same size ! (Put 0 if the atom is real but not detected)")
		sys.exit()

	#Number of existing atoms, and number of atoms detected in g or e (but NOT '0')
	N_existing_atoms = 0
	N_detected_only_ge = 0
	for i in range(len(Detection)):
		if Detection[i] == 'e' or Detection[i] == 'g':
			N_detected_only_ge += 1
	for i in range(len(Ideal_detection)):
		if Ideal_detection[i] == 'e' or Ideal_detection[i] == 'g':
			N_existing_atoms += 1
	
	#Initialization
	proba = 1
	
	#Detection_efficiency
	for i in range(len(Ideal_detection)):
		if Ideal_detection[i] == '0' and (Detection[i] == 'g' or Detection[i] == 'e'):
			return 0 # It's impossible to detect an atom which does not exist !
	
	if N_existing_atoms == 0 and N_detected_only_ge == 0:
		return 1
	elif N_existing_atoms < N_detected_only_ge: 
		return 0 # It's impossible to detect more atoms than real atoms !
	else:
		N_diff = N_existing_atoms - N_detected_only_ge
		proba = proba * (args['epsilon'] ** N_detected_only_ge) * ((1 - args['epsilon']) ** N_diff)
		
	#Detection errors
	for i in range(len(Ideal_detection)):
		if Ideal_detection[i] == 'e' and Detection[i] == 'e':
			proba = proba * (1 - args['eta_e'])
		if Ideal_detection[i] == 'g' and Detection[i] == 'g':
			proba = proba * (1 - args['eta_g'])
		if Ideal_detection[i] == 'e' and Detection[i] == 'g':
			proba = proba * args['eta_e']
		if Ideal_detection[i] == 'g' and Detection[i] == 'e':
			proba = proba * args['eta_g']

	return proba


def PossibleDetections(Nreal,args):
	# This function returns a list with every possible detection.
	# Nreal is the number of real atoms.
	
	# Create the list with all possible detection results (so the length of this list is equal to the number of real atoms)
	Possible_detections = ['0','e','g']
	Detection_results = {} # This dictionary allows to store intermediate list of results
	if Nreal == 0:
		return {}
	Detection_results[0] = [[]]
	k = 0
	for i in range(Nreal):
		k += 1
		Detection_results[k] = []
		for result in Detection_results[k-1]:
			for detection in Possible_detections:
				Detection_results[k] += [result + [detection]]
	
	return Detection_results[k]


def ProbasPossibleDetections(Nr,args):
	# This function returns a dictionary with every possible detection associated to their respective probabilities.
	All_detection_results_with_probas = {}
	All_detection_results = PossibleDetections(len(Nr),args)
	for result in All_detection_results:	
		All_detection_results_with_probas[tuple(result)] = Detector_imperfections(list(Nr),result,args)
	
	return All_detection_results_with_probas


def KrausOperator(Ideal_detection,args):
	# Create a Kraus operator corresponding to the list 'Ideal_Detection"
	Tensor_components = []
	for atom in Ideal_detection:
		if atom == 'e':
			Tensor_components += [qt.basis(2,0)*qt.basis(2,0).dag()] # State for an atom ideally detected in e
		if atom == 'g':
			Tensor_components += [qt.basis(2,1)*qt.basis(2,1).dag()] # State for an atom ideally detected in g
	Tensor_components += [qt.qeye(args['NH1']),qt.qeye(args['NH2'])]  # State for both cavities
		
	if False:
		qt.hinton(qt.tensor(Tensor_components)) 
		plt.suptitle("Kraus operator for the detection " + str(Ideal_detection))
		plt.show()
	
	return qt.tensor(Tensor_components)


def Partial_trace(rho,Ideal_detection,args):
	# This function returns the quantity Tr(Sq * rho * Sq.dag()), where Sq is a Kraus operator, and Tr the partial trace on the atoms.
	rho_bis = KrausOperator(Ideal_detection,args) * rho * KrausOperator(Ideal_detection,args).dag()
	#print (rho_bis)
	N_atoms = len(Ideal_detection)
	partial_rho = rho_bis.ptrace((N_atoms,N_atoms+1))
	return partial_rho 
	
	
def SelectDetection(Nreal, Detection, args):
	# This function returns a list with all possible configurations of detection compatible with 'Detection', knowing Nreal
	# Nreal is the number of real atoms, and Detection the list of detected atoms
	Selected_Detection = []
	All_possible_Detections = PossibleDetections(Nreal, args)
	Detection = list(Detection)
	Ng = 0
	Ne = 0
	for i in range(len(Detection)):
		if Detection[i] == 'g':
			Ng += 1
		elif Detection[i] == 'e':
			Ne += 1
		else:
			print ("There must not be anything else than 'e' and 'g' in Detection")
			sys.exit()
			
	# Keep detections with the same number of atoms in g and in e, whatever the order of detection, and whatever the number of non detected atom
	for detection in All_possible_Detections:
		Ng_bis = 0
		Ne_bis = 0
		for i in range(len(detection)):
			if detection[i] == 'g':
				Ng_bis += 1
			elif detection[i] == 'e':
				Ne_bis += 1
		if Ng_bis == Ng and Ne_bis == Ne:
			Selected_Detection += [detection]
						
	return Selected_Detection
				















