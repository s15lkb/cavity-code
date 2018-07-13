Prediction of the state of the cavities


parameters.py
In this file, specify values for different useful parameters related to preparation, interaction between atoms and cavities, and detection.
All parameters are stored in a dictionary, args, which is called by almost all the other files.

state_prediction.py
Use the function RhoFinal(Nrealmax, Detection, args) to plot and return a density matrix which represents the states of the cavities, given :
-	the maximum number of real atoms considered Nrealmax (so the code will calculate density matrices from 0 to Nrealmax atoms, and sum them with the proper weights)
-	the sequence of Detection (post-selection)
The Quantum_map(Nreal, Detection, args)  function is similar to RhoFinal : it just calculates one density matrix for a given real number of atoms (without normalization).

If you have any problem, I suggest you to change debug = False to debug = True, at the beginning of the code.

plotDM.py
This code allows to plot several density matrices corresponding to different values of one varying parameter.
Then, it calculates the best fidelity between the ideal state |10 > + exp(i*phi) |01 > and the simulated density matrix (it finds the value of phi which maximises the fidelity), for each value of the varying parameter.
Finally, it shows the fidelity as a function of the varying parameter. 

At the beginning of the code, specify the number of plots you want, the varying parameter, its range of values, and the maximum number of real atoms considered. Then, run the code. 

preparation.py
This code calculates the probabilities to prepare a bunch with Ng atoms in the fundamental state and Ne atoms in the excited state, as a function of the mean number of atoms in each bunch, mu, and the preparation errors : eta_prep_g and eta_prep_e
Don’t forget « from parameters.py import args » in the Console.

cavities_interaction.py
This code calculates the density matrix after interaction with both cavities. It plots the number of photons in each cavity as well as the probabilities for the atoms to be in the excited state, as functions of time.
Don’t forget « from parameters.py import args » in the Console.

detection.py
This codes models detector imperfections by taking into account :
- the detector efficiency (corresponding to the number of detected atoms vs the number of real atoms) : epsilon
- the detection errors (on atom states) : eta_e & eta_g
The function Detector_imperfections is based on the table p.58 in Xingxing ZHOU's thesis.

This code also calculates possible detections (with their probabilities) given a specific number of real atoms, and several functions useful for the computation of quantum_maps, done in the file state_prediction.py .
Don’t forget « from parameters.py import args » in the Console.

Other files
The other files are called by cavities_interaction.py : they all deal with interaction between atoms and cavities.




