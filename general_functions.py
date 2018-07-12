""" 
GENERAL FUNCTIONS

Author : Gautier Depambour

Some useful general functions.

"""
import numpy as np

def Factoriel(n):
	factoriel = 1
	i = n
	while i > 1:
		factoriel = factoriel * i
		i = i - 1
	return factoriel
	
def Poisson(mu, n):
	return np.exp(-mu) * (mu**n) / Factoriel(n)
	
def DistributionPoisson(mu,n): #for n in 0,1,2
	if n>2:
		print ("n doit etre inferieur ou egal a 2 !")
	return Poisson(mu,n) / (Poisson(mu,0) + Poisson(mu,1) +Poisson(mu,2) )
	
def ListePoisson(mumin, mumax, pas, n):
	ListeMu = np.linspace(mumin, mumax, pas)
	LPoisson = []
	for mu in ListeMu:
		LPoisson += [Poisson(mu,n)]
	return LPoisson

def Binome(n):
	if n == 0:
		return [1]
	Coeffs =[1,1]
	for i in range(n-1):
		Coeffsmilieu = []
		for j in range(1,len(Coeffs)):
			Coeffsmilieu += [Coeffs[j-1] + Coeffs[j]]
		Coeffs = [1] + Coeffsmilieu + [1]
	return Coeffs
	
def CoeffBinome(n,k):
	ListeCoeffs = Binome(n)
	return ListeCoeffs[k]
