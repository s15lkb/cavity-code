# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 17:25:56 2018

@author: valentin.metillon
"""

import qutip as qt
from numpy import sqrt

def contrast(rho):
    #Ramsey contrast for the final state of an atom that was prepared in e+g
    x_m = (rho*qt.sigmax()).tr().real
    y_m = (rho*qt.sigmay()).tr().real
    return(sqrt(x_m**2+y_m**2))

    