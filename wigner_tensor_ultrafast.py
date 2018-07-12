# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 09:05:18 2017

@author: valentin.metillon
"""

from numpy import meshgrid, sqrt, array, linspace, conj, exp, tensordot, real, zeros, shape, pi, cos, sin, diag, append, loadtxt
import matplotlib.pyplot as plt 
import matplotlib as mpl
from matplotlib import cm
from scipy import sparse
from scipy.sparse import lil_matrix


    








def wigner_iterative_basis(NH,xvec,yvec, square_grid = True):
    #Computes iteratively the values Wmn(alpha) for |m><n| on the chosen grid (alpha = xvec+1j*yvec)
    #W is normalized to be equal to the parity of the photon number
    if square_grid == True:
        X, Y = meshgrid(xvec, yvec)
        A = (X + 1.0j * Y)
    
    Wlist = array([zeros(shape(A), dtype=complex) for k in range(NH**2)])
    Wlist[0] = exp(-2.0 * abs(A) ** 2)
    
        
    for n in range(1, NH):
        Wlist[n] = (2.0 * A * Wlist[n - 1]) / sqrt(n)
    
    for m in range(1, NH):
        Wlist[m*NH+m] = (2 * conj(A) * Wlist[(m-1)*NH+m]  - sqrt(m) * Wlist[(m-1)*NH+m-1]) / sqrt(m)
        
        for n in range(m+1 , NH):
            Wlist[m*NH+n] = (2 * conj(A) * Wlist[(m-1)*NH+n]  - sqrt(n) * Wlist[(m-1)*NH+n-1]) / sqrt(m)
            
    return Wlist



def wigner_2m(rho, NHa, NHb, xa, ya, xb, yb, rho_min = 10**-6):
    #Computes the 2-mode Wigner function of rho using the functions for |m><n| in the Fock basis 
    #computed by wigner_iterative_basis.
    #rho_min allows to remove the small terms of rho to reduce computation time.
    Nxa = len(xa)    
    Nya = len(ya)
    Nxb = len(xb)
    Nyb = len(yb)
    
    
    rho.data= lil_matrix(rho.data)
    for i in range(rho.shape[0]):
        for j in range(rho.shape[1]):
            if abs(rho[i,j])<rho_min:
                rho.data[i,j]= 0      #Only the terms bigger than 10^-6 in rho are kept for the computation
    
    
    a,b=sparse.find(rho.data)[0],sparse.find(rho.data)[1]   #indices for which the terms of the truncated matrix are non-zero
    
    W = tensordot(zeros((Nya,Nxa)),zeros((Nyb,Nxb)),0)
    
    Wa =  wigner_iterative_basis(NHa, xa, ya)    #Builds the wigner functions for the |i><j| basis for Ha, on the chosen values of x and y
    Wb =  wigner_iterative_basis(NHb, xb, yb)    #Builds the wigner functions for the |k><l| basis for Hb, on the chosen values of x and y
    for p in range(len(a)):
        i = a[p]//NHb  #The terms of the matrix have the form |i><j|X|k><l|
        j = b[p]//NHb    
        k = a[p]% NHb
        l = b[p]% NHb   
        
        if i==j and k==l:
            W+= real(rho[NHb*i+k,NHb*i+k]*tensordot(Wa[i*NHa+i],Wb[k*NHb+k],0))   #diagonal terms
        else:
            if j>i:
                if l>k:    #Non diagonal terms are computed as 2*real(rho_kl W^|k><l|), using rho_kl W^|k><l| = conj (rho_lk W^|l><k|)
                    W+= 2*real(rho[NHb*i+k,NHb*j+l]*tensordot(Wa[i*NHa+j],Wb[k*NHb+l],0))
                else:
                    W+= 2*real(rho[NHb*i+k,NHb*j+l]*tensordot(Wa[i*NHa+j],Wb[l*NHb+k].conjugate(),0)) 
                    #wigner_iterative_basis doesn't compute the W^|k><l| if k>l, in order to spare time
            elif j==i:  #if i>j, the terms are already taken into account by the trick above
                if l>k: #if k>l, the terms are already taken into account by the trick above
                    W+= 2*real(rho[NHb*i+k,NHb*j+l]*tensordot(Wa[i*NHa+j],Wb[k*NHb+l],0))
    return(W)
            



def plot_wigner(W, xa, ya, xb, yb):
    #Create the subplots for the 2-mode wigner function W calculated on grid (xa, ya, xb, yb).
    #To plot a subset of the function, cut the lists and arrays first.
    Nxa = len(xa)
    Nya = len(ya)
    Nxb = len(xb)
    Nyb = len(yb)
    fig, axes = plt.subplots(Nya, Nxa, figsize=(3, 3))
    
    nrm = mpl.colors.Normalize(-1, 1)
    
    K = zeros((Nyb,Nxb))
    K[0,0]=1    
    K[0,1]=-1
    
    plt1 = axes[0,0].contourf(xb, yb, K[:,:], 100,cmap=cm.RdBu,norm=nrm)	
    cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
    plt.colorbar(plt1, cax=cax, norm=nrm)
    
    for k in range(Nya):
        for l in range(Nxa):
            #axes[k,l].pcolor(xb,xb,W[k,l,:,:])
            plt1 = axes[k,l].contourf(xb, yb, W[Nya-1-k,l,:,:], 100,cmap=cm.RdBu,norm=nrm)
            axes[k,l].get_xaxis().set_visible(False)
            axes[k,l].get_yaxis().set_visible(False)








#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
### for Palatino and other serif fonts use:
##rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
#
#plt.rcParams["text.latex.preamble"] = [r"\usepackage{amsfonts}"]
#
#
#NHa = 10
#NHb = 10
#
#Nxa = 11
#Nxb = 51
#Nya = 11
#Nyb = 51
#xa = linspace(-1.5,1.5,Nxa)
#xb = linspace(-1.5,1.5,Nxb)
#ya = linspace(-1.5,1.5,Nya)
#yb = linspace(-1.5,1.5,Nyb)
#

#xa = array([sqrt(3/8)])
#ya = array([0])
#
#theta = linspace(0,2*pi,201)
#xb = cos(theta)
#yb = sin(theta)




#rho_min = 10**-6  #Si |rho_ijkl|<rho_min, le terme correspondant dans la Wigner n'est pas calculé.
#
#
#
#
#
#na = 0
#nb = 1
#alpha = 2
#
#
#phi = pi/2
#def mnnm_state(m,n,phi,NHa,NHb):
#    return (qt.tensor(qt.fock(NHa,m),qt.fock(NHb,n))+exp(1j*phi)*qt.tensor(qt.fock(NHa,n),qt.fock(NHb,m)))/sqrt(2)
#
##psi_p = (qt.tensor(qt.fock(NHa,na),qt.fock(NHb,nb))+qt.tensor(qt.fock(NHa,nb),qt.fock(NHb,na)))/sqrt(2)
##psi_m = (qt.tensor(qt.fock(NHa,na),qt.fock(NHb,nb))-qt.tensor(qt.fock(NHa,nb),qt.fock(NHb,na)))/sqrt(2)
##psi_phi = (qt.tensor(qt.fock(NHa,na),qt.fock(NHb,nb))+exp(1j*phi)*qt.tensor(qt.fock(NHa,nb),qt.fock(NHb,na)))/sqrt(2)
#
#psi_p = mnnm_state(1,0,0,NHa,NHb)
#psi_m = mnnm_state(1,0,pi,NHa,NHb)
#psi_phi = mnnm_state(1,0,phi,NHa,NHb)
#
#a=loadtxt('rhoML_real.dat')
#b=loadtxt('rhoML_imag.dat')
#rho_igor = qt.Qobj(a+1j*b)
##psi = (qt.tensor(qt.coherent(NHa,alpha),qt.coherent(NHb,alpha))+qt.tensor(qt.coherent(NHa,-alpha),qt.coherent(NHb,-alpha)))
##psi = psi/psi.norm()
#rho_p = qt.ket2dm(psi_p)
#rho_m = qt.ket2dm(psi_m)
#rho_phi = qt.ket2dm(psi_phi)
#rho_coh = .5*(rho_p-rho_m)
#rho_0 = qt.ket2dm(qt.tensor(qt.fock(NHa,0),qt.fock(NHb,0)))
#rho_1 = qt.ket2dm(qt.tensor(qt.fock(NHa,0),qt.fock(NHb,1)))
#rho_mix = .5*(rho_p+rho_m)
###State with Schmidt rank 3
#psi = (qt.tensor(qt.fock(NHa,0),qt.fock(NHb,0))+qt.tensor(qt.fock(NHa,1),qt.fock(NHb,1))+\
#qt.tensor(qt.fock(NHa,2),qt.fock(NHb,2)))/sqrt(3)
#rho = qt.ket2dm(psi)

###Separable state
#rho = qt.ket2dm(qt.tensor(qt.fock(NHa,1),qt.fock(NHb,1)))


###Plot full wigner

#W = wigner_2m(rho_p, NHa, NHb, xa, ya, xb, yb)
#W = wigner_2m(rho_m, NHa, NHb, xa, ya, xb, yb)
#W = wigner_2m(rho_phi, NHa, NHb, xa, ya, xb, yb)
#W = wigner_2m(rho_coh, NHa, NHb, xa, ya, xb, yb)
#plot_wigner(W,xa,ya,xb,yb)





###Fringes for several values

#fig , ax = plt.subplots()
#rs = linspace(0,1,11)
#for i,r in enumerate(rs): 
#    xa = array([r])
#    W = wigner_2m(rho_p, NHa, NHb, xa, ya, r*xb, r*yb)
#    ax.plot(theta,diag(W[0,0,:,:]), label = "r = "+str(r))
#ax.legend(loc='lower right')


###Contrast w.r.t. radius
#rs = linspace(0,2,101)
#C = zeros(len(rs))
#for i,r in enumerate(rs): 
#    xa = array([r])
#    W = wigner_2m(rho_p, NHa, NHb, xa, ya, r*xb, r*yb)
#    C[i]= diag(W[0,0,:,:]).max()-diag(W[0,0,:,:]).min()
#fig , ax = plt.subplots()
#ax.plot(rs,C, label= r'Contrast as a function of $|\alpha| = |\beta|$')
#ax.legend()

###radial Wigner

#Nxa = 101
#Nxb = 101
#Nya = 1
#Nyb = 1
#xa = linspace(-1.5,1.5,Nxa)
#xb = linspace(-1.5,1.5,Nxb)
#ya = array([0])
#yb = array([0])
#
#W_igor = wigner_2m(rho_igor, 5, 5,xa,ya,xb,yb)
#W_p = wigner_2m(rho_p, NHa, NHb,xa,ya,xb,yb)
#W_m = wigner_2m(rho_m, NHa, NHb,xa,ya,xb,yb)
#W_phi = wigner_2m(rho_phi, NHa, NHb,xa,ya,xb,yb)
#W_0 = wigner_2m(rho_0, NHa, NHb,xa,ya,xb,yb)
#W_1 = wigner_2m(rho_1, NHa, NHb,xa,ya,xb,yb)
#W_mix = wigner_2m(rho_mix, NHa, NHb,xa,ya,xb,yb)
#t_end =time()
#
#fig , ax = plt.subplots()
#ax.plot(xb,diag(W_igor[0,:,0,:]),'m',label=r'Etat reconstruit')
#ax.plot(xb,diag(W_p[0,:,0,:]),'b',label=r'$\left|10\right\rangle+\left|01\right\rangle$')
#ax.plot(xb,diag(W_m[0,:,0,:]),'r',label=r'$\left|10\right\rangle-\left|01\right\rangle$')
#ax.plot(xb,diag(W_0[0,:,0,:]),'c',label=r'$\left|00\right\rangle$')
#ax.plot(xb,(4*xb**2-1)*exp(-2*xb**2),'black', label = r'wigner monomode $\left|1\right\rangle$')
#ax.plot(xb,diag(W_phi[0,:,0,:]),'g',label=r'$\left|10\right\rangle+e^{i\phi} \left|01\right\rangle$')
#ax.plot(xb,diag(W_1[0,:,0,:]),'g',label=r'$\left|01\right\rangle$')
#ax.plot(xb,diag(W_mix[0,:,0,:]),'g',label=r'mixture')
#x_m = array([.5])
#ax.plot([x_m,x_m],[wigner_2m(rho_m, NHa, NHb,x_m,ya,x_m,yb)[0,0,0,0],wigner_2m(rho_p, NHa, NHb,x_m,ya,x_m,yb)[0,0,0,0]],'b-')
#ax.set_title(r'Fonctions de Wigner pour $\alpha = \beta \in \mathbb{R}$')
#plt.legend()	
#
#plt.show()