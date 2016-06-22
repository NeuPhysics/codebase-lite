import numpy as np
from scipy.integrate import odeint
from scipy.integrate import ode
import matplotlib.pylab as plt
import csv
import time



endpoint = 10000000; # integration range
dx = 10.0; # step size
lam0 = 0.845258; # in unit of omegam, omegam = 3.66619*10^-17
dellam = np.array([0.00003588645221954444, 0.06486364865874367]); # deltalambda/omegam
ks = [1.0,1.0/90]; # two k's
thm = 0.16212913985547778; # theta_m

psi0, x0 = [1.0, 0.0, 0.0, 0.0], 0 # initial condition


def hamiltonian(x, deltalambda, k, thetam):

    return [[ 0,   0.5* np.sin(2*thetam) * ( deltalambda[0] * np.sin(k[0]*x) + deltalambda[1] * np.sin(k[1]*x) ) * np.exp( 1.0j * ( - x - np.cos(2*thetam) * (  ( deltalambda[0]/k[0] * np.cos(k[0]*x) + deltalambda[1]/k[1] * np.cos(k[1]*x) ) )  ) )     ],   [ 0.5* np.sin(2*thetam) * ( deltalambda[0] * np.sin(k[0]*x) + deltalambda[1] * np.sin(k[1]*x) ) * np.exp( -1.0j * ( - x - np.cos(2*thetam) * ( deltalambda[0] /k[0] * np.cos(k[0]*x) + deltalambda[1] /k[1] * np.cos(k[1]*x) )  ) ), 0 ]]   # Hamiltonian for double frequency

def hamiltonian4(x, deltalambda, k, thetam):
    hr = np.array(hamiltonian(x, deltalambda, k, thetam)).real;
    hi = np.array(hamiltonian(x, deltalambda, k, thetam)).imag;


    return [[hi[0][0],hi[0][1],hr[0][0],hr[0][1]], [hi[1][0],hi[1][1],hr[1][0],hr[1][1]], [- hr[0][0], - hr[0][1], hi[0][0], hi[0][1]],  [- hr[1][0], - hr[1][1], hi[1][0], hi[1][1]] ]

def sysdpsidt(psi, x, deltalambda, k, thetam):

    return np.dot(hamiltonian4(x, deltalambda, k, thetam), [psi[0], psi[1], psi[2], psi[3]])

def sysjac(psi, x, deltalambda, k, thetam):

    return hamiltonian4(x, deltalambda, k, thetam)

def integral_tol(total_error_needed,totalrange, stepsize): # tolenrance of the integral that we require
    return total_error_needed*stepsize/totalrange

rtol_req = integral_tol(1e-4,endpoint,dx)
xlin = np.linspace(0, endpoint, np.floor(endpoint/dx) )
solodeint = odeint(sysdpsidt, psi0, xlin, args = (dellam,ks,thm), full_output = 1,rtol=rtol_req)

prob0=solodeint[0][:,0]**2+solodeint[0][:,2]**2
prob1=solodeint[0][:,1]**2+solodeint[0][:,3]**2


np.save("assets/two-freq-real-ize-xlin-"+str(endpoint)+"-"+str(dx),xlin)
np.save("assets/two-freq-real-ize-prob0-"+str(endpoint)+"-"+str(dx),prob0)
np.save("assets/two-freq-real-ize-prob1-"+str(endpoint)+"-"+str(dx),prob1)

print solodeint[1]['message']
print solodeint[0][-1], xlin[-1]
print "step size is"+str(dx)
print "rtol_req is"+str(rtol_req)
