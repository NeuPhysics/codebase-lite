import numpy as np
from scipy.integrate import odeint
from scipy.integrate import ode
import matplotlib.pylab as plt


# Parameters are shared by ALL methods in this notebook

endpoint = 1000000; # integration range
dx = 10.0; # step size
lam0 = 0.845258; # in unit of omegam, omegam = 3.66619*10^-17
dellam = np.array([0.00003588645221954444, 0.06486364865874367]); # deltalambda/omegam
ks = [1.0,1.0/90]; # two k's
thm = 0.16212913985547778; # theta_m

savestep = 1;

### Real System
psi40, x40 = [1.0, 0.0, 0.0, 0.0], 0 # initial condition

xlin4 = np.arange(dx,endpoint+1*dx, dx)

psi4 = np.zeros([len(xlin4)  , 4])


xlin4save = np.zeros(len(xlin4)/savestep);
psi4save = np.zeros([len(xlin4save)  , 5])



#########################
# Make the equation all Real
#########################

def hamiltonian(x, deltalambda, k, thetam):

    return np.array( [[ 0,   0.5* np.sin(2*thetam) * ( deltalambda[0] * np.sin(k[0]*x) + deltalambda[1] * np.sin(k[1]*x) ) * np.exp( 1.0j * ( - x - np.cos(2*thetam) * (  ( deltalambda[0]/k[0] * np.cos(k[0]*x) + deltalambda[1]/k[1] * np.cos(k[1]*x) ) )  ) )     ],   [ 0.5* np.sin(2*thetam) * ( deltalambda[0] * np.sin(k[0]*x) + deltalambda[1] * np.sin(k[1]*x) ) * np.exp( -1.0j * ( - x - np.cos(2*thetam) * ( deltalambda[0] /k[0] * np.cos(k[0]*x) + deltalambda[1] /k[1] * np.cos(k[1]*x) )  ) ), 0 ]]  ) # Hamiltonian for double frequency


def hamiltonian4(x, deltalambda, k, thetam):

    hr = np.array(hamiltonian(x, deltalambda, k, thetam)).real;
    hi = np.array(hamiltonian(x, deltalambda, k, thetam)).imag;

    return np.array([[hi[0][0],hi[0][1],hr[0][0],hr[0][1]], [hi[1][0],hi[1][1],hr[1][0],hr[1][1]], [- hr[0][0], - hr[0][1], hi[0][0], hi[0][1]],  [- hr[1][0], - hr[1][1], hi[1][0], hi[1][1]] ] )

def sysdpsidt(x, psi, deltalambda, k, thetam):

    return np.dot(hamiltonian4(x, deltalambda, k, thetam), [psi[0], psi[1], psi[2], psi[3]])


atol_req = 1e-8

sol4 = ode(sysdpsidt).set_integrator('dopri5', atol=1e-8)
sol4.set_initial_value(psi40, x40).set_f_params(dellam,ks,thm)

flag4 = 0
flag4save = 0

while sol4.successful() and sol4.t < endpoint:
    sol4.integrate(xlin4[flag4])
    if np.mod(flag4,savestep)==0:
        psi4save[flag4save] = [sol4.t, sol4.y[0],sol4.y[1],sol4.y[2],sol4.y[3]]

        with open(r"assets/ode-dopri5-range-"+str(endpoint)+"-step-"+str(dx)+"-atol-"+str(atol_req)+".csv", 'a') as f_handle:
            np.savetxt(f_handle, psi4save[flag4save])

        flag4save = flag4save + 1
    flag4 = flag4 + 1
    #   print sol.t, sol.y





print "successful"
print "calculation range is "+str(endpoint)
print "step size is"+str(dx)
print "atol_req is"+str(atol_req)
