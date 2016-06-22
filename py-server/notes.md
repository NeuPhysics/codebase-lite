## List of Files and Outputs

The purpose of this file it to keep an record of what has been done in each file and what to expect in the outputs.

### Two Freq

```
lam0 = 0.845258; # in unit of omegam, omegam = 3.66619*10^-17
dellam = np.array([0.00003588645221954444, 0.06486364865874367]); # deltalambda/omegam
ks = [1.0,1.0/90]; # two k's
thm = 0.16212913985547778; # theta_m
psi0, x0 = [1.0, 0.0, 0.0, 0.0], 0 # initial condition
```

1. `two-freq-real-ize-eq-solver.py` is to solve the system which has been rewritten as real matrix functions.


#### Output

1. `two-freq-real-ize-prob0-1000-1.0.npy` and `two-freq-real-ize-prob1-1000-1.0.npy` are the output as a test of `two-freq-real-ize-eq-solver.py`.
2. `two-freq-real-ize-eq-solver-1e7-step-10-with-tol-func-input-1e-4.out` is the result for range 1e7, stepsize 10, relative tol rtol set using tol function which has input 1e-4 (the total maximum error required)

