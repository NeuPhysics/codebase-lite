#####################################
# This is a test module for
# generating fig on server
#####################################
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import numpy as np
import matplotlib.pylab as plt


xlin_supernova_save = np.load("assets/two-freq-real-ize-xlin-100000000-1.0.npy")
prob0_supernova_save = np.load("assets/two-freq-real-ize-prob0-100000000-1.0.npy")
prob1_supernova_save = np.load("assets/two-freq-real-ize-prob1-100000000-1.0.npy")

range_save = 100000000
stepsize_figure_save = 1.0
rtol_figure_save = 1e-10



fig = plt.figure()
plt.plot(xlin_supernova_save, prob1_supernova_save + prob0_supernova_save,'r-')
plt.title("Total Probability (step="+str(stepsize_figure_save)+", rtol= "+ str(rtol_figure_save) +")",fontsize=20)
plt.xlabel("$\hat x$",fontsize=20)
plt.ylabel("Probability",fontsize=20)
fig.savefig('assets/two-freq-real-ize-total-'+str(range_save)+'-'+str(stepsize_figure_save)+'-'+str(rtol_figure_save)+'.png', dpi=fig.dpi)

