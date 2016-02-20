#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

plt.figure()



A=np.genfromtxt('../Figures/data/sigma_v.txt',usecols=[0,1])
x=A[:,0]
y=A[:,1]

plt.plot(x,y)


#xlabel(r"Boost velocity (km s$^{-1}$)",fontsize=14)
#ylabel(r"$\chi^2_b/N_b$ in direction of minimum $\chi^2_a$",fontsize=14)
plt.xlim([min(x),1e5])
#plt.ylim([1e-20,1e-9])

plt.xscale('log')
plt.yscale('log')


plt.savefig("../Figures/Figures/sigma_v.eps")