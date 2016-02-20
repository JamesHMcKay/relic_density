#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

plt.figure()



A=np.genfromtxt('../Figures/integrand.txt',usecols=[0,1])

s=A[:,0]
integrand=A[:,1]

plt.plot(s,integrand,color='black')


#xlabel(r"Boost velocity (km s$^{-1}$)",fontsize=14)
#ylabel(r"$\chi^2_b/N_b$ in direction of minimum $\chi^2_a$",fontsize=14)
plt.xlim([min(s),max(s)])
#plt.ylim([min(sigma_v(s_plot)),max(sigma_v(s_plot))])

plt.xscale('log')
plt.yscale('log')


plt.savefig("../Figures/Figures/integrand.eps")