#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

plt.figure()



A=np.genfromtxt('../Figures/data/thermal_av.txt',usecols=[0,1,2])
x=A[:,0]
y=A[:,1]
y_interp=A[:,2];
plt.plot(x,y)
plt.plot(x,y_interp,'--')

#xlabel(r"Boost velocity (km s$^{-1}$)",fontsize=14)
#ylabel(r"$\chi^2_b/N_b$ in direction of minimum $\chi^2_a$",fontsize=14)
#plt.xlim([min(x),max(y)])
#plt.ylim([min(sigma_v(s_plot)),max(sigma_v(s_plot))])

#plt.xscale('log')
#plt.yscale('log')


plt.savefig("../Figures/Figures/thermal_av.eps")