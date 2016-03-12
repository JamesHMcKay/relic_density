#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

plt.figure()



A=np.genfromtxt('../Figures/data/sigma_v.txt',usecols=[0,1,2])
x=A[:,0]
y=A[:,1]
y2=A[:,2]

plt.plot(x,y)
plt.plot(x,y2)


#
#ul=1e6;
#
#fit=interpolate.interp1d(x,y,kind='linear')
#
#
#plt.xlim([min(x),ul])
#plt.ylim([fit(ul)/10,max(y)])
#
#plt.xlabel("$s$",fontsize=18)
#plt.ylabel("$\sigma v_{rel}$",fontsize=18)
#

plt.xscale('log')
plt.yscale('log')


plt.savefig("../Figures/Figures/sigma_v.eps")