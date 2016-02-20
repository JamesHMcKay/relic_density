#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

m_h=[90,100,110,120,130,140,150,180,200,250,300,400,500,600,800,1000]
gamma=[2.2e-3,2.46e-3,2.82e-3,3.47e-3,4.87e-3,8.12e-3,1.73e-2,6.31e-1,1.43,4.04,8.43,29.2,68,123,304,647]

plt.figure()



A=np.genfromtxt('../Figures/interp_linear.txt',usecols=[0,1])

mh_fit=A[:,0]
gamma_linear=A[:,1]

A=np.genfromtxt('../Figures/interp_cubic.txt',usecols=[0,1])
gamma_cubic=A[:,1]

plt.plot(mh_fit,gamma_cubic,color='black')

plt.plot(mh_fit,gamma_linear,'--',color='blue')

plt.plot(m_h,gamma,'x')


#xlabel(r"Boost velocity (km s$^{-1}$)",fontsize=14)
#ylabel(r"$\chi^2_b/N_b$ in direction of minimum $\chi^2_a$",fontsize=14)
#plt.xlim([90,1000])
#plt.ylim([min(sigma_v(s_plot)),max(sigma_v(s_plot))])

plt.yscale('log')
plt.xscale('log')


plt.savefig("test.eps")