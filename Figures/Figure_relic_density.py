#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy import interpolate
from scipy import integrate
from scipy import special

plt.figure()



A=np.genfromtxt('../Figures/data/relic_density.txt',usecols=[0,1,2])
lambda_hs=log10(A[:,0])
ms=A[:,1]
rd=log10(A[:,2]/0.11)

xi = linspace(min(ms), max(ms))
yi = linspace(min(lambda_hs), max(lambda_hs))
zi = griddata((ms, lambda_hs), rd, (xi[None,:], yi[:,None]),method='cubic')

print zi

levels=[-3,-1,0,5,5,6,7,8,9,10]
CS=plt.contour(xi, yi, zi,levels=levels)

levels2=[0,10000]

plt.contourf(xi,yi,zi,levels=levels2)



plt.clabel(CS, levels,inline=1,fmt='%1.1f',fontsize=14)


plt.xlabel("$M_S$ (GeV)",fontsize=18)
plt.ylabel("$\log_{10}\lambda_{hs}$",fontsize=18)


plt.ylim([min(lambda_hs),max(lambda_hs)])
plt.yticks([0,-1,-2,-3])
#plt.yscale('log')
#plt.colorbar()


plt.savefig("../Figures/Figures/relic_density.eps")





