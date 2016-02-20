#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special
#def lamcutoff(LB):
#    return -0.065/(1-0.01*log(246.221/LB))

# can assume that h_eff=g_eff

h_eff=100
g_eff=100
dh_eff_dt=0
Pi=3.14159265359
M_pl=1e19

lambda_hs=0.1
v_0=246
m_s=100
Mh=125
m_h=[90,100,110,120,130,140,150,180,200,250,300,400,500,600,800,1000]
gamma=[2.2e-3,2.46e-3,2.82e-3,3.47e-3,4.87e-3,8.12e-3,1.73e-2,6.31e-1,1.43,4.04,8.43,29.2,68,123,304,647]

Gamma_h_func=interpolate.interp1d(m_h,gamma,kind='linear')

def Gamma_h(input):
  if (input>min(m_h)) & ( input<max(m_h)):
    return Gamma_h_func(input)
  else:
    return 0;



def Dh2(s):
    return 1/((s-Mh**2)**2+(Mh**2)*Gamma_h(Mh)**2)
def sigma_v(s):
    return (2*(lambda_hs**2)*(v_0**2)/(sqrt(s)))*(abs(Dh2(s)))*Gamma_h(sqrt(s))



#for i in range(1,10):
#  y=20.5+log((10^26)*cross_section)+log(m)-0.5*log(g_star)

# calculate cross-section

def integrand(s,m_s,T):
    return s*sqrt(s-4*m_s**2)*special.kn(1,sqrt(s)/T)*sigma_v(s)/ (16 * T * (m_s**4) * (special.kn(4,m_s/T)**2))

def x(T):
    return m_s/T;

def cross_section(T):
    return integrate.quad(integrand ,4*m_s**2,1000000,args=(m_s,T))

def g(T):
    return 100 #(h_eff/sqrt(g_eff)) * (1 + (T/(3*h_eff)) * dh_eff_dt)


def Z(x):
    return sqrt(Pi/45)*((m_s*M_pl)/(x**2)) * (sqrt(100))*(cross_section(m_s/x)[0])

def Y_eq(x):
    return special.kn(2,x)*(45/(4*Pi**4))*(x**2)/(h_eff)

def dY_eq_hat(x):
    return (special.kvp(2,x,1)*(2*45/(4*Pi**4))*(x**2)/(h_eff) + 2*Y_eq(x)/x+Y_eq(x)) * exp(x)

def hat(x,F):
    return F*exp(x)


x_star=0.1
deltaf=0.618


for i in range(1,10):
  x_star=log( (deltaf*(2+deltaf)/(1+deltaf)) * (Z(x_star)*(hat(x_star,Y_eq(x_star))**2))/ (hat(x_star,Y_eq(x_star))-dY_eq_hat(x_star)) )
  print x_star
x_sol=4
for i in range(1,10):
 x_sol=20.5+log((10**26)*2.2e-26)+log(m_s)-0.5*log(g_eff)+0.5*log(x_sol)-log(x_sol-1.5)
 print x_sol

#
#print cross_section(100)
#
#plt.figure()
#
#s_plot=np.linspace(4*m_s**2,1000000,100)
##plt.vlines(Mh**2,min(sigma_v(s_plot)),max(sigma_v(s_plot)))
#T_plot=np.linspace(90,900,100)
#
#cs=zeros(size(T_plot))
#for i in range(1,100):
#  cs[i]=cross_section(T_plot[i])[0]
#print 'cross section', cross_section(100)[0]
#
#plt.plot(T_plot,cs,color='black')
##plt.plot(v,2.3*ones(size(v)),'--',color='black')
#
##xlabel(r"Boost velocity (km s$^{-1}$)",fontsize=14)
##ylabel(r"$\chi^2_b/N_b$ in direction of minimum $\chi^2_a$",fontsize=14)
##plt.xlim([90,1000])
##plt.ylim([min(sigma_v(s_plot)),max(sigma_v(s_plot))])
#
#plt.yscale('log')
##plt.xscale('log')
#
#
#plt.savefig("test.eps")


# integration for final result


Y_f=(1+deltaf)*Y_eq(x_star)

def A_f(input):
    return integrate.quad(Z,input,input*100)

Y_today=Y_f/(1+Y_f*A_f(x_star)[0])

print Y_today


#for i in range(1,10000):
#  print m_s/float(i), cross_section(m_s/float(i))[0]






