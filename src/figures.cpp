#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <fstream>

#include "cross_section.hpp"
#include "interp.hpp"
#include "figures.hpp"

using namespace std;



void Figures::plot_thermal_av()
{
//Cross_section cross_section(m_s);
Relic_density relic_density(m_s);
//plot integrand as a function of s
ofstream myfile_integrand;
myfile_integrand.open ("../Figures/data/thermal_av.txt");
//double m_t=173;
//double lower=4*pow(m_s,2);
//double upper=1e15;//pow(1000,2)*1e40;

double lower = 0.5;
double upper = 2;


relic_density.thermal_average_make_interp(lower*0.9,upper*1.1,5);

Relic_density relic_density2(m_s);

cout << "test value = " << relic_density.calc_cross_section(lower) << endl;

double T=0;
for (int i=0;i<10;i++)
{
T=(float(i)/10)*(upper-lower)+lower;
myfile_integrand << T << " " << relic_density.calc_cross_section(T) << " " << relic_density2.calc_cross_section(T) << endl;
}
myfile_integrand.close();

system("python ../Figures/Figure_thermal_av.py");

}

void Figures::plot_sigma_v(double T)
{
//Cross_section cross_section(m_s);
//plot integrand as a function of s
ofstream myfile_integrand;
myfile_integrand.open ("../Figures/data/sigma_v.txt");
//double m_t=173;
double lower=4*pow(m_s,2);
double upper=1e15;//pow(1000,2)*1e40;
double s,p;

cs_func func(T,m_s);

for (int i=0;i<400;i++)
{
p=(log10(upper)-log10(lower))/400*float(i);
s=pow(10,p)+lower;
//myfile_integrand << s << " " << cross_section.cs_integral(s,T)<<" "<< func(s) << endl;

myfile_integrand << s << " " << cross_section.sigma_v(s) <<" "<< cross_section.sigma_v(s) << endl;
}
myfile_integrand.close();

system("python ../Figures/Figure_sigma_v.py");

}


void Figures::plot_Z()
{
ofstream myfile_integrand;
myfile_integrand.open ("../Figures/data/Z.txt");
Z_func func_z(m_s);

double x_upper=1e3, x_lower=20;
double T_lower=m_s/x_upper,T_upper=m_s/x_lower;


Relic_density rd(m_s);


rd.thermal_average_make_interp(T_lower,T_upper,4);

Z_func func_z2(m_s, rd.T_m, rd.thermal_av_m);



double x;
//double m_t=173;
for (int i=0;i<100;i++)
{
x=(x_upper-x_lower)*(float(i)/float(100))+x_lower;

myfile_integrand << x<< " " << func_z(x) << " " << func_z2(x) << endl;
//cout << pow(10,x) << " " <<func_z(pow(10,x)) << endl;
}
myfile_integrand.close();
system("python ../Figures/Figure_Z.py");
}

void Figures::gamma_h()
{
std::vector<double> x  ={90,100,110,120,130,140,150,180,200,250,300,400,500,600,800,1000};
std::vector<double> y  ={2.2e-3,2.46e-3,2.82e-3,3.47e-3,4.87e-3,8.12e-3,1.73e-2,6.31e-1,1.43,4.04,8.43,29.2,68,123,304,647};

Linear_interp myfunc(x,y);
Poly_interp myfunc_poly(x,y,4);
double m_h=444,gamma,gamma_2;

gamma=myfunc.interp(m_h);
gamma_2=myfunc_poly.interp(m_h);
cout<< "linear: " <<gamma << endl;
cout<< "polynomial: " <<gamma_2 << endl;


ofstream myfile_linear;
ofstream myfile_cubic;
myfile_linear.open ("../Figures/data/interp_linear.txt");
myfile_cubic.open ("../Figures/data/interp_cubic.txt");


for (int i=1;i<400;i++)
{
m_h=float(i)*(990-100)/400+90;

myfile_linear << m_h << " " << myfunc.interp(m_h) << endl;
myfile_cubic  << m_h << " " << myfunc_poly.interp(m_h) << endl;

}

myfile_linear.close();
myfile_cubic.close();

}

