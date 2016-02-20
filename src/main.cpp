#include "cross_section.hpp"
#include "interp.hpp"
#include "figures.hpp"
#include "integrate.hpp"

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

using namespace std;

int main()
{
double m_s=100;
Cross_section cross_section(m_s);
Relic_density relic_density(m_s);


simple func(10);
int m=30;
//cout << "func = " << func(100) << endl;
double val;

//Trapzd<simple> s(func,0,30);
//for(int j=1;j<=m+1;j++) val=s.next();
//cout << "integral = " << val << endl;
// test interpolator


//

double s=4*pow(m_s,2), T_=0.0001;


cs_func func2(T_,m_s);

cout << "integral = " << qtrap(func2, s,40001.2,1e-5) << endl;

cout<< "func = " << func2(s) << " " << func2(40001.2) << endl;

//cout << "x_f = " << relic_density.x_f() << endl;

//

//cout << "true bessel = " << boost::math::cyl_bessel_k(2,100) << endl;
//
//cout << "asymptotic approximation bessel = " << cross_section.bessel_k_asymtotic(2,100) << endl;
//
//cout << "exp of log_asymptotic approximation bessel = " << exp(cross_section.log_bessel_k_asymtotic(2,100)) << endl;
//



//double x_f=23;//cross_section.x_f();

//cout << "CALCULATING A_F " << endl;
//cout << "A(x_f) = " << cross_section.A(x_f) << endl;



//cout << "Generating figures" << endl;

Figures figures(cross_section, relic_density);
//figures.plot_Z();

//cout<< "thermal average is = " << cross_section.Z(10000000000) << endl;

figures.plot_sigma_v(0.0001);
//figures.plot_Z();

//double Y=cross_section.Y_today(23);//cross_section.x_f());
//double H_0=70;

//double rho_crit=1.05375e-5;// 3*(pow(H_0,2))/(8*Pi*G);
//
//double T=2.73*8.62e-14;
//double g=100;
//double Pi=3.14;
//double s_0=2890;//(2*pow(Pi,2)/45)*g*pow(T,3);


//cout << "mass fraction = " <<  Y*m_s*s_0/rho_crit << endl;

}