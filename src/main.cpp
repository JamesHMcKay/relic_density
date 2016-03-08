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


double mass_frac(double lambda_hs, double m_s)
{
Relic_density relic_density(m_s,lambda_hs);
double Y=relic_density.Y_today(relic_density.x_f());//cross_section.x_f());
double rho_crit=1.05375e-5;// 3*(pow(H_0,2))/(8*Pi*G);  units: h^2 GeV /( c^2 cm^3)
double s_0=2890;//(2*pow(Pi,2)/45)*g*pow(T,3);
return  Y*m_s*s_0/rho_crit;
}




void plot_rd()
{
int x_pts=20,y_pts=20;
double m_lower=50, m_upper=70;
double lam_lower_log=-4, lam_upper_log=0;
double lambda_hs,ms;

ofstream myfile;
myfile.open ("../Figures/data/relic_density.txt");


double lam_step= (lam_upper_log-lam_lower_log)/y_pts;

double m_step= (m_upper-m_lower)/x_pts;

for (int i=0; i<20 ;i++)
{
 lambda_hs=pow(10,-lam_step*float(i));
 
  for (int j=0; j<20 ; j++)
  {
   ms=m_step*float(j)+m_lower;
   
   cout << "calculating point " << lambda_hs << " " << ms << endl;
   
   
   myfile << lambda_hs << " " << ms << " " << mass_frac(lambda_hs,ms) << endl;
   
  }
}
myfile.close();


}













int main()
{
double m_s=2000,lambda_hs=0.3;
Cross_section cross_section(m_s,lambda_hs);
Relic_density relic_density(m_s,lambda_hs);


simple func(10);
int m=30;
//cout << "func = " << func(100) << endl;
double val;

plot_rd();


//Trapzd<simple> s(func,0,30);
//for(int j=1;j<=m+1;j++) val=s.next();
//cout << "integral = " << val << endl;
// test interpolator


//

double s=4*pow(m_s,2), T_=0.0001;


//cs_func func2(T_,m_s);

//cout << "integral = " << qtrap(func2, s,40001.2,1e-5) << endl;





//cout << "x_f = " << relic_density.x_f() << endl;


//Z_func func_z(m_s);

//double t;

//cin >> t ;

//cout << "cs integral = " << relic_density.calc_cross_section(t) << " at " << t <<  endl;
//cs_func func3(t,m_s);
//cout << "cs func at s = " << s+1 << " is " << func3(s+1) << endl;
//cout << "cs func at 1.001 s = " << func3(1.001*s) << endl;

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

//Figures figures(cross_section, relic_density);
//figures.plot_Z();

//figures.plot_thermal_av();

//cout<< "thermal average is = " << cross_section.Z(10000000000) << endl;

//figures.plot_sigma_v(0.00045);

//cout<< "integral is (thermal av) = " << relic_density.calc_cross_section(0.00045) << endl;

//figures.plot_Z();

double Y=relic_density.Y_today(relic_density.x_f());//cross_section.x_f());
double rho_crit=1.05375e-5;// 3*(pow(H_0,2))/(8*Pi*G);  units: h^2 GeV /( c^2 cm^3)
double s_0=2890;//(2*pow(Pi,2)/45)*g*pow(T,3);
cout << "mass fraction = " <<  Y*m_s*s_0/rho_crit << endl;

}