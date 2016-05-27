#include "rd.hpp"
#include "interp.hpp"
#include "figures_rd.hpp"
#include "integrate.hpp"
#include "data.hpp"
#include "models.hpp"
#include "rd.cpp"
#include "figures_rd.cpp"
#include "RK4.hpp"
#include "mcmc.hpp"

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

using namespace std;


double mass_frac(Data data)
{
double Y;
Relic_density<SingletDM_RD> relic_density(data);
try
{
Y = relic_density.Y_today(relic_density.x_f());
}
catch (const char* msg)
{
cout << msg << endl;
return 0;
}
double rho_crit = 1.05375e-5;
double s_0=2890;
return  Y*data.M_s*s_0/rho_crit;
}




void plot_rd()
{
int pts=10;
int x_pts=pts,y_pts=pts;
double m_lower=50, m_upper=70;
double lam_lower_log=-4, lam_upper_log=0;
Data data;

ofstream myfile;
myfile.open ("../Figures/data/relic_density.txt");


double lam_step= (lam_upper_log-lam_lower_log)/y_pts;

double m_step= (m_upper-m_lower)/x_pts;

for (int i=0; i<pts ;i++)
{
 data.Lam_hs=pow(10,-lam_step*float(i));
 
  for (int j=0; j<pts ; j++)
  {
   data.M_s=m_step*float(j)+m_lower;
   
   cout << "calculating point " << data.Lam_hs << " " << data.M_s << endl;
   
   
   myfile << data.Lam_hs << " " << data.M_s << " " << mass_frac(data) << endl;
   
  }
}
myfile.close();


}



int main(int argc, char* argv[])
{

Data data(argc,argv);

//plot_rd();

//Figures_RD<MDM_RD> figures(data);
//figures.plot_Z();

//figures.plot_thermal_av();

//figures.plot_sigma_v(10);

//figures.plot_Z();

//data.method=0;
cout << "Oh2 (method 0) = " << mass_frac(data) << endl;
//data.method=1;
//cout << "Oh2 (method 1) = " << mass_frac(data) << endl;

// test mcmc scanner

MCMC mcmc(data);

mcmc.scanner();



// test the RK4 solver with simple test function
//simple_DE F(3.0);
//RK4<simple_DE> rk(F,1.0,20.0,2.0);
//cout << "solution is = " << rk.solve_RK4() << endl;
}