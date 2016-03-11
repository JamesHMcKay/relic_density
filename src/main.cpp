#include "rd.hpp"
#include "interp.hpp"
//#include "figures.hpp"
#include "integrate.hpp"
#include "data.hpp"
#include "models.hpp"
#include "rd.cpp"

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

using namespace std;


double mass_frac(Data data)
{
Relic_density<SingletDM> relic_density(data);
double Y=relic_density.Y_today(relic_density.x_f());
double rho_crit=1.05375e-5;
double s_0=2890;
return  Y*data.M_s*s_0/rho_crit;
}




void plot_rd()
{
int x_pts=20,y_pts=20;
double m_lower=50, m_upper=70;
double lam_lower_log=-4, lam_upper_log=0;
Data data;

ofstream myfile;
myfile.open ("../Figures/data/relic_density.txt");


double lam_step= (lam_upper_log-lam_lower_log)/y_pts;

double m_step= (m_upper-m_lower)/x_pts;

for (int i=0; i<20 ;i++)
{
 data.Lam_hs=pow(10,-lam_step*float(i));
 
  for (int j=0; j<20 ; j++)
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

//Figures figures(cross_section, relic_density,data);
//figures.plot_Z();

//figures.plot_thermal_av();

//figures.plot_sigma_v(0.00045);

//figures.plot_Z();

cout << "Oh2 = " << mass_frac(data) << endl;

}