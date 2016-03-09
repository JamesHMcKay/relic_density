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



int main(int argc, char* argv[])
{
double lambda_hs=0.1,M_s=100;double param [2];

std::string name [2]; int i=0;

if (argc==1)
{
cout << "Please enter a data file in the format ./main input.txt " << endl;
return 0;
}
else
{
std::ifstream input(argv[1]);
std::string line;


while(getline(input, line)) {
      if (!line.length() || line[0] == '#')
         continue;
      std::istringstream iss(line);
      iss>> name[i] >> param[i];
      i=i+1;
   }
}

for (int n=0;n<i+1;n++)
{
if (name[n]=="lambda_hs")
{
lambda_hs=param[n];
}
if (name[n]=="M_s")
{
M_s=param[n];
}
}

cout << "lambda_hs = " << lambda_hs << endl;
cout << "M_s = " << M_s << endl;




//Figures figures(cross_section, relic_density);
//figures.plot_Z();

//figures.plot_thermal_av();

//figures.plot_sigma_v(0.00045);

//figures.plot_Z();


cout << "Oh2 = " << mass_frac(lambda_hs,M_s) << endl;

}