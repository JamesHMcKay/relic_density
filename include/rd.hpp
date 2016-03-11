#ifndef RD_H
#define RD_H

#include "interp.hpp"
#include "integrate.hpp"
#include "singletdm.hpp"
#include "data.hpp"

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

using namespace std;

class Relic_density
{

  private:

  int make_interp;
  Data data;
  public:
  std::vector<double> T_m,thermal_av_m;
  Relic_density (){make_interp=0;}  // defualt constructor
  Relic_density(Data _data)
  {
      data=_data;
      make_interp=0;
   } //constructor
//  double Z(double x);
  
  double Yeq(double x);
  
  double dYeq(double x);
  
  double hat(double func,double x);
  
  double x_f();
  
  double Y_today(double x_f);
  
  double A(double x_f);
  
  double calc_cross_section(double T);
  
  void thermal_average_make_interp(double T_lower, double T_upper,int pts);


};


struct Z_func {

Z_func(Data data) : data(data) {use_interp=0;}

Z_func(Data data ,std::vector<double> T,std::vector<double> thermal_av) : data(data)
,T_m(T), thermal_av_m(thermal_av) {use_interp=1;}

double operator()(double x) {


double g_eff=10,cs;   /// check this !!!


if (use_interp==1)
{
Poly_interp myfunc(T_m,thermal_av_m,4);
cs = myfunc.interp(data.M_s/x);
}
else
{
Relic_density rd(data);
cs=rd.calc_cross_section(data.M_s/x);
}
return pow(data.Pi/float(45),0.5)*((data.M_s*data.M_pl)/pow(x,2)) * (pow(g_eff,0.5))*cs;

}//cs_integral(x,T);}

private:
double T;
Data data;
int use_interp;
std::vector<double> T_m, thermal_av_m;
};





#endif
