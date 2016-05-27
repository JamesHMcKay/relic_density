#ifndef RD_H
#define RD_H

#include "interp.hpp"
#include "integrate.hpp"
#include "models.hpp"
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


template<typename Model>
class Relic_density
{

  private:
  int make_interp=0;
  Data data;
  bool fast_mode;
  int method=1; // default approximation method for solving Boltzman equation
  public:
  std::vector<double> T_m,thermal_av_m;
  Relic_density (){}  // defualt constructor
  Relic_density(Data _data)
  {
      fast_mode=_data.fast_mode;
      method = _data.method;
      Model model;
      data = model.set_dm_mass(_data);
   } //constructor
  double Z(double x);
  
  double Yeq(double x);
  
  double dYeq(double x);
  
  double hat(double func,double x);
  
  double x_f();
  
  double Y_today(double x_f);
  
  double A(double x_f);
  
  double calc_cross_section(double T);

  double calc_cross_section_fast(double T);
  
  double calc_cross_section_slow(double T);
  
  double solve_boltzman_Y_today(double x_f);
  
  void thermal_average_make_interp(double T_lower, double T_upper,int pts);

};

template <class Model>
struct Z_func {

Z_func(Data data) : data(data) {use_interp=0;}

Z_func(Data data ,std::vector<double> T,std::vector<double> thermal_av) : data(data)
,T_m(T), thermal_av_m(thermal_av) {use_interp=1;}

double operator()(double x) {


double g_eff=10,cs;   /// check this !!! need to use a proper g effective from somewhere


if (use_interp==1)
{
Poly_interp myfunc(T_m,thermal_av_m,4);
cs = myfunc.interp(data.M_dm/x);
}
else
{
Relic_density<Model> rd(data);
cs=rd.calc_cross_section(data.M_dm/x);
}
return pow(data.Pi/float(45),0.5)*((data.M_dm*data.M_pl)/pow(x,2)) * (pow(g_eff,0.5))*cs;

}

private:
double T;
Data data;
int use_interp;
std::vector<double> T_m, thermal_av_m;
};



template <class Model>
struct Boltz_func {

Boltz_func(Data data ,std::vector<double> T,std::vector<double> thermal_av) : data(data)
,T_m(T), thermal_av_m(thermal_av) {}

double operator()(double x,double Y) {


double g_eff=10,cs;   /// check this !!! need to use a proper g effective from somewhere

Relic_density<Model> rd(data);

Poly_interp myfunc(T_m,thermal_av_m,4);
cs = myfunc.interp(data.M_dm/x);
double Zx =  pow(data.Pi/float(45),0.5)*((data.M_dm*data.M_pl)/pow(x,2)) * (pow(g_eff,0.5))*cs;

if (data.alpha == 0)
{
return Zx * ( pow(rd.Yeq(x),2) - pow(Y,2) );
}
else
{
return Zx * ( (1-data.alpha)*pow(rd.Yeq(x),2) + data.alpha*Y*rd.Yeq(x) - pow(Y,2) );
}
}

private:
double T;
Data data;
std::vector<double> T_m, thermal_av_m;
};





#endif
