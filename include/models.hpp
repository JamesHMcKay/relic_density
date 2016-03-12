#ifndef MODELS_H
#define MODELS_H

#include "data.hpp"
#include "interp.hpp"
#include "integrate.hpp"
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

using namespace std;

class SingletDM_RD
{
  private:
  Data data;
  public:
  SingletDM_RD (){}  // defualt constructor
  SingletDM_RD(Data _data)
  {
    data=_data;
   } //constructor
  
  double gamma_h(double x);
  
  double Dh2(double s);
  
  double sigma_v_VV(double s);
  
  double sigma_v_ff(double s);
  
  double sigma_v_hh(double s);
  
  double cs_integral(double s, double T);

  double sigma_v(double s);
  
  double quick_cs(double T){return 1e-8;}//(1e-12)*(1e-28)*(3e10)*(10*10);};
  
  Data set_dm_mass(Data data){data.M_dm=data.M_s;return data;};
  
};

class SingletDMZ3_RD
{
  private:
  Data data;
  public:
  SingletDMZ3_RD (){}  // defualt constructor
  SingletDMZ3_RD(Data _data)
  {
    data=_data;
   } //constructor
  
  double gamma_h(double x);
  
  double Dh2(double s);
  
  double sigma_v_VV(double s);
  
  double sigma_v_ff(double s);
  
  double sigma_v_hh(double s);
  
  double cs_integral(double s, double T);

  double sigma_v(double s);
  
  double quick_cs(double T){return (1e-12)*(1e-28)*(3e8);};

  Data set_dm_mass(Data data){data.M_dm=data.M_s;return data;};
  
  
};


class MDM_RD
{
  private:
  Data data;
  public:
  MDM_RD (){}  // defualt constructor
  MDM_RD(Data _data)
  {
    data=_data;
   } //constructor
  
  double quick_cs(double T);
  
  double cs_integral(double s, double T);

  Data set_dm_mass(Data data){data.M_dm=data.M_chi;return data;};
  
  
};


template <class M>
struct cs_func {
cs_func(double T, Data data) : T(T) , data(data) {}
double operator()(double x) const {
M cs(data);
return cs.cs_integral(x,T);
}//cs_integral(x,T);}

private:
double T;
Data data;
};



#endif
