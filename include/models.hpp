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

class SingletDM
{
  private:
  Data data;
  public:
  SingletDM (){}  // defualt constructor
  SingletDM(Data _data)
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
  
  
};

class SingletDMZ3
{
  private:
  Data data;
  public:
  SingletDMZ3 (){}  // defualt constructor
  SingletDMZ3(Data _data)
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
