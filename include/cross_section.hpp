#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include "interp.hpp"
#include "integrate.hpp"

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

using namespace std;

class Cross_section
{
  private:
  double Mh;

  double lambda_hs,v_0;
  double Pi=3.14;
  double M_Z=91, M_W=80 , alpha_s=0.1184, mt=173;
  public:
  double m_s;
  Cross_section (){}  // defualt constructor
  Cross_section(double Ms)
  {
      m_s=Ms;
      lambda_hs=0.1;
      v_0=246;
      Mh=125.66;
   } //constructor

  
  double calc_cross_section(double T);
  
  double gamma_h(double x);
  
  double Dh2(double s);
  
  
  double sigma_v_VV(double s);
  
  double sigma_v_ff(double s);
  
  double sigma_v_hh(double s);
  
  double cs_integral(double s, double T);
//
//  double integrator(double lower, double upper, double T); // use logarithmic spacing
//  double integrator_2(double lower, double upper, double T); // linear spacing

  double sigma_v(double s);
  
  double bessel_k_asymtotic(int v, double z)
  { // approximation for asymptotic Bessel function from http://userpages.umbc.edu/~dfrey1/ench630/modbes.pdf
   double mu=4*pow(v,2);
    
   double result = pow(Pi/(2*z),0.5)*exp(-z) * ( 1 + (mu-1)/( float(8) * z) + (mu-1)*(mu-9)/ ( 2 * pow( 8*z, 2 ) ) + (mu-1)*(mu-9)*(mu-25)/ ( 3*2*pow(8*z,3) ) );
  
   return result;
  };
  
  
    double log_bessel_k_asymtotic(int v, double z)
  {
  
   if (z>200)
   {
   return 0.5*log(Pi/(2*z))-z;
   }
   else
   {
   return log(boost::math::cyl_bessel_k(v,z));
   }
  };
  
};



struct cs_func {
cs_func(double T, double Ms) : T(T) , Ms(Ms) {}
double operator()(double x) const {
Cross_section cs(Ms);
return cs.cs_integral(x,T);
}//cs_integral(x,T);}

private:
double T, Ms;
};



class Relic_density
{

  private:
  double Mh;

  double lambda_hs,v_0;
  double Pi=3.14;
  double M_Z=91, M_W=80 , alpha_s=0.1184, mt=173;
  public:
  double m_s;
  Relic_density (){}  // defualt constructor
  Relic_density(double Ms)
  {
      m_s=Ms;
      lambda_hs=0.1;
      v_0=246;
      Mh=125.66;
   } //constructor



  double Z(double x);
  
  double Yeq(double x);
  
  double dYeq(double x);
  
  double hat(double func,double x);
  
  double x_f();
  
  double Y_today(double x_f);
  
  double A(double x_f);
  
  double integrator_A(double lower, double upper);
  
  double integrator_A_2(double lower, double upper);
  
  double calc_cross_section(double T)
  {
  double s=4*pow(m_s,2);

  cs_func func(T,m_s);

  return qtrap(func, s,1e6,1e-2);

  
  };

};



#endif
