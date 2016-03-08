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

  double v_0;
  double Pi=3.14;
  double M_Z=91, M_W=80 , alpha_s=0.1184, mt=173;
  public:
  double m_s,lambda_hs;
  Cross_section (){}  // defualt constructor
  Cross_section(double Ms,double Lam_hs)
  {
      m_s=Ms;
      lambda_hs=Lam_hs;
      v_0=246;
      Mh=125.66;
   } //constructor
  
  double gamma_h(double x);
  
  double Dh2(double s);
  
  double sigma_v_VV(double s);
  
  double sigma_v_ff(double s);
  
  double sigma_v_hh(double s);
  
  double cs_integral(double s, double T);

  double sigma_v(double s);
  
  double bessel_k_asymtotic(int v, double z)
  { // approximation for asymptotic Bessel function from http://userpages.umbc.edu/~dfrey1/ench630/modbes.pdf
   double mu=4*pow(v,2);
   double result;
   if (z>200)
   {
    
   result = pow(Pi/(2*z),0.5)*exp(-z) * ( 1 + (mu-1)/( float(8) * z) + (mu-1)*(mu-9)/ ( 2 * pow( 8*z, 2 ) ) + (mu-1)*(mu-9)*(mu-25)/ ( 3*2*pow(8*z,3) ) );
   }
   else
   {
   result = boost::math::cyl_bessel_k(v,z);
   }
   
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
cs_func(double T, double Ms,double lambda_hs) : T(T) , Ms(Ms) , lambda_hs(lambda_hs) {}
double operator()(double x) const {
Cross_section cs(Ms,lambda_hs);
return cs.cs_integral(x,T);
}//cs_integral(x,T);}

private:
double T, Ms, lambda_hs;
};



class Relic_density
{

  private:
  double Mh;

  int make_interp;
  double lambda_hs,v_0;
  double Pi=3.14;
  double M_Z=91, M_W=80 , alpha_s=0.1184, mt=173;
  public:
  double m_s;
  std::vector<double> T_m,thermal_av_m;
  Relic_density (){make_interp=0;}  // defualt constructor
  Relic_density(double Ms,double Lam_hs)
  {
      m_s=Ms;
      lambda_hs=Lam_hs;
      v_0=246;
      Mh=125.66;
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

Z_func(double Ms,double lambda_hs) : Ms(Ms), lambda_hs(lambda_hs) {use_interp=0;}

Z_func(double Ms ,double lambda_hs ,std::vector<double> T,std::vector<double> thermal_av) : Ms(Ms),lambda_hs(lambda_hs)
,T_m(T), thermal_av_m(thermal_av) {use_interp=1;}

double operator()(double x) {

double M_pl=1.2e19;
double g_eff=10,cs;


if (use_interp==1)
{
Poly_interp myfunc(T_m,thermal_av_m,4);
cs = myfunc.interp(Ms/x);
}
else
{
Relic_density rd(Ms,lambda_hs);
cs=rd.calc_cross_section(Ms/x);
}

double Pi=3.14;
return pow(Pi/float(45),0.5)*((Ms*M_pl)/pow(x,2)) * (pow(g_eff,0.5))*cs;

}//cs_integral(x,T);}

private:
double T, Ms,lambda_hs;
int use_interp;
std::vector<double> T_m, thermal_av_m;
};





#endif
