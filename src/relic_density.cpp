#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <fstream>

#include "cross_section.hpp"
#include "interp.hpp"

using namespace std;



//double Relic_density::Z(double x)
//{
//double M_pl=1e19;
//double g_eff=100;
//return pow(Pi/float(45),0.5)*((m_s*M_pl)/pow(x,2)) * (pow(g_eff,0.5))*(calc_cross_section(m_s/x));
//}

double Relic_density::Yeq(double x)
{
double h_eff=100;
return boost::math::cyl_bessel_k(2,x)*(45/(4*pow(Pi,4)))*pow(x,2)/(h_eff);
}

double Relic_density::dYeq(double x)
{
double h_eff=100;
return boost::math::cyl_bessel_k_prime(2,x)*(45/(4*pow(Pi,4)))*pow(x,2)/(h_eff)  +  (float(2)/x)*Yeq(x);
}

double Relic_density::hat(double func, double x)
{
return exp(x)*func;
}


double Relic_density::x_f()
{
// calculate x value of free-out

double deltaf=0.618;
double x_star=0.1;
Z_func func(m_s);



for (int i=0;i<10;i++)
{
x_star=log( (deltaf*(2+deltaf)/(1+deltaf)) * (  func(x_star)* pow(hat(Yeq(x_star),x_star),2) )/ (hat(Yeq(x_star),x_star) -hat(Yeq(x_star)+dYeq(x_star),x_star) ));
//cout<< "x* is = " << x_star << endl;
//cout<< "corresponing to T = " << m_s/x_star << endl;
}

cout<< "x* is = " << x_star << endl;

return x_star;
}





double Relic_density::Y_today(double x_f)
{
double Y_f=(1+0.618)*Yeq(x_f);


double result=Y_f/(1+Y_f*A(x_f));

return result;
}


double Relic_density::A(double x_f)
{
  double x_upper=1e5, x_lower=x_f;
  double T_lower=m_s/x_upper,T_upper=m_s/x_lower;
  
  
  thermal_average_make_interp(T_lower*0.9,T_upper*1.1,5);

  Z_func func(m_s,T_m,thermal_av_m);
  
  double I=qtrap(func, x_lower,x_upper,1e-3);
  
  return I;

//Z_func func(m_s);
//
//return x_f*func(x_f);

}


double Relic_density::calc_cross_section(double T)
{
  if (make_interp==1)
  {
  
  Poly_interp myfunc(T_m,thermal_av_m,4);

  return myfunc.interp(T);
  }
  else
  {

  double s=4*pow(m_s,2)+1;
  cs_func func(T,m_s);
  double ul=1e8;

// first need to determine appropriate upper limit for integration range, such that function is non-zero there
  if (func(s)==0){ return 0;}
  
  while (func(ul)==0)
  {
  ul=(ul-s)/2+s;
  //cout << "ul = " << ul << " f = " << func(ul) << endl;
  }
   cout << "integrating from s to " << ul << "for T = " << T << endl;
  
  double result;
  
  ul=112617;
  
  if (s<2*pow(Mh,2))  // split integral into two parts, one around the Higgs resonance and one for the tail
  {
  double mid_pt=2*pow(Mh,2);
  cout << " (s, mid) = " << s << " " << mid_pt << endl;
  
  double int_lower=qtrap(func,s,mid_pt,1e-4);
  double int_upper=qtrap(func,mid_pt,ul,1e-4);
  result=int_lower+int_upper;
  }
  else
  {
  result=qtrap(func,s,ul,1e-4);
  }
  
  
  return result;
  }
}


void Relic_density::thermal_average_make_interp(double T_lower, double T_upper,int pts)
{
  std::vector<double> T(pts);
  std::vector<double> thermal_av(pts);
  double step=(T_upper-T_lower)/float(pts);
  T[0]=T_lower;
  for (int n=0;n<pts;n++)
  {
  thermal_av[n]=calc_cross_section(T[n]);
  T[n+1]=T[n]+step;
  cout << "T = " << T[n] << "thermal av = " << thermal_av[n] << endl;

  }
  
  T_m=T;
  thermal_av_m=thermal_av;
  make_interp=1;
}


