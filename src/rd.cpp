#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <fstream>

#include "rd.hpp"
#include "interp.hpp"
#include "singletdm.hpp"


using namespace std;


template<typename Model>
double Relic_density<Model>::Yeq(double x)
{
double h_eff=100;
return boost::math::cyl_bessel_k(2,x)*(45/(4*pow(data.Pi,4)))*pow(x,2)/(h_eff);
}

template<typename Model>
double Relic_density<Model>::dYeq(double x)
{
double h_eff=100;
return boost::math::cyl_bessel_k_prime(2,x)*(45/(4*pow(data.Pi,4)))*pow(x,2)/(h_eff)  +  (float(2)/x)*Yeq(x);
}


template<typename Model>
double Relic_density<Model>::hat(double func, double x)
{
return exp(x)*func;
}

template<typename Model>
double Relic_density<Model>::x_f()
{
// calculate x value of free-out

double deltaf=0.618;
double x_star=15,x_prev=0;
Z_func<Model> func(data);



for (int i=0;i<10;i++)

while (abs(x_star-x_prev)>0.1)
{
x_star=log( (deltaf*(2+deltaf)/(1+deltaf)) * (  func(x_star)* pow(hat(Yeq(x_star),x_star),2) )/ (hat(Yeq(x_star),x_star) -hat(Yeq(x_star)+dYeq(x_star),x_star) ));
x_prev=x_star;
//cout<< "x* is = " << x_star << endl;
//cout<< "corresponing to T = " << m_s/x_star << endl;
}

cout<< "x* is = " << x_star << endl;

return x_star;
}



template<typename Model>
double Relic_density<Model>::Y_today(double x_f)
{
double Y_f=(1+0.618)*Yeq(x_f);


double result=Y_f/(1+Y_f*A(x_f));

return result;
}

template<typename Model>
double Relic_density<Model>::A(double x_f)
{
  double x_upper=1e5, x_lower=x_f;
  double T_lower=data.M_s/x_upper,T_upper=data.M_s/x_lower;
  
  
  thermal_average_make_interp(T_lower*0.9,T_upper*1.1,5);

  Z_func<Model> func(data,T_m,thermal_av_m);
  
  double I=qtrap(func, x_lower,x_upper,1e-3);
  
  return I;

//Z_func func(m_s);
//
//return x_f*func(x_f);

}

template<typename Model>
double Relic_density<Model>::calc_cross_section(double T)
{
  if (make_interp==1)
  {
  
  Poly_interp myfunc(T_m,thermal_av_m,5);

  return myfunc.interp(T);
  }
  else
  {
  // find ul such that f(ul)=1e-30
  
  double s=4*pow(data.M_s,2)+1;
  cs_func<Model> f(T,data);
  double a=s*10.0;
  double delta=1;
  double fa=f(a);
  double val=1e-30;
  double fdelta=f(a+delta);
  double diff=fa-fdelta;
  double b,bb=-1;
    if (diff<0)
    {
    
    cout << " gradient is positive, going over the hill " << endl;
    
    while (diff<0)
    {
    a=a*10;
    fa=f(a);
    fdelta=f(a+delta);
    diff=fa-fdelta;
    }
    }

    if (fa>val)
    {
//      cout << " fa is too large, descending downhill " << endl;
      if (f(a+delta)<val)
      {
      // desired point is between a and delta, so set b = a + delta
      bb=a+delta;
      }
      else
      {
//        cout << "looking for f < val "<< endl;
      // need to go beyond a+delta to find point
      // here we will find a value such that f < val
        double b=delta+a,fb=f(a+delta);
        while (fb>val)  // find a and b around root
        {
        a=b;
        b=b*10;
        bb=b;
//        cout << "attempt = " << b << " f = " << f(b) << endl;
        fb=f(b);
        }
        
      }
    }
    else
    {
//      cout << " fa is too small, going uphill" << endl;
      if (f(a-delta)>val)
      {
      // desired point is between a and delta, so set b = a + delta
      bb=a+delta;
      }
      else
      {
      // need to go below a-delta to find point
      // but we our bounded below by s so just use this
      bb=a;
      a=s;
      // now have lower (a) and upper (b) bound for x such that f(x)=val
      }
    }
    b=bb;
 //   cout << "root bounded between (a,b) = " << a << " " << b << endl;
 //   cout << "root bounded between (f(a),f(b)) = " << f(a) << " " << f(b) << endl;
  

    fa=f(a);
    double ft=val*1e10,m;
    // now root find between a and b
    while (abs(log10(ft)-log10(val))>5)
    {
      m=(b-a)/2+a;
//      cout << " mid point = " << m << endl;
      if (f(m)<val)
      {
      b=m;
      }
      else
      {
        a=m;
      }
      double fm=f(m);
      if (fm!=0)
      {
      ft=fm;
      }
    }

   double ul=m;
  

//   cout << "integrating from s to " << ul << " for T = " << T << endl;
  
  double result;
  
//  if (s<2*pow(Mh,2))  // split integral into two parts, one around the Higgs resonance and one for the tail
//  {
//  double mid_pt=2*pow(Mh,2);
//  cout << " (s, mid) = " << s << " " << mid_pt << endl;
//  
//  double int_lower=qtrap(f,s,mid_pt,1e-4);
//  double int_upper=qtrap(f,mid_pt,ul,1e-4);
//  result=int_lower+int_upper;
//  }
//  else
//  {
//  result=qtrap(f,s,ul,1e-4);
//  }
   result=qtrap(f,s,ul,1e-4);
  
  return result;
  }
}

template<typename Model>
void Relic_density<Model>::thermal_average_make_interp(double T_lower, double T_upper,int pts)
{
  std::vector<double> T(pts);
  std::vector<double> thermal_av(pts);
  double step=(T_upper-T_lower)/float(pts);
  T[0]=T_lower;
  for (int n=0;n<pts;n++)
  {
  thermal_av[n]=calc_cross_section(T[n]);
  T[n+1]=T[n]+step;
  }
  T_m=T;
  thermal_av_m=thermal_av;
  make_interp=1;
}


