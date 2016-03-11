#ifndef BESSEL_H
#define BESSEL_H

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

using namespace std;

class Bessel
{
  private:
  double Pi=3.14;
  public:
  Bessel (){}
  
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
#endif
