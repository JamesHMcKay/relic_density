#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

#include "models.hpp"
#include "interp.hpp"
#include "bessel.hpp"

using namespace std;

///  Cross section integral class functions //


double MDM_RD::cs_integral(double s, double T)
{

double Pi=data.Pi;
double g2=data.g2;
double M=data.M_chi;

Bessel bessel;// class for asymptotic bessel function definitions


double cp = (float(1215)/(8*Pi))*pow(g2,4);
double cs = (float(1035)/(8*Pi))*pow(g2,4);

double x=s/pow(M,2);
double beta=pow(1-4/x,0.5);



double A =exp( bessel.log_bessel_k_asymtotic(1,pow(s,0.5)/T) + log (pow(s,0.5)*(T/(64*pow(Pi,4))))  );

//cout << "A = " << A << endl;

//cout << "beta = " << beta << endl;

return (cs*beta+cp*pow(beta,3))*A;

}



double MDM_RD::quick_cs(double T)
{
double Pi=data.Pi;
double g2=data.g2;
double M=data.M_chi;

Bessel bessel;// class for asymptotic bessel function definitions

double cp = (float(1215)/(8*Pi))*pow(g2,4);
double cs = (float(1035)/(8*Pi))*pow(g2,4);

double result = ( cs + 0.5*(3*T/M)*(cp+0.5*cs) );

double A = (M*pow(T,3)*exp(-2*M/T)) / (32*pow(Pi,3));

return result*A;

}