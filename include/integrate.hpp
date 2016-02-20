#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "interp.hpp"

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

#include "cross_section.hpp"

using namespace std;


struct simple {
simple(double t) : t(t) {}
//simple(){}
double operator()(double x) const { return t*exp(-x);}

private:
double t;
};


struct Quadrature{
//abstract base class
int n; // current level or refinement
virtual double next() = 0 ;
};

template<class T>
struct Trapzd : Quadrature {
T &func;
double a, b , s;

Trapzd() {}; // constructor
Trapzd(T &funcc, const double aa, const double bb) :
func ( funcc) , a(aa) , b(bb) {n=0;}

double next()
{
  double x, tnm, sum , del;
  int it,j;
  n++;
  if (n==1)
  {
    return (s=0.5*(b-a)*(func(a)+func(b)));
  }
  else
  {
    for (it=1, j=1; j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0, j=0;j<it;j++,x+=del) sum +=func(x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
  
}

};


template<class T>
double qtrap(T &func, const double a, const double b, const double eps=1.0e-10)
{
const int JMAX=30;
double s,olds=0.0;
Trapzd<T> t(func,a,b);
for (int j=0;j<JMAX;j++) {
    s=t.next();
    if (j > 5)
if (abs(s-olds) < eps*abs(olds) ||
(s == 0.0 && olds == 0.0)) return s;
olds=s; }
cout<< "too many steps " << endl;//throw("Too many steps in routine qtrap");
}








#endif