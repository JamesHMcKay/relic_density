#ifndef INTERP_H
#define INTERP_H

#include <vector>
#include <complex>
#include <iostream>
using namespace std;

struct Base_interp
{
  int n,mm, jsav, cor, dj;
  const double *xx, *yy;
  
  Base_interp(std::vector<double> &x, const double *y, int m) // constructor
  : n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y)
  {
  
  dj = MIN(1,(int)pow((double)n,0.25));
  }
  
  template<class T>
  inline const T &MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}

  inline float MIN(const double &a, const float &b)
        {return b < a ? (b) : float(a);}

  inline float MIN(const float &a, const double &b)
        {return b < a ? float(b) : (a);}
  
  
  template<class T>
  inline const T &MAX(const T &a, const T &b)
          {return b > a ? (b) : (a);}

  inline float MAX(const double &a, const float &b)
          {return b > a ? (b) : float(a);}

  inline float MAX(const float &a, const double &b)
          {return b > a ? float(b) : (a);}


  
  
  double interp(double x)
  {
    int jlo = cor ? hunt(x) : locate(x);
    return rawinterp(jlo, x);
  }
  
  int locate(const double x);
  int hunt(const double x);
  
  double virtual rawinterp(int jlo, double x) = 0; // derived class, provide this as the interpolation class


  
  
};


struct Linear_interp : Base_interp
//Piecewise linear interpolation object. Construct with x and y vectors, then call interp for interpolated values.
{
Linear_interp(std::vector<double> &xv, std::vector<double> &yv) : Base_interp(xv,&yv[0],2) {}
double rawinterp(int j, double x) {
if (xx[j]==xx[j+1]) return yy[j]; //Table is defective, but we can recover.
else return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
}
};


struct Poly_interp : Base_interp
//Polynomial interpolation object. Construct with x and y vectors, and the number M of points to be used locally (polynomial order plus one), then call interp for interpolated values.
{
double dy;
Poly_interp(std::vector<double> &xv, std::vector<double> &yv, int m)
: Base_interp(xv,&yv[0],m), dy(0.) {} double rawinterp(int jl, double x);
};

#endif













  
  