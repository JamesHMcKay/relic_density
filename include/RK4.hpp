#ifndef RK4_H
#define RK4_H


#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

#include "rd.hpp"

using namespace std;



struct simple_DE {
simple_DE(double t) : t(t) {}
//simple(){}
double operator()(double x,double y) const { return t*y*sin(x);}

private:
double t;
};

template<class T>
struct RK4 {
T &func;
double xi, xf, yi;

RK4() {}; // constructor
RK4(T &funcc, const double xii, const double xff, const double yii ) :
func ( funcc) , xi(xii) , xf(xff), yi(yii) {}


double RK4_step(double x,double Y, double h)
{

double k1 =h* func(x,Y);

double k2 = h* func(x+0.5*h, Y+0.5*k1);

double k3 = h* func(x+0.5*h, Y+0.5*k2);

double k4 = h* func(x+h,  Y+k3 ) ;

return Y + 0.166666667 * (k1+2*k2+2*k3+k4);
}

double solve_RK4()
{
ofstream myfile;
myfile.open ("../Figures/data/RK4.txt");
// the ODE to solve is dY/dx = F(x,Y) where F(Y,x) = Z(x) * (Y_eq^2(x)-Y^2(x))
double tol=1e-6;
double x = xi;
double Y=yi;
double h=0.5;
double k1,k2,k3,k4;
double Yf;
double Y_1=Y;
while ( x < xf )
{

double Y_orig=Y;
double x_orig=x;

// perform one step of RK4 method

Y_1 = RK4_step(x,Y,h);

// perform two steps of RK4 method

h=0.5*h;

Y =  RK4_step(x,Y,h);
x = x + h;

double Y_2 = RK4_step(x,Y,h);
x = x+h; // x is now at next step value
// end of second step of RK4
h=h*2;

// now need to adjust step size
double error = (Y_2 - Y_1);


h = 0.9 * (h) * pow(tol / abs(error),0.2);  // multiply h by two since we halved it above

cout << "error = " << error << " h = " << h << endl;

if (abs(error)< tol)
{
//// step successful, move onto next point
//cout << "step successful" << Y << endl;

Y=Y_2+error/15.0;
myfile << x << " " << Y << " " << Y_1 << endl;
// leave x as is, it's already at the next point
}
else
{
// step not successful, go back and use new h
Y=Y_orig;
x=x_orig; // use saved value, since h has changed
}

Yf=Y;

}
myfile.close();


return Yf;


}
  
  

};








#endif