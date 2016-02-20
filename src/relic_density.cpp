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



double Relic_density::Z(double x)
{
double M_pl=1e19;
double g_eff=100;
return pow(Pi/float(45),0.5)*((m_s*M_pl)/pow(x,2)) * (pow(g_eff,0.5))*(calc_cross_section(m_s/x));
}

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

for (int i=0;i<10;i++)
{
x_star=log( (deltaf*(2+deltaf)/(1+deltaf)) * (  Z(x_star)* pow(hat(Yeq(x_star),x_star),2) )/ (hat(Yeq(x_star),x_star) -hat(Yeq(x_star)+dYeq(x_star),x_star) ));

cout<< "x* is = " << x_star << endl;
cout<< "corresponing to T = " << m_s/x_star << endl;
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

double integrand=1;
double p=0,s;
double lower_cutoff=x_f;
while (integrand>0)
{
p=p+1;
s=pow(10,p)+lower_cutoff;
//cout<< "s = " << s << endl;
integrand=Z(s);

cout << "Z = " << integrand << "p= " << p << endl;
}
cout << "upper limit set = " << s << endl;

return integrator_A(600,1e12)+integrator_A_2(x_f,600);

}



double Relic_density::integrator_A(double lower, double upper)
{
    double fa,fb,l=0,h,l_old=0;
    double e=1;
    int it;
    it=0;
    double a=lower,b=upper;
    fa=Z(lower);
    fb=Z(upper);
 //   cout << "fa = " << fa << endl;
  //  cout << "fb = " << fb << endl;
    h=(b-a)/2;
    l=(h/2)*(fa+fb);
    //printf ("First estimate is %f \n",l);
    //double midpoint=(a+b)/2;
    int n;
    n=1;
    int pts;
    pts=100;
    int s;
    double c,c_prev;
    s=1;
  
    while(e>0.3)
    {
        l_old=l;
        
       // printf ("Iteration number %i \n",it);

        //l=l/h;
        //h=(b-a)/(pts-1); // h is the length of each section in the integration
        //l=l*h;
        l=0;
        c_prev=log10(a);
        for (int m=1; m<2*s+1; m +=2)
        {
            c=m*(log10(b)-log10(a))/(2*s)+log10(a);
            h=pow(10,c)-pow(10,c_prev);
            l=l+h*Z(pow(10,c));
            c_prev=c;
          
            //cout<< "sampling at " << pow(10,c) << " obtained value = " << cs_integral(pow(10,c),T) << " h = " << h << endl;
          

        }
      //  printf ("updated value of integral %f \n",l);
        cout << "updated value of integral = " << l << "prev  = " << l_old << endl;
        pts=((pts-2)*2)+pts;
        n=n+1;
        s=pow(2,n-1);
        pts=2*s+1;
       // printf ("number of points %i \n",2*s+1);
        e=abs((l-l_old)/l);
        it=it+1;
    }  // set the desired accuracy here <<<<<<<<<<<<<
    
    
    
    return l;
}


double Relic_density::integrator_A_2(double lower, double upper)
{
    double fa,fb,l=0,h,l_old=0;
    double e=1;
    int it;
    it=0;
    double a=lower,b=upper;
    fa=Z(lower);
    fb=Z(upper);
 //   cout << "fa = " << fa << endl;
  //  cout << "fb = " << fb << endl;
    h=(b-a)/2;
    l=(h/2)*(fa+fb);
    //printf ("First estimate is %f \n",l);
    //double midpoint=(a+b)/2;
    int n;
    n=1;
    int pts;
    pts=100;
    int s;
    double c,c_prev;
    s=1;
  
    while(e>0.1)
    {
        l_old=l;
        
       // printf ("Iteration number %i \n",it);

        //l=l/h;
        //h=(b-a)/(pts-1); // h is the length of each section in the integration
        //l=l*h;
        l=0;
        c_prev=(a);
        for (int m=1; m<2*s+1; m +=2)
        {
            c=m*((b)-(a))/(2*s)+(a);
            h=c-c_prev;//pow(10,c)-pow(10,c_prev);
            l=l+h*Z(c);
            c_prev=c;
          
            //cout<< "sampling at " << pow(10,c) << " obtained value = " << cs_integral(pow(10,c),T) << " h = " << h << endl;
          

        }
      //  printf ("updated value of integral %f \n",l);
        cout << "updated value of integral = " << l << "prev  = " << l_old << endl;
        pts=((pts-2)*2)+pts;
        n=n+1;
        s=pow(2,n-1);
        pts=2*s+1;
       // printf ("number of points %i \n",2*s+1);
        e=abs((l-l_old)/l);
        it=it+1;
    }  // set the desired accuracy here <<<<<<<<<<<<<
    
    
    
    return l;
}

