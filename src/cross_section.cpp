#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

#include "cross_section.hpp"
#include "interp.hpp"

using namespace std;

double Cross_section::calc_cross_section(double T)
{
//// need to find when the integral converges, this can be done by testing the integrand and finding where it falls below a critical value
//double p;
//double integrand=1;
//double lower_cutoff=4*pow(m_s,2);
//double s;
//
//p=log10(lower_cutoff);
//
//while (integrand>0)
//{
//p=p+0.01;
//s=pow(10,p)+lower_cutoff;
////cout<< "s = " << s << endl;
//integrand=cs_integral(s,T);
////cout << "integrand = " << integrand << endl;
//}
////cout << "lower cs limit set = " << s << " lower cutoff = " << lower_cutoff<<  endl;
//
//double result;
//
//if (cs_integral(s=pow(10,p-0.1)+lower_cutoff,T)==0)
//{
//result = integrator_2(4*pow(m_s,2),s,T);
//}
//else
//{
//result = integrator(4*pow(m_s,2),s,T);
//}
//return result;



}
///  Cross section integral class functions //


double Cross_section::cs_integral(double s, double T)
{
double a= s*pow(s-4*pow(m_s,2),0.5)* sigma_v(s);
double b= (16 * T * pow(m_s,4));
return exp( log (a/b)+ log_bessel_k_asymtotic(1,pow(s,0.5)/T)- 2 * log_bessel_k_asymtotic(2,m_s/T) );
//double result_1=s*pow(s-4*pow(m_s,2),0.5)*boost::math::cyl_bessel_k(1,pow(s,0.5)/T)*sigma_v(s);
//double result_2=(16 * T * pow(m_s,4) * pow(boost::math::cyl_bessel_k(2,m_s/T),2));
}

double Cross_section::sigma_v_VV(double s)
{
double result=0;
double delta_W=1,delta_Z=0.5;
double x_W=pow(M_W,2)/s, x_Z=pow(M_Z,2)/s;
double v_W=pow(1-4*x_W,0.5) , v_Z=pow(1-4*x_Z,0.5);
if (pow(s,0.5)<2*M_W)
{
result=0;
}
else
{
result = pow(lambda_hs,2)*s*delta_W*v_W*Dh2(s)*(1-4*x_W+12*pow(x_W,2))/(8*Pi);
}

if (pow(s,0.5)<2*M_Z)
{
result=result+0;
}
else
{
result =result + pow(lambda_hs,2)*s*delta_Z*v_Z*Dh2(s)*(1-4*x_Z+12*pow(x_Z,2))/(8*Pi);
}
return result;

}


double Cross_section::sigma_v_ff(double s)
{
double result;
if (pow(s,0.5)<2*mt)
{
result=0;
}
else
{
// just do top quark for now
double X_t=3*(1+(1.5*log(pow(mt,2)/s)+float(9)/4)*float(4)*alpha_s/(float(3)*Pi));
double v_t=pow(1-4*pow(mt,2)/s,0.5);
result=pow(lambda_hs*mt,2) * X_t* pow(v_t,3)* Dh2(s)  /(4*Pi);
}

return result;
}

double Cross_section::sigma_v_hh(double s)
{
double result;
if (pow(s,0.5)<2*Mh)
{
result=0;
}
else
{
double v_s=pow(1-4*pow(m_s,2)/s,0.5);
double v_h=pow(1-4*pow(Mh,2)/s,0.5);

double tp=pow(m_s,2)+pow(Mh,2)-0.5*s*(1-v_s*v_h);
double tm=pow(m_s,2)+pow(Mh,2)-0.5*s*(1+v_s*v_h);

double a_R=1+3*pow(Mh,2)*(s-pow(Mh,2))*Dh2(s);
double a_I=3*pow(Mh,3)*pow(s,0.5)*gamma_h(Mh)*Dh2(s);

double prefactor=(pow(lambda_hs,2)/(16*Pi*pow(s,2)*v_s));

double a =  (pow(a_R,2)+pow(a_I,2))*s*v_s*v_h;
double b =  4 * lambda_hs * pow(v_0,2)* (a_R- lambda_hs*pow(v_0,2)/(s-2*pow(Mh,2))) * log(abs( (pow(m_s,2)-tp)/(pow(m_s,2)-tm)));
double c =  2 * pow(lambda_hs,2)*pow(v_0,4)*s*v_s*v_h/ ( (pow(m_s,2)-tm)*(pow(m_s,2)-tp));
result = prefactor*(a+b+c);
}

return result;

}


double Cross_section::sigma_v(double s)
{
double sigv;

if (pow(s,0.5)<1000 && pow(s,0.5)>90)
{
sigv= (2*pow(lambda_hs,2)*pow(v_0,2)/(pow(s,0.5)))*(abs(Dh2(s)))*gamma_h(pow(s,0.5));
}
else if (pow(s,0.5)>=1000)
{
sigv=sigma_v_ff(s)+sigma_v_VV(s)+sigma_v_hh(s);

}
else
{
sigv=0;
cout << "using low estimate for gamma" << endl;
cout<< "requested s value is = " << s << endl;
}
return sigv;
}
double Cross_section::Dh2(double s)
{
return 1/(pow((s-pow(Mh,2)),2)+(pow(Mh,2))*pow(gamma_h(Mh),2));
}


double Cross_section::gamma_h(double x)
{
std::vector<double> mh_input  ={90,100,110,120,130,140,150,180,200,250,300,400,500,600,800,1000};
std::vector<double> gamma_input  ={2.2e-3,2.46e-3,2.82e-3,3.47e-3,4.87e-3,8.12e-3,1.73e-2,6.31e-1,1.43,4.04,8.43,29.2,68,123,304,647};

double gamma;
if (x<1000 && x>90)
{

Poly_interp myfunc(mh_input,gamma_input,4);
gamma=myfunc.interp(x);
}
else{
gamma=0;
cout<< "gamma requested but outside of range at " << x << endl;
return 0;
}
return gamma;
}

// end of cross section integral functions //



//
//
//
//double Cross_section::integrator(double lower, double upper, double T)
//{
//    double fa,fb,l=0,h,l_old=0;
//    double e=1;
//    int it;
//    it=0;
//    double a=lower,b=upper;
//    fa=cs_integral(lower, T);
//    fb=cs_integral(upper, T);
// //   cout << "fa = " << fa << endl;
//  //  cout << "fb = " << fb << endl;
//    h=(b-a)/2;
//    l=(h/2)*(fa+fb);
//    //printf ("First estimate is %f \n",l);
//    //double midpoint=(a+b)/2;
//    int n;
//    n=1;
//    int pts;
//    pts=100;
//    int s;
//    double c,c_prev;
//    s=1;
//  
//    while(e>1e-3)
//    {
//        l_old=l;
//        
//       // printf ("Iteration number %i \n",it);
//
//        //l=l/h;
//        //h=(b-a)/(pts-1); // h is the length of each section in the integration
//        //l=l*h;
//        l=0;
//        c_prev=log10(a);
//        for (int m=1; m<2*s+1; m +=2)
//        {
//            c=m*(log10(b)-log10(a))/(2*s)+log10(a);
//            h=pow(10,c)-pow(10,c_prev);
//            l=l+h*cs_integral(pow(10,c),T);
//            c_prev=c;
//          
//            //cout<< "sampling at " << pow(10,c) << " obtained value = " << cs_integral(pow(10,c),T) << " h = " << h << endl;
//          
//
//        }
//      //  printf ("updated value of integral %f \n",l);
//       // cout << "updated value of integral = " << l << "prev  = " << l_old << endl;
//        pts=((pts-2)*2)+pts;
//        n=n+1;
//        s=pow(2,n-1);
//        pts=2*s+1;
//       // printf ("number of points %i \n",2*s+1);
//        e=abs((l-l_old)/l);
//        it=it+1;
//    }  // set the desired accuracy here <<<<<<<<<<<<<
//    
//    
//    
//    return l;
//}
//
//double Cross_section::integrator_2(double lower, double upper, double T)
//{
//    double fa,fb,l=0,h,l_old=0;
//    double e=1;
//    int it;
//    it=0;
//    double a=lower,b=upper;
//
//    //printf ("First estimate is %f \n",l);
//    //double midpoint=(a+b)/2;
//    int n;
//    n=1;
//    int pts;
//    pts=100;
//    int s; double f_upper;double b_old;
//    double c,c_prev;
//    s=1;
//    while (f_upper==0)
//    {
//    b_old=b;
//    b=(b-a)/2+a;
//    f_upper=cs_integral(b,T);
//    cout << "a , b = " << a << "  " << b << endl;
//    cout << "f _upper = " << f_upper << endl;
//    if (abs(a-b)==0)
//    {
//    f_upper=1;
//    }
//    }
//    if (f_upper==1)
//    {
//    return 0;
//    }
//  
//    f_upper=1;
//    double del=(b-a)/2;
//    while (f_upper>0)
//    {
//    b=b+del;
//    f_upper=cs_integral(b,T);
//    cout << "a , b = " << a << "  " << b << endl;
//    cout << "f _upper = " << f_upper << endl;
//    }
//  
//    fa=cs_integral(lower, T);
//    fb=cs_integral(upper, T);
//    //   cout << "fa = " << fa << endl;
//    //  cout << "fb = " << fb << endl;
//    h=(b-a)/2;
//    l=(h/2)*(fa+fb);
//  
//    while(e>1e-3)
//    {
//        l_old=l;
//        l=0;
//        c_prev=a;
//        for (int m=1; m<2*s+1; m +=2)
//        {
//            c=m*(b-a)/(2*s)+a;
//            h=c-c_prev;
//            l=l+h*cs_integral(c,T);
//            c_prev=c;
//        }
//        pts=((pts-2)*2)+pts;
//        n=n+1;
//        s=pow(2,n-1);
//        pts=2*s+1;
//       // printf ("number of points %i \n",2*s+1);
//        if (l>0)
//        {
//        e=abs((l-l_old)/l);
//        cout << "e = "  << e << endl;
//        }
//        else
//        {
//        e=1;
//        }
//        it=it+1;
//    }  // set the desired accuracy here <<<<<<<<<<<<<
//    cout << "integral_2 is = " << l << endl;
//    return l;
//}
//
