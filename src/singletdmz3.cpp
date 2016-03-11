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


double SingletDMZ3::cs_integral(double s, double T)
{
Bessel bessel;// class for asymptotic bessel function definitions
double a= s*pow(s-4*pow(data.M_s,2),0.5)* sigma_v(s);
double b= (16 * T * pow(data.M_s,4));
return exp( log (a/b)+ bessel.log_bessel_k_asymtotic(1,pow(s,0.5)/T)- 2 * bessel.log_bessel_k_asymtotic(2,data.M_s/T) );

////cout<< "calculating integrand at (s , T) = " << s << " " << T << endl;
//
//double result_1=s*pow(s-4*pow(m_s,2),0.5)*bessel_k_asymtotic(1,pow(s,0.5)/T)*sigma_v(s);
//double result_2=(16 * T * pow(m_s,4) * pow(bessel_k_asymtotic(2,m_s/T),2));
//
//return result_1/result_2;
}




// vector bosons
double SingletDMZ3::sigma_v_VV(double s)
{
double result=0;
double delta_W=1,delta_Z=0.5;
double x_W=pow(data.M_w,2)/s, x_Z=pow(data.M_z,2)/s;
double v_W=pow(1-4*x_W,0.5) , v_Z=pow(1-4*x_Z,0.5);

if (pow(s,0.5)<2*data.M_w)
{
result=0;
}
else
{
result = pow(data.Lam_hs,2)*s*delta_W*v_W*Dh2(s)*(1-4*x_W+12*pow(x_W,2))/(8*data.Pi);
}

if (pow(s,0.5)<2*data.M_z)
{
result=result+0;
}
else
{
result = result + pow(data.Lam_hs,2)*s*delta_Z*v_Z*Dh2(s)*(1-4*x_Z+12*pow(x_Z,2))/(8*data.Pi);
}

return result;
}

// fermions
double SingletDMZ3::sigma_v_ff(double s)
{
double result;
if (pow(s,0.5)<2*data.M_t)
{
result=0;
}
else
{
// just do top quark for now
double X_t=3*(1+(1.5*log(pow(data.M_t,2)/s)+float(9)/4)*float(4)*data.alpha_s/(float(3)*data.Pi));
double v_t=pow(1-4*pow(data.M_t,2)/s,0.5);
result=pow(data.Lam_hs*data.M_t,2) * X_t* pow(v_t,3)* Dh2(s)  /(4*data.Pi);
}
return result;
}



// Higgs
double SingletDMZ3::sigma_v_hh(double s)
{
double result;
double m_s=data.M_s,Mh=data.M_h;
double lambda_hs=data.Lam_hs,v_0=data.v0;
if (pow(s,0.5)<2*data.M_h)
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

double prefactor=(pow(lambda_hs,2)/(16*data.Pi*pow(s,2)*v_s));

double a =  (pow(a_R,2)+pow(a_I,2))*s*v_s*v_h;
double b =  4 * lambda_hs * pow(v_0,2)* (a_R- lambda_hs*pow(v_0,2)/(s-2*pow(Mh,2))) * log(abs( (pow(m_s,2)-tp)/(pow(m_s,2)-tm)));
double c =  2 * pow(lambda_hs,2)*pow(v_0,4)*s*v_s*v_h/ ( (pow(m_s,2)-tm)*(pow(m_s,2)-tp));
result = prefactor*(a+b+c);
}

return result;

}



// total sigma_v
// still need to add invisible contribution to propogator when ms < mh/2 (see comment on end of page 2 of Cline et al. 2013)

double SingletDMZ3::sigma_v(double s)
{
double sigv;

if (pow(s,0.5)<1000 && pow(s,0.5)>90)
{
sigv = (2*pow(data.Lam_hs,2)*pow(data.v0,2)/(pow(s,0.5)))*(abs(Dh2(s)))*gamma_h(pow(s,0.5))+sigma_v_hh(s);
}
else if (pow(s,0.5)>=1000)
{
sigv=sigma_v_ff(s)+sigma_v_VV(s)+sigma_v_hh(s);
}
else
{
sigv=0;
cout<< "requested s value is too low = " << s << endl;
}

return sigv;
}


double SingletDMZ3::Dh2(double s)
{
return 1/(pow((s-pow(data.M_h,2)),2)+(pow(data.M_h,2))*pow(gamma_h(data.M_h),2));
}





double SingletDMZ3::gamma_h(double x)
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



