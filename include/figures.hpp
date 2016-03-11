#ifndef FIGURES_H
#define FIGURES_H

#include "interp.hpp"
#include "singletdm.hpp"
#include "data.hpp"

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

#include "rd.hpp"

using namespace std;

class Figures
{
private:
SingletDM cross_section;
Relic_density relic_density;
double m_s,lambda_hs;
Data data;
public:


Figures(SingletDM _cross_section, Relic_density _rd,Data _data)
{
  cross_section=_cross_section;
  relic_density=_rd;
 // m_s=cross_section.m_s;
 // lambda_hs=cross_section.lambda_hs;
 m_s=_data.M_s;
 lambda_hs=_data.Lam_hs;
 data=_data;
}


void plot_thermal_av();
void plot_sigma_v(double T);
void plot_Z();
void gamma_h();


  
};

#endif