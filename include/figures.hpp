#ifndef FIGURES_H
#define FIGURES_H

#include "interp.hpp"

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

#include "cross_section.hpp"

using namespace std;

class Figures
{
private:
Cross_section cross_section;
Relic_density relic_density;
double m_s,lambda_hs;
public:


Figures(Cross_section _cross_section, Relic_density _rd)
{
  cross_section=_cross_section;
  relic_density=_rd;
  m_s=cross_section.m_s;
  lambda_hs=cross_section.lambda_hs;
}


void plot_thermal_av();
void plot_sigma_v(double T);
void plot_Z();
void gamma_h();


  
};

#endif