#ifndef FIGURES_H
#define FIGURES_H

#include "interp.hpp"
#include "models.hpp"
#include "data.hpp"

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

#include "rd.hpp"

using namespace std;

template <class Model>
class Figures_RD
{
private:
Relic_density<Model> relic_density;
double m_s,lambda_hs;
Data data;
public:

Figures_RD(Data _data)
{
//  relic_density=_rd;
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