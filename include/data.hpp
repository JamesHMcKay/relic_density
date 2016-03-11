// structure to hold all data and constants
// convention for masses is M_x (<upper case>_<lower case>) in all parts of the program


#ifndef DATA_H
#define DATA_H

#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <iostream>

#include <fstream>

#include "singletdm.hpp"

using namespace std;


struct Data
{
  public:
  
  
  // default SM parameters
  // masses
  double M_h=125.66;
  double M_z=91;
  double M_w=80;
  double M_t=173.15;
  //double Lam_h;
  //double mu_h;
  
  // constants
  double Pi=3.14159;
  double v0=246;
  double alpha_s=0.1184;
  double M_pl=1.2e19;
  
  // BSM parameters
  
  double Lam_hs=0.1;
  double Lam_s=0;
  
  double M_s=100;
  
  // constructor
  Data (){};
  
  Data(int argc, char* argv[]) {
  double param [10];
  std::string name [10]; int i=0;

  if (argc==1)
  {
  cout << "Please enter a data file in the format ./main input.txt, using default values " << endl;
  }
  else
  {
  std::ifstream input(argv[1]);
  std::string line;
  while(getline(input, line)) {
        if (!line.length() || line[0] == '#')
           continue;
        std::istringstream iss(line);
        iss>> name[i] >> param[i];
    
        i=i+1;
     }
  }
  for (int n=0;n<i+1;n++)
  {
  if (name[n]=="lambda_hs")
  {
  Lam_hs=param[n];
  }
  if (name[n]=="M_s")
  {
  M_s=param[n];
  }
  if (name[n]=="M_t")
  {
  M_t=param[n];
  }
  if (name[n]=="M_h")
  {
  M_h=param[n];
  }
  if (name[n]=="M_w")
  {
  M_w=param[n];
  }
  if (name[n]=="M_z")
  {
  M_z=param[n];
  }


  }
  }
  
  
  
  
  
  
  
  
  
  
  
};
#endif



