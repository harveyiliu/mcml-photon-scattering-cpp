/***********************************************************
 *	main program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
  Compile command:
    g++ -std=c++11 mcml_main.cpp mcml_model.cpp -o mcml
  Input arguments:
    modelName = predefined model name, current set of models are:
      1. BARE_DERMIS (800-nm) - default
      2. TYPE_II_SKIN (800-nm)
    N = number of photons to be used in the Monte-Carlo simulation
  Run command example:
    ./mcml TYPE_II_SKIN 1000 &
 ****/

#include <iostream>
#include "mcml_model.h"

using namespace std;

/***********************************************************
 
 ****/
int main(int argc, char * argv[]) 
{
  // Check the number of parameters
  if (argc < 3) {
    // Tell the user how to run the program
    cerr << "Usage: " << argv[0] << " MODEL NAME" << " PHOTON NUMBER" << endl;
    return 1;
  }
  string modelName(argv[1]);
  string numPhotonsSetStr(argv[2]);
  string::size_type sz;   // alias of size_t
  
  long numPhotonsSet = stol (numPhotonsSetStr, &sz);

  MCMLModel model;
  model.SelectMCMLModel(modelName);
  model.DoOneRun(numPhotonsSet);
  model.SumScaleResult();
  
  cout << "\n" << "numPhotons: " << model.numPhotons << "\n";
  cout << "Rsp: " << model.Rsp << "\n";
  cout << "Rd: " << model.Rd << "\n";
  cout << "A: " << model.A << "\n";
  cout << "Tt: " << model.Tt << endl;

  model.FreeMCMLModel();
  return(0);
}
