/***********************************************************
 *	main program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/

#include <iostream>
#include "mcml_model.h"

using namespace std;

/***********************************************************
 
 ****/
int main(void) 
{
  MCMLModel model;
  model.SelectMCMLModel(ModelInput::TYPE_II_SKIN);

  model.DoOneRun(100);
  model.SumScaleResult();
  
  cout << "numPhotons: " << model.numPhotons << "\n";
  cout << "Rsp: " << model.Rsp << "\n";
  cout << "Rd: " << model.Rd << "\n";
  cout << "A: " << model.A << "\n";
  cout << "Tt: " << model.Tt << endl;

  model.FreeMCMLModel();
  return(0);
}
