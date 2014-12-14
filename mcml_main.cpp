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
  Photon photon(model.layerObj);
  photon.RunOnePhoton(&model);
  
  cout << "A_rz[0][0]: " << model.A_rz[0][0] << "\n";
  cout << "A_rz[nr-1][nz-1]: " << model.A_rz[model.nr-1][model.nz-1] << endl;
  model.FreeMCMLModel();
  return(0);
}
