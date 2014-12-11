/***********************************************************
 *	main program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/

#include <iostream>
#include "mcml_class.h"

using namespace std;

/***********************************************************
 
 ****/
int main(void) 
{
  MCMLModel skin;
  skin.SelectMCMLModel(ModelInput::TYPE_II_SKIN);
  cout << "dz: " << skin.dz << "\n";
  cout << "da: " << skin.da << endl;
  skin.FreeMCMLModel();
  return(0);
}
