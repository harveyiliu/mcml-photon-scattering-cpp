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
      3. CORNEA (1060-nm)
      4. EYE_ANTERIOR (1060-nm)
    N = number of photons to be used in the Monte-Carlo simulation
    beamType = predefined beam type,
                    1. FLAT
                    2. GAUSSIAN
    P = incident beam power/energy (W or J)
    W = beam width (mm)
  Run command example:
    ./mcml BARE_DERMIS 1000 FLAT 1.0 10.0 &
 ****/

#include <iostream>
#include <ctime>
#include "mcml_model.h"
#include "mcml_conv.h"

using namespace std;

/***********************************************************
 
 ****/
int main(int argc, char * argv[]) 
{
  // Check the number of parameters
  if (argc < 4) {
    // Tell the user how to run the program
    cerr << "Usage: " << argv[0] << " MODEL NAME" << " PHOTON NUMBER" 
        << "CONV NAME" << endl;
    return 1;
  }
  string modelName(argv[1]);
  string numPhotonsSetStr(argv[2]);
  string beamType(argv[3]);      // beam type
  string PStr(argv[4]);         // beam power/energy [W or J]
  string WStr(argv[5]);         // beam width [mm]
    
  string::size_type sz;   // alias of size_t
  
  long numPhotonsSet = stol (numPhotonsSetStr, &sz);
  double P = stod (PStr, &sz);
  double W = stod (WStr, &sz);
  double R = 0.1*W/2;                 // beam radius [cm]

  clock_t start;
  double runTime;

  start = clock();

  MCMLModel model;
  model.SelectMCMLModel(modelName);
  model.DoOneRun(numPhotonsSet);
  model.SumScaleResult();
  
  cout.precision(2);
  cout << "\n" << "numPhotons: " << model.numPhotons << "\n" << fixed;
  cout << "Specular reflection Rsp (%): " << 100*model.Rsp << "\n";
  cout << "Diffused reflection Rd (%): " << 100*model.Rd << "\n";
  cout << "Absorption A (%): " << 100*model.A << "\n";
  cout << "Transmission Tt (%): " << 100*model.Tt << endl;

  MCMLConv conv;
  conv.SelectMCMLConv(model, beamType, P, R);
  conv.RunConv();
  double halfMaxDepth = conv.CenterHalfMaxDepth();
  double halfMaxWidth = conv.SurfaceHalfMaxWidth();

  // duration [sec.]
  runTime = (clock() - start)/(double) CLOCKS_PER_SEC;

  cout.precision(3);
  cout << "\n" << "Center half max depth (mm): " << 10*halfMaxDepth << "\n";
  cout << "Surface half max width (mm): " << 10*halfMaxWidth << "\n";
  cout.precision(1);
  cout << "Total run time (sec.): " << runTime << endl;

  conv.FreeMCMLConv();
  return(0);
}
