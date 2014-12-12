/***********************************************************
 *	class setup program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/


#include "mcml_class.h"


void Medium::SelectMedium (Medium::MediumName mediumName) {
  switch (mediumName) {
    case Medium::TYPE_II_EPIDERMIS:
      n = 1.3;
      mua = 5.0;
      mus = 200.0;
      g = 0.70;
      break;
    case Medium::DERMIS:
      n = 1.4;
      mua = 0.26;
      mus = 137.0;
      g = 0.90;
      break;
    default:
      n = 1.4;
      mua = 0.26;
      mus = 137.0;
      g = 0.90;
  }
}




void LayerStruct::SelectLayerStruct (LayerStruct::LayerName layerName) {
  
  this->FreeLayerStruct();    // free previously set array members if needed
  switch (layerName) {
    case LayerStruct::TYPE_II_SKIN:
      nIn = 1.0;      // incidence medium index
      nOut = 1.4;     // exit medium index
      numLayers = 2;  // number of layers
      layer = new Medium [numLayers];
      layer[0].SelectMedium(Medium::TYPE_II_EPIDERMIS);
      layer[1].SelectMedium(Medium::DERMIS);
      layerThickness = new double [numLayers];
      layerThickness[0] = 0.006;
      layerThickness[1] = 0.3;
      break;

    case LayerStruct::BARE_DERMIS:
      nIn = 1.0;
      nOut = 1.4;
      numLayers = 1;
      layer = new Medium [numLayers];
      layer[0].SelectMedium(Medium::DERMIS);
      layerThickness = new double [numLayers];
      layerThickness[0] = 0.3;
      break;

    default:
      nIn = 1.0;
      nOut = 1.4;
      numLayers = 1;
      layer = new Medium [numLayers];
      layer[0].SelectMedium(Medium::DERMIS);
      layerThickness = new double [numLayers];
      layerThickness[0] = 0.3;     
  }
  layerZ[0] = new double [numLayers];
  layerZ[1] = new double [numLayers];
  cosCrit[0] = new double [numLayers];
  cosCrit[1] = new double [numLayers];
  int z = 0;
  double n1, n2;
  for (int i=0; i < numLayers; i++) {
    layerZ[0][i] = z;
    layerZ[1][i] = z + layerThickness[i];
    z = layerZ[1][i];
    // calculate the critical angle cosines for each layer
    // crticial angle at top interface of the current layer
    n1 = layer[i].n;
    if (i == 0) {
      n2 = nIn;
    } else {
      n2 = layer[i-1].n;
    }
    if (n1 > n2) {
      cosCrit[0][i] = sqrt(1.0 - n2*n2/(n1*n1));
    } else {
      cosCrit[0][i] = 0.0;
    }
    // crticial angle at bottom interface of the current layer
    if ((i+1) == numLayers) {
      n2 = nOut;
    } else {
      n2 = layer[i+1].n;
    }
    if (n1 > n2) {
      cosCrit[1][i] = sqrt(1.0 - n2*n2/(n1*n1));
    } else {
      cosCrit[1][i] = 0.0;
    }
  }
}


void LayerStruct::FreeLayerStruct () {
  if (layer != NULL) {
    delete[] layer;
    layer = NULL;
  }
  if (layerThickness != NULL) {
    delete[] layerThickness;
    layerThickness = NULL;
  }
  if (layerZ[0] != NULL) {
    delete[] layerZ[0];
    layerZ[0] = NULL;
    delete[] layerZ[1];
    layerZ[1] = NULL;
  }
  if (cosCrit[0] != NULL) {
    delete[] cosCrit[0];
    cosCrit[0] = NULL;
    delete[] cosCrit[1];
    cosCrit[1] = NULL;
  }
}


Photon::Photon (LayerStruct layerObj, double rSpecular) {
  x = 0.0;
  y = 0.0;
  z = 0.0;
  ux = 0.0;
  uy = 0.0;
  uz = 1.0;
  w = 1.0 - rSpecular;
  dead = false;
  layer = 0;
  s = 0;
  sleft = 0;
  
  // take care of the case when the first layer is glass
  if ((layerObj.layer[0].mua == 0) && (layerObj.layer[0].mus == 0)) {
    layer = 1;      // skip to next layer
    z = layerObj.layerZ[0][1];  // use z0 from the next layer
  }
}


void ModelInput::SelectModelInput (ModelInput::ModelInputName modelInputName,
    long numPhotonsSet) {
  
  this->FreeModelInput();
  switch (modelInputName) {
    case LayerStruct::TYPE_II_SKIN:
      dz = 20e-4;
      dr = 20e-4;
      nz = 10;
      nr = 50;
      na = 10;
      layerObj.SelectLayerStruct(LayerStruct::TYPE_II_SKIN);
      break;
    
    case LayerStruct::BARE_DERMIS:
      dz = 100e-4;
      dr = 100e-4;
      nz = 30;
      nr = 50;
      na = 10;
      layerObj.SelectLayerStruct(LayerStruct::BARE_DERMIS);
      break;

    default:
      dz = 100e-4;
      dr = 100e-4;
      nz = 30;
      nr = 50;
      na = 10;
      layerObj.SelectLayerStruct(LayerStruct::BARE_DERMIS);
  }
  numPhotons = numPhotonsSet;
  Wth = WEIGHT;
  da = 0.5*PI/na; 
}


void ModelInput::FreeModelInput () {
  layerObj.FreeLayerStruct();  
}



void MCMLModel::SelectMCMLModel (ModelInput::ModelInputName modelInputName,
    long numPhotonsSet) {
  
  this->FreeMCMLModel();
  this->SelectModelInput (modelInputName, numPhotonsSet);

  Rsp = 0.0;
  Rd = 0.0;
  A = 0.0;
  Tt = 0.0;
   
  Rd_ra = new double * [nr];
  for (int i = 0; i < nr; i++) {
    Rd_ra[i] = new double [na];
    std::fill (Rd_ra[i], Rd_ra[i]+na-1, 0);   // initialize array with 0  
  }

  Rd_r = new double [nr];
  std::fill (Rd_r, Rd_r+nr-1, 0);
  Rd_a = new double [na];
  std::fill (Rd_a, Rd_a+na-1, 0);

  A_rz = new double * [nr];
  for (int i = 0; i < nr; i++) {
    A_rz[i] = new double [nz];
    std::fill (A_rz[i], A_rz[i]+nz-1, 0);  
  }

  A_z = new double [nz];
  std::fill (A_z, A_z+nz-1, 0);
  A_l = new double [2+layerObj.numLayers];
  std::fill (A_l, A_l+layerObj.numLayers+1, 0);

  Tt_ra = new double * [nr];
  for (int i = 0; i < nr; i++) {
    Tt_ra[i] = new double [na];
    std::fill (Tt_ra[i], Tt_ra[i]+na-1, 0);    
  }

  Tt_r = new double [nr];
  std::fill (Tt_r, Tt_r+nr-1, 0);
  Tt_a = new double [na];
  std::fill (Tt_a, Tt_a+na-1, 0);

}


void MCMLModel::FreeMCMLModel () {
  if (Rd_ra != NULL) {
    for (int i = 0; i < nr; i++) 
      delete[] Rd_ra[i];
  }
  Rd_ra = NULL;
  if (Rd_r != NULL) 
    delete[] Rd_r;
  Rd_r = NULL;
  if (Rd_a != NULL) 
    delete[] Rd_a;
  Rd_a = NULL;
  if (A_rz != NULL) {
    for (int i = 0; i < nr; i++) 
      delete[] A_rz[i];
  }
  A_rz = NULL;
  if (A_z != NULL) 
    delete[] A_z;
  A_z = NULL;
  if (A_l != NULL) 
    delete[] A_l;
  A_l = NULL;
  if (Tt_ra != NULL) {
    for (int i = 0; i < nr; i++) 
      delete[] Tt_ra[i];
  }
  Tt_ra = NULL;
  if (Tt_r != NULL) 
    delete[] Tt_r;
  Tt_r = NULL;
  if (Tt_a != NULL) 
    delete[] Tt_a;
  Tt_a = NULL;  
  this->FreeModelInput();  
}



