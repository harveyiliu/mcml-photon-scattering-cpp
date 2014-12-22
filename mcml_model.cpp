/***********************************************************
 *	class setup program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/

#include "mcml_model.h"

#define STANDARDTEST 0
  /* testing program using fixed rnd seed. */

#define PARTIALREFLECTION 0     
  /* 1=split photon, 0=statistical reflection. */

#define COSZERO (1.0-1.0e-12)	
  /* cosine of about 1e-6 rad. */

#define COS90D  1.0e-6		
  /* cosine of about 1.57 - 1e-6 rad. */

double RFresnel(double n1, double n2, double ca1, double * ca2Ptr);
double SpinTheta(double g);
float ran3(int * idum);
double RandomNum();

void Medium::SelectMedium (Medium::MediumName mediumName) {
  switch (mediumName) {
    case Medium::AIR:
      n = 1.0;
      mua = 0.0;
      mus = 0.0;
      g = 1.0;
      break;
    case Medium::DERMIS:    // 800-nm wavelngth
      n = 1.4;
      mua = 0.26;
      mus = 137.0;
      g = 0.90;
      break;
    case Medium::TYPE_II_EPIDERMIS:   // 800-nm
      n = 1.3;
      mua = 5.0;
      mus = 200.0;
      g = 0.70;
      break;    
    case Medium::CORNEA:    // 1060-nm
      n = 1.376;
      mua = 0.157;    // from Bovine NIR paper 2011
      mus = 1.064;
      g = 0.90;
      break;
    case Medium::AQUEOUS_HUMOR:   // 1060-nm
      n = 1.336;
      mua = 0.78;     // from Optical properties of ocular tissues
      mus = 0.60;
      g = 0.99;
      break;
    case Medium::LENS:   // 1060-nm
      n = 1.406;
      mua = 0.206;     // from Bovine NIR paper 2011
      mus = 1.131;
      g = 0.90;
      break;
    case Medium::VITREOUS_HUMOR:   // 1060-nm
      n = 1.337;
      mua = 0.48;     // from Optical properties of ocular tissues
      mus = 0.25;
      g = 0.99;
      break;
    case Medium::RETINA:   // 1060-nm
      n = 1.358;
      mua = 2.745;     // from Bovine NIR paper 2011
      mus = 71.50;
      g = 0.70;
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
    case LayerStruct::BARE_DERMIS:
      numLayers = 1;
      layer = new Medium [numLayers+2];
      layer[0].SelectMedium(Medium::AIR);
      layer[1].SelectMedium(Medium::DERMIS);
      layer[2].SelectMedium(Medium::DERMIS);
      layerThickness = new double [numLayers];
      layerThickness[0] = 0.3;
      break; 
    case LayerStruct::TYPE_II_SKIN:
      numLayers = 2;  // number of layers
      layer = new Medium [numLayers+2];
      layer[0].SelectMedium(Medium::AIR);     // incidence medium
      layer[1].SelectMedium(Medium::TYPE_II_EPIDERMIS);
      layer[2].SelectMedium(Medium::DERMIS);
      layer[3].SelectMedium(Medium::DERMIS);  //exit medium
      layerThickness = new double [numLayers];
      layerThickness[0] = 0.006;
      layerThickness[1] = 0.3;
      break;
    case LayerStruct::CORNEA:
      numLayers = 1;
      layer = new Medium [numLayers+2];
      layer[0].SelectMedium(Medium::AIR);
      layer[1].SelectMedium(Medium::CORNEA);
      layer[2].SelectMedium(Medium::AQUEOUS_HUMOR);
      layerThickness = new double [numLayers];
      layerThickness[0] = 0.0449;
      break;
    case LayerStruct::EYE_ANTERIOR:
      numLayers = 3;
      layer = new Medium [numLayers+2];
      layer[0].SelectMedium(Medium::AIR);
      layer[1].SelectMedium(Medium::CORNEA);
      layer[2].SelectMedium(Medium::AQUEOUS_HUMOR);
      layer[3].SelectMedium(Medium::LENS);
      layer[4].SelectMedium(Medium::VITREOUS_HUMOR);
      layerThickness = new double [numLayers];
      layerThickness[0] = 0.0449;
      layerThickness[1] = 0.2794;
      layerThickness[2] = 0.4979;
      break;
    default:
      numLayers = 1;
      layer = new Medium [numLayers+2];
      layer[0].SelectMedium(Medium::AIR);
      layer[1].SelectMedium(Medium::DERMIS);
      layer[2].SelectMedium(Medium::DERMIS);
      layerThickness = new double [numLayers];
      layerThickness[0] = 0.3;     
  }
  layerZ[0] = new double [numLayers+2];   // include the incidence and exit
  layerZ[1] = new double [numLayers+2];
  cosCrit[0] = new double [numLayers+2];
  cosCrit[1] = new double [numLayers+2];
  int z = 0;
  double n1, n2;
  for (int i=1; i <= numLayers; i++) {
    layerZ[0][i] = z;
    layerZ[1][i] = z + layerThickness[i-1];
    z = layerZ[1][i];
    // calculate the critical angle cosines for each layer
    // crticial angle at top interface of the current layer
    n1 = layer[i].n;
    n2 = layer[i-1].n;
    if (n1 > n2)
      cosCrit[0][i] = sqrt(1.0 - n2*n2/(n1*n1));
    else
      cosCrit[0][i] = 0.0;
    // crticial angle at bottom interface of the current layer
    n2 = layer[i+1].n;
    if (n1 > n2)
      cosCrit[1][i] = sqrt(1.0 - n2*n2/(n1*n1));
    else
      cosCrit[1][i] = 0.0;
  }
}


void LayerStruct::FreeLayerStruct () {
  if (layer != nullptr) {
    delete[] layer;
    layer = nullptr;
  }
  if (layerThickness != nullptr) {
    delete[] layerThickness;
    layerThickness = nullptr;
  }
  if (layerZ[0] != nullptr) {
    delete[] layerZ[0];
    layerZ[0] = nullptr;
    delete[] layerZ[1];
    layerZ[1] = nullptr;
  }
  if (cosCrit[0] != nullptr) {
    delete[] cosCrit[0];
    cosCrit[0] = nullptr;
    delete[] cosCrit[1];
    cosCrit[1] = nullptr;
  }
}


double LayerStruct::CalcRSpecular () {
/* direct reflections from the 1st and 2nd layers. */

  double temp, r1, r2;

  temp = (layer[0].n - layer[1].n)/(layer[0].n + layer[1].n);
  r1 = temp*temp;
  if ((layer[1].mua == 0.0) && (layer[1].mus == 0.0)) {
    // glass layer.
    temp = (layer[1].n - layer[2].n)/(layer[1].n + layer[2].n);
    r2 = temp*temp;
    r1 = r1 + (1 - r1)*(1 - r1)*r2/(1 - r1*r2);
  } 
  return r1;
}


void ModelInput::SelectModelInput (ModelInput::ModelInputName modelInputName) {
  
  this->FreeModelInput();
  switch (modelInputName) {
    case LayerStruct::BARE_DERMIS:
      dz = 100e-4;
      dr = 100e-4;
      nz = 30;
      nr = 50;
      na = 10;
      layerObj.SelectLayerStruct(LayerStruct::BARE_DERMIS);
      break;   
    case LayerStruct::TYPE_II_SKIN:
      dz = 20e-4;
      dr = 20e-4;
      nz = 10;
      nr = 50;
      na = 10;
      layerObj.SelectLayerStruct(LayerStruct::TYPE_II_SKIN);
      break;
    case LayerStruct::CORNEA:
      dz = 10e-4;
      dr = 10e-4;
      nz = 100;
      nr = 50;
      na = 10;
      layerObj.SelectLayerStruct(LayerStruct::CORNEA);
      break;
    case LayerStruct::EYE_ANTERIOR:
      dz = 20e-4;
      dr = 20e-4;
      nz = 500;
      nr = 250;
      na = 10;
      layerObj.SelectLayerStruct(LayerStruct::EYE_ANTERIOR);
      break;
    default:
      dz = 100e-4;
      dr = 100e-4;
      nz = 30;
      nr = 50;
      na = 10;
      layerObj.SelectLayerStruct(LayerStruct::BARE_DERMIS);
  }
  Wth = WEIGHT;
  da = 0.5*PI/na; 
}


void ModelInput::FreeModelInput () {
  layerObj.FreeLayerStruct();  
}



void MCMLModel::SelectMCMLModel (std::string modelName) {
  
  this->FreeMCMLModel();

  if (modelName.compare("BARE_DERMIS") == 0)
    this->SelectModelInput (ModelInput::BARE_DERMIS);
  else if (modelName.compare("TYPE_II_SKIN") == 0)
    this->SelectModelInput (ModelInput::TYPE_II_SKIN);
  else if (modelName.compare("CORNEA") == 0)
    this->SelectModelInput (ModelInput::CORNEA);
  else if (modelName.compare("EYE_ANTERIOR") == 0)
    this->SelectModelInput (ModelInput::EYE_ANTERIOR);
  else
    this->SelectModelInput (ModelInput::BARE_DERMIS);    

  numPhotons = 0;
  Rsp = layerObj.CalcRSpecular();
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
  if (Rd_ra != nullptr) {
    for (int i = 0; i < nr; i++) 
      delete[] Rd_ra[i];
  }
  Rd_ra = nullptr;
  if (Rd_r != nullptr) 
    delete[] Rd_r;
  Rd_r = nullptr;
  if (Rd_a != nullptr) 
    delete[] Rd_a;
  Rd_a = nullptr;
  if (A_rz != nullptr) {
    for (int i = 0; i < nr; i++) 
      delete[] A_rz[i];
  }
  A_rz = nullptr;
  if (A_z != nullptr) 
    delete[] A_z;
  A_z = nullptr;
  if (A_l != nullptr) 
    delete[] A_l;
  A_l = nullptr;
  if (Tt_ra != nullptr) {
    for (int i = 0; i < nr; i++) 
      delete[] Tt_ra[i];
  }
  Tt_ra = nullptr;
  if (Tt_r != nullptr) 
    delete[] Tt_r;
  Tt_r = nullptr;
  if (Tt_a != nullptr) 
    delete[] Tt_a;
  Tt_a = nullptr;  
  this->FreeModelInput();  
}



void MCMLModel::DoOneRun (long numPhotonsSet) {
  Photon photon;
  for (long i = 0; i < numPhotonsSet; i++) {
    photon.RunOnePhoton(this);
  }
}


void MCMLModel::SumScaleResult () {
  // Get 1D & 0D results.
  Sum2DRd();
  Sum2DA();
  Sum2DTt();  
  ScaleRdTt();
  ScaleA();
}


void MCMLModel::Sum2DRd () {
  // Get 1D array elements by summing the 2D array elements.
  short ir,ia;
  double sum; 

  for (ir = 0; ir < nr; ir++) {
    sum = 0.0;
    for (ia = 0; ia < na; ia++)
      sum += Rd_ra[ir][ia];
    Rd_r[ir] = sum;
  }
  
  for (ia = 0; ia < na; ia++) {
    sum = 0.0;
    for (ir = 0; ir < nr; ir++)
      sum += Rd_ra[ir][ia];
    Rd_a[ia] = sum;
  }
  
  sum = 0.0;
  for (ir = 0; ir < nr; ir++)
    sum += Rd_r[ir];
  Rd = sum;
}


short MCMLModel::IzToLayer (short iz) {
/* Return the index to the layer according to the index
  to the grid line system in z direction (Iz).
  Use the center of box. */
  
  short i;

  i = 1;     	// index to layer.
  while ((iz+0.5)*dz >= layerObj.layerZ[1][i] && i < layerObj.numLayers)
    i++;

  return i;
}


void MCMLModel::Sum2DA () {
// Get 1D array elements by summing the 2D array elements.
  short iz,ir;
  double sum;

  for (iz = 0; iz < nz; iz++) {
    sum = 0.0;
    for (ir = 0; ir < nr; ir++)
      sum += A_rz[ir][iz];
    A_z[iz] = sum;
  } 
  sum = 0.0;
  for (iz = 0; iz < nz; iz++) {
    sum += A_z[iz];
    A_l[IzToLayer(iz)] += A_z[iz];
  }
  A = sum;
}


void MCMLModel::Sum2DTt () {
// Get 1D array elements by summing the 2D array elements.
  short ir,ia;
  double sum;
  
  for (ir = 0; ir < nr; ir++) {
    sum = 0.0;
    for (ia = 0; ia < na; ia++)
      sum += Tt_ra[ir][ia];
    Tt_r[ir] = sum;
  }
  
  for (ia = 0; ia < na; ia++) {
    sum = 0.0;
    for (ir = 0; ir < nr; ir++)
      sum += Tt_ra[ir][ia];
    Tt_a[ia] = sum;
  }
  
  sum = 0.0;
  for (ir = 0; ir < nr; ir++)
    sum += Tt_r[ir];
  Tt = sum;  
}

void MCMLModel::ScaleRdTt () {
/* Scale Rd and Tt properly.
  "a" stands for angle alpha.
  Scale Rd(r,a) and Tt(r,a) by
  (area perpendicular to photon direction)
  x(solid angle)x(No. of photons).
  or
  [2*PI*r*dr*cos(a)]x[2*PI*sin(a)*da]x[No. of photons]
  or
  [2*PI*PI*dr*da*r*sin(2a)]x[No. of photons]
  Scale Rd(r) and Tt(r) by
  (area on the surface)x(No. of photons).
  Scale Rd(a) and Tt(a) by
  (solid angle)x(No. of photons). */

  short ir,ia;
  double scale1, scale2;
  
  scale1 = 4.0*PI*PI*dr*sin(da/2)*dr*numPhotons;

  // The factor (ir+0.5)*sin(2a) to be added.

  for (ir = 0; ir < nr; ir++) {  
    for (ia = 0; ia < na; ia++) {
      scale2 = 1.0/((ir+0.5)*sin(2.0*(ia+0.5)*da)*scale1);
      Rd_ra[ir][ia] *= scale2;
      Tt_ra[ir][ia] *= scale2;
    }
  }
  
  scale1 = 2.0*PI*dr*dr*numPhotons;  
  // area is 2*PI*[(ir+0.5)*dr]*dr. 
  // ir+0.5 to be added.

  for (ir = 0; ir < nr; ir++) {
    scale2 = 1.0/((ir+0.5)*scale1);
    Rd_r[ir] *= scale2;
    Tt_r[ir] *= scale2;
  }
  
  scale1  = 2.0*PI*da*numPhotons;
  // solid angle is 2*PI*sin(a)*da. sin(a) to be added.

  for (ia = 0; ia < na; ia++) {
    scale2 = 1.0/(sin((ia+0.5)*da)*scale1);
    Rd_a[ia] *= scale2;
    Tt_a[ia] *= scale2;
  }
  
  scale2 = 1.0/numPhotons;
  Rd *= scale2;
  Tt *= scale2;
}


void MCMLModel::ScaleA () {
/* Scale absorption arrays properly.
  Scale A_rz */

  short iz,ir;
  short il;
  double scale1;
        
  scale1 = 2.0*PI*dr*dr*dz*numPhotons;	
  //volume is 2*pi*(ir+0.5)*dr*dr*dz. 
  // ir+0.5 to be added.
  for (iz = 0; iz < nz; iz++)
    for (ir = 0; ir < nr; ir++)
      A_rz[ir][iz] /= (ir+0.5)*scale1;
  
  // Scale A_z.
  scale1 = 1.0/(dz*numPhotons);
  for (iz = 0; iz < nz; iz++)
    A_z[iz] *= scale1;
  
  // Scale A_l. Avoid int/int.
  scale1 = 1.0/numPhotons;	
  for (il = 0; il < layerObj.numLayers+2; il++)
    A_l[il] *= scale1;
  
  A *=scale1;
}



void Photon::Reset (MCMLModel * model) {
  x = 0.0;
  y = 0.0;
  z = 0.0;
  ux = 0.0;
  uy = 0.0;
  uz = 1.0;
  w = 1.0 - model->Rsp;
  dead = false;
  layer = 1;
  s = 0;
  sleft = 0;
  
  // take care of the case when the first layer is glass
  if ((model->layerObj.layer[1].mua == 0.0) && 
      (model->layerObj.layer[1].mus == 0.0)) {
    layer = 2;      // skip to next layer
    z = model->layerObj.layerZ[0][2];  // use z0 from the next layer
  }
}


void Photon::RunOnePhoton (MCMLModel * model) {
// run a single photon scattering till completion
  Reset(model);
  while (!dead) {
    HopDropSpin(model);
  }
  model->numPhotons += 1;
}

void Photon::HopDropSpin(MCMLModel * model) {
  if ((model->layerObj.layer[layer].mua == 0.0) &&
      (model->layerObj.layer[layer].mus == 0.0))
    // glass layer
    HopInGlass(model);
  else
    HopDropSpinInTissue(model);
  if ((w < model->Wth) && dead)
    Roulette();
}


void Photon::HopInGlass(MCMLModel * model) {
/* Move the photon packet in glass layer.
  Horizontal photons are killed because they will
  never interact with tissue again. */

  if (uz == 0.0)
  // horizontal photon in glass is killed
    dead = true;
  else {
    StepSizeInGlass(model);
    Hop();
    CrossOrNot(model);
  }
}


void Photon::HopDropSpinInTissue(MCMLModel * model) {
/* Set a step size, move the photon, drop some weight, 
  choose a new photon direction for propagation.
  When a step size is long enough for the photon to 
  hit an interface, this step is divided into two steps. 
  First, move the photon to the boundary free of 
  absorption or scattering, then decide whether the 
  photon is reflected or transmitted.
  Then move the photon in the current or transmission 
  medium with the unfinished stepsize to interaction 
  site.  If the unfinished stepsize is still too long, 
  repeat the above process. */

  StepSizeInTissue(model);
  if (HitBoundary(model)) {
    Hop();      // move to boundary plane
    CrossOrNot(model);
  } else {
    Hop();
    Drop(model);
    Spin(model->layerObj.layer[layer].g);
  }
}


void Photon::StepSizeInGlass(MCMLModel * model) {
/* If uz != 0, return the photon step size in glass, 
  Otherwise, return 0.
  The step size is the distance between the current 
  position and the boundary in the photon direction.
  Make sure uz !=0 before calling this function. */

  // Stepsize to the boundary.
  double dl_b;
  if (uz > 0.0)
    dl_b = (model->layerObj.layerZ[1][layer] - z)/uz;
  else if (uz < 0.0)
    dl_b = (model->layerObj.layerZ[0][layer] - z)/uz;
  else
    dl_b = 0.0;
  
  s = dl_b;

}

void Photon::StepSizeInTissue(MCMLModel * model) {
/* Pick a step size for a photon packet when it is in 
  tissue.
  If the member sleft is zero, make a new step size 
  with: -log(rnd)/(mua+mus).
  Otherwise, pick up the leftover in sleft.
  Layer is the index to layer.
  In_Ptr is the input parameters. */

  double mua, mus;

  mua = model->layerObj.layer[layer].mua;
  mus = model->layerObj.layer[layer].mus;

  if (sleft == 0.0)       // make a new step
    s = -log(RandomNum())/(mua + mus);
  else {                       // take the leftover
    s = sleft/(mua + mus);
    sleft = 0.0;
  }
}


void Photon::Hop() {
// Move the photon s away in the current layer of medium.
  x += s*ux;
  y += s*uy;
  z += s*uz;
}


void Photon::CrossOrNot(MCMLModel * model) {
  if (uz < 0.0)
    CrossUpOrNot(model);
  else
    CrossDnOrNot(model);
}


void Photon::CrossUpOrNot(MCMLModel * model) {
/* Decide whether the photon will be transmitted or 
  reflected on the upper boundary (uz<0) of the current 
  layer.
  If "layer" is the first layer, the photon packet will 
  be partially transmitted and partially reflected if 
  PARTIALREFLECTION is set to 1,
  or the photon packet will be either transmitted or 
  reflected determined statistically if PARTIALREFLECTION 
  is set to 0.
  Record the transmitted photon weight as reflection.
  If the "layer" is not the first layer and the photon 
  packet is transmitted, move the photon to "layer-1".
  Update the photon parmameters. */

  double uz0 = uz;
  double r = 0.0;     // reflectance
  double uz1;	// cosines of transmission alpha. always positive
  short layer0 = layer;
  double ni = model->layerObj.layer[layer0].n;
  double nt = model->layerObj.layer[layer0-1].n;
  
  
  // Get r.
  if (-uz0 <= model->layerObj.cosCrit[0][layer0])
    r = 1.0;		      // total internal reflection
  else
    r = RFresnel(ni, nt, -uz0, &uz1);
  
#if PARTIALREFLECTION
  if ((layer0 == 1) && (r < 1.0)) {      // partially transmitted
    uz = -uz1;          // transmitted photon
    RecordR(model, r);
    uz = -uz0;           // reflected photon
  } else if (RandomNum() > r) {    // transmitted to layer-1
    layer -= 1;
    ux *= ni/nt;
    uy *= ni/nt;
    uz = -uz1;
  } else			      		// reflected
    uz = -uz0;
#else
  if (RandomNum() > r) {       // transmitted to layer-1
    if (layer0 == 1) {
      uz = -uz1;
      RecordR(model, 0.0);
      dead = true;
    } else {
      layer -= 1;
      ux *= ni/nt;
      uy *= ni/nt;
      uz = -uz1;
    }
  } else 						// reflected
    uz = -uz0;
#endif
}


void Photon::CrossDnOrNot(MCMLModel * model) {
/* Decide whether the photon will be transmitted  or be 
  reflected on the bottom boundary (uz>0) of the current 
  layer.
  If the photon is transmitted, move the photon to 
  "layer+1". If "layer" is the last layer, record the 
  transmitted weight as transmittance. See comments for 
  CrossUpOrNot.
  Update the photon parmameters. */
        
  double uz0 = uz;        // z directional cosine
  double r = 0.0;             // reflectance
  double uz1;	          // cosines of transmission alpha. always positive
  short layer0 = layer;
  double ni = model->layerObj.layer[layer0].n;
  double nt = model->layerObj.layer[layer+1].n;

  
  // Get r
  if (uz0 <= model->layerObj.cosCrit[1][layer]) 
    r = 1.0;		// total internal reflection
  else
    r = RFresnel(ni, nt, uz0, &uz1);
  
#if PARTIALREFLECTION
  if ((layer0 == model->layerObj.numLayers) && (r < 1.0)) {
    uz = uz1;
    RecordT(model, r);
    uz = -uz0;
  } else if (RandomNum() > r) {     // transmitted to layer+1
    layer += 1;
    ux *= ni/nt;
    uy *= ni/nt;
    uz = uz1;
  } else 						// reflected
    uz = -uz0;
#else
  if (RandomNum() > r) {	      // transmitted to layer+1
    if (layer0 == model->layerObj.numLayers) {
      uz = uz1;
      RecordT(model, 0.0);
      dead = true;
    } else {
      layer += 1;
      ux *= ni/nt;
      uy *= ni/nt;
      uz = uz1;
    }
  } else 						// reflected
    uz = -uz0;
#endif
}



bool Photon::HitBoundary(MCMLModel * model) {
/* Check if the step will hit the boundary.
  Return 1 if hit boundary.
  Return 0 otherwise.
  If the projected step hits the boundary, the members
  s and sleft of Photon_Ptr are updated. */

  double dl_b;    // length to boundary
  double mut;
  bool hit;
 
// Distance to the boundary.
  if (uz > 0.0)
    dl_b = (model->layerObj.layerZ[1][layer] - z)/uz;	    // dl_b>0
  else if (uz < 0.0)
    dl_b = (model->layerObj.layerZ[0][layer] - z)/uz;    // dl_b>0
  
  if ((uz != 0.0) && (s > dl_b)) {
    // not horizontal & crossing
    mut = model->layerObj.layer[layer].mua + model->layerObj.layer[layer].mus;
    sleft = (s - dl_b)*mut;
    s = dl_b;
    hit = true;
  } else
    hit = false;
  
  return hit;
}


void Photon::Drop(MCMLModel * model) {
/* Drop photon weight inside the tissue (not glass).
  The photon is assumed not dead.
  The weight drop is dw = w*mua/(mua+mus).
  The dropped weight is assigned to the absorption array 
  elements. */
  
  double dwa;		  // absorbed weight.
  short  iz, ir;	    // index to z & r.
  double mua, mus;
 
  // compute array indices
  iz = (short)(z/model->dz);
  if (iz > (model->nz - 1))
    iz = model->nz - 1;
  
  ir = (short)(sqrt(x*x + y*y)/model->dr);
  if (ir > (model->nr - 1))
    ir = model->nr - 1;
  
  // update photon weight.
  mua = model->layerObj.layer[layer].mua;
  mus = model->layerObj.layer[layer].mus;
  dwa = w * mua/(mua+mus);
  w -= dwa;
  
  // assign dwa to the absorption array element.
  model->A_rz[ir][iz] += dwa;
}


void Photon::Spin(double g) {
/* Choose a new direction for photon propagation by 
  sampling the polar deflection angle theta and the 
  azimuthal angle psi.
  Note:
  theta: 0 - pi so sin(theta) is always positive 
  feel free to use sqrt() for cos(theta). 
  psi:   0 - 2pi 
  for 0-pi  sin(psi) is + 
  for pi-2pi sin(psi) is - */
  
  double cost, sint;	/* cosine and sine of the */
						/* polar deflection angle theta. */
  double cosp, sinp;	/* cosine and sine of the */
						/* azimuthal angle psi. */
  double ux0 = ux;
  double uy0 = uy;
  double uz0 = uz;
  double psi;
  
  cost = SpinTheta(g);
  sint = sqrt(1.0 - cost*cost);
  // sqrt() is faster than sin().

  psi = 2.0*PI*RandomNum();       // spin psi 0-2pi
  cosp = cos(psi);
  if (psi < PI)
    sinp = sqrt(1.0 - cosp*cosp);
    // sqrt() is faster than sin().
  else
    sinp = -sqrt(1.0 - cosp*cosp);	
  
  if (fabs(uz0) > COSZERO) {       // normal incident.
    ux = sint*cosp;
    uy = sint*sinp;
    uz = cost*SIGN(uz0);      // SIGN() is faster than division.
  }	else {       // regular incident.
    double temp = sqrt(1.0 - uz0*uz0);
    ux = sint*(ux0*uz0*cosp - uy0*sinp)/temp + ux0*cost;
    uy = sint*(uy0*uz0*cosp + ux0*sinp)/temp + uy0*cost;
    uz = -sint*cosp*temp + uz0*cost;
  }
}


void Photon::RecordR(MCMLModel * model, double refl) {
/* Record the photon weight exiting the first layer(uz<0), 
  no matter whether the layer is glass or not, to the 
  reflection array.
  Update the photon weight as well. */
  
  short  ir, ia;	    // index to r & angle.

  ir = (short)(sqrt(x*x + y*y)/model->dr);
  if (ir > (model->nr - 1))
    ir = model->nr - 1;
  
  ia = (short)(acos(-uz)/model->da); 
  if (ia > (model->na - 1))
    ia = model->na - 1;
  
  // assign photon to the reflection array element.
  model->Rd_ra[ir][ia] += w*(1.0 - refl);
  w *= refl;
}


void Photon::RecordT(MCMLModel * model, double refl) {
/* Record the photon weight exiting the last layer(uz>0), 
  no matter whether the layer is glass or not, to the 
  transmittance array.
  Update the photon weight as well. */

  short  ir, ia;	    // index to r & angle.
 
  ir = (short)(sqrt(x*x + y*y)/model->dr);
  if (ir > (model->nr - 1))
    ir = model->nr - 1;
  
  ia = (short)(acos(uz)/model->da);
  if (ia > (model->na - 1))
    ia = model->na - 1;
  
  // assign photon to the transmittance array element.
  model->Tt_ra[ir][ia] += w*(1.0 - refl); 
  w *= refl;
}


void Photon::Roulette() {
/* The photon weight is small, and the photon packet tries 
  to survive a roulette. */

  if (w == 0.0)	
    dead = true;
  else if (RandomNum() < CHANCE)    // survived the roulette.
    w /= CHANCE;
  else 
    dead = true;
}



double RFresnel(double n1, double n2, double ca1, double * ca2Ptr) {
/* Compute the Fresnel reflectance.
  Make sure that the cosine of the incident angle a1
  is positive, and the case when the angle is greater 
  than the critical angle is ruled out.
  Avoid trigonometric function operations as much as
  possible, because they are computation-intensive.
    n1 - incident refractive index
    n2 - transmit refractive index
    ca1 - cosine of the incident angle. 0<a1<90 degrees.
    ca2Ptr - pointer to the cosine of the transmission angle. a2>0
*/
  
  double r;
  
  if (n1 == n2) {			  	// matched boundary
    *ca2Ptr = ca1;
    r = 0.0;
  } else if (ca1 > COSZERO) {     // normal incident
    *ca2Ptr = ca1;
    r = (n2-n1)/(n2+n1);
    r *= r;
  } else if (ca1 < COS90D) {      // very slant
    *ca2Ptr = 0.0;
    r = 1.0;
  } else {           // general	
    double sa1, sa2;
    // sine of the incident and transmission angles
    double ca2;

    sa1 = sqrt(1.0 - ca1*ca1);
    sa2 = n1*sa1/n2;
    if (sa2 >= 1.0) {
      // double check for total internal reflection
      *ca2Ptr = 0.0;
      r = 1.0;
    } else {
      double cap, cam;    /* cosines of the sum ap or
                            difference am of the two
                            angles. ap = a1+a2
                            am = a1 - a2 */
      double sap, sam;	  // sines.
     
      *ca2Ptr = ca2 = sqrt(1.0 - sa2*sa2);     
      cap = ca1*ca2 - sa1*sa2;     // c+ = cc - ss
      cam = ca1*ca2 + sa1*sa2;     // c- = cc + ss
      sap = sa1*ca2 + ca1*sa2;     // s+ = sc + cs
      sam = sa1*ca2 - ca1*sa2;     // s- = sc - cs
      r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam);
    }
  } 
  return r;
}




double SpinTheta(double g) {
/* Choose (sample) a new theta angle for photon propagation
  according to the anisotropy.
  If anisotropy g is 0, then
  cos(theta) = 2*rand-1.
  otherwise
  sample according to the Henyey-Greenstein function.
  Returns the cosine of the polar deflection angle theta. */

  double temp, cost;
 
  if (g == 0.0) 
    cost = 2*RandomNum() - 1;
  else {
    temp = (1 - g*g)/(1 - g + 2*g*RandomNum());
    cost = (1 + g*g - temp*temp)/(2*g);
    if (cost < -1)
      cost = -1.0;
    else if (cost > 1)
      cost = 1.0;
  }
  return cost;
}


/***********************************************************
 *	A random number generator from Numerical Recipes in C.
 ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

float ran3(int *idum)
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  
  if (*idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC



/***********************************************************
 *	Generate a random number between 0 and 1.  Take a 
 *	number as seed the first time entering the function.  
 *	The seed is limited to 1<<15.  
 *	We found that when idum is too large, ran3 may return 
 *	numbers beyond 0 and 1.
 ****/
double RandomNum() {
  static bool first_time=true;
  static int idum;	/* seed for ran3. */
  
  if(first_time) {
#if STANDARDTEST /* Use fixed seed to test the program. */
    idum = - 1;
#else
    idum = -(int)time(NULL)%(1<<15);
	  /* use 16-bit integer as the seed. */
#endif
    ran3(&idum);
    first_time = 0;
    idum = 1;
  }
  
  return( (double)ran3(&idum) );
  
  
}
