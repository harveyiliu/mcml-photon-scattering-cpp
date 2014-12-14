/***********************************************************
 *	class setup program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/


#include "mcml_model.h"

#define PARTIALREFLECTION 0     
  /* 1=split photon, 0=statistical reflection. */

#define COSZERO (1.0-1.0E-12)	
  /* cosine of about 1e-6 rad. */

#define COS90D  1.0E-6		
  /* cosine of about 1.57 - 1e-6 rad. */

double RFresnel(double n1, double n2, double ca1, double * ca2Ptr);
double SpinTheta(double g);

void Medium::SelectMedium (Medium::MediumName mediumName) {
  switch (mediumName) {
    case Medium::AIR:
      n = 1.0;
      mua = 0.0;
      mus = 0.0;
      g = 0.0;
      break;
    case Medium::DERMIS:
      n = 1.4;
      mua = 0.26;
      mus = 137.0;
      g = 0.90;
      break;
    case Medium::TYPE_II_EPIDERMIS:
      n = 1.3;
      mua = 5.0;
      mus = 200.0;
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





void ModelInput::SelectModelInput (ModelInput::ModelInputName modelInputName,
    long numPhotonsSet) {
  
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




Photon::Photon (LayerStruct layerObj, double rSpecular) {
  x = 0.0;
  y = 0.0;
  z = 0.0;
  ux = 0.0;
  uy = 0.0;
  uz = 1.0;
  w = 1.0 - rSpecular;
  dead = false;
  layer = 1;
  s = 0;
  sleft = 0;
  
  // take care of the case when the first layer is glass
  if ((layerObj.layer[1].mua == 0.0) && (layerObj.layer[1].mus == 0.0)) {
    layer = 2;      // skip to next layer
    z = layerObj.layerZ[0][2];  // use z0 from the next layer
  }
}


void Photon::RunOnePhoton (MCMLModel * model) {
// run a single photon scattering till completion
  while (!dead) {
    HopDropSpin(model);
  }
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

  double mua, mus, rnd;

  std::uniform_real_distribution<double> unif(0, 1);
  std::default_random_engine re;
  mua = model->layerObj.layer[layer].mua;
  mus = model->layerObj.layer[layer].mus;

  if (sleft == 0.0) {      // make a new step
    rnd = unif(re);
    s = -log(rnd)/(mua + mus);
  } else {                       // take the leftover
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

  std::uniform_real_distribution<double> unif(0, 1);
  std::default_random_engine re;
  
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
  } else if (unif(re) > r) {    // transmitted to layer-1
    layer -= 1;
    ux *= ni/nt;
    uy *= ni/nt;
    uz = -uz1;
  } else			      		// reflected
    uz = -uz0;
#else
  if (unif(re) > r) {       // transmitted to layer-1
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

  std::uniform_real_distribution<double> unif(0, 1);
  std::default_random_engine re;
  
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
  } else if (unif(re) > r) {     // transmitted to layer+1
    layer += 1;
    ux *= ni/nt;
    uy *= ni/nt;
    uz = uz1;
  } else 						// reflected
    uz = -uz0;
#else
  if (unif(re) > r) {	      // transmitted to layer+1
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

  std::uniform_real_distribution<double> unif(0, 1);
  std::default_random_engine re;
  
  cost = SpinTheta(g);
  sint = sqrt(1.0 - cost*cost);
  // sqrt() is faster than sin().

  psi = 2.0*PI*unif(re);       // spin psi 0-2pi
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
  
  std::uniform_real_distribution<double> unif(0, 1);
  std::default_random_engine re;

  if (w == 0.0)	
    dead = true;
  else if (unif(re) < CHANCE)    // survived the roulette.
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

  std::uniform_real_distribution<double> unif(0, 1);
  std::default_random_engine re;
 
  if (g == 0.0) 
    cost = 2*unif(re) - 1;
  else {
    temp = (1 - g*g)/(1 - g + 2*g*unif(re));
    cost = (1 + g*g - temp*temp)/(2*g);
    if (cost < -1)
      cost = -1.0;
    else if (cost > 1)
      cost = 1.0;
  }
  return cost;
}
