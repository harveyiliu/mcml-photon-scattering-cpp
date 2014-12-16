/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	Monte Carlo simulation of photon distribution in 
 *	multi-layered turbid media in ANSI Standard C.
 ****
 *	Starting Date:	10/1991.
 *	Current Date:	6/1992.
 *
 *	Lihong Wang, Ph. D.
 *	Steven L. Jacques, Ph. D.
 *	Laser Biology Research Laboratory - 17
 *	M.D. Anderson Cancer Center
 *	University of Texas
 *	1515 Holcombe Blvd.
 *	Houston, TX 77030
 *	USA
 *
 *	This program was based on:
 *	(1) The Pascal code written by Marleen Keijzer and 
 *	Steven L. Jacques in this laboratory in 1989, which
 *	deals with multi-layered turbid media.
 *
 *	(2) Algorithm for semi-infinite turbid medium by 
 *	S.A. Prahl, M. Keijzer, S.L. Jacques, A.J. Welch, 
 *	SPIE Institute Series Vol. IS 5 (1989), and by 
 *	A.N. Witt, The Astrophysical journal Supplement
 *	Series 35, 1-6 (1977).
 *	
 *	Major modifications include:
 *		. Conform to ANSI Standard C.
 *		. Removal of limit on number of array elements, 
 *		  because arrays in this program are dynamically 
 *		  allocated. This means that the program can accept 
 *		  any number of layers or gridlines as long as the 
 *		  memory permits.
 *		. Avoiding global variables whenever possible.  This
 *		  program has not used global variables so far.
 *		. Grouping variables logically using structures.
 *		. Top-down design, keep each subroutine clear & 
 *		  short.	
 *		. Reflectance and transmittance are angularly 
 *		  resolved.
 ****
 *	General Naming Conventions:
 *	Preprocessor names: all capital letters, 
 *		e.g. #define PREPROCESSORS
 *	Globals: first letter of each word is capital, no 
 *		underscores, 
 *		e.g. short GlobalVar;
 *	Dummy variables:  first letter of each word is capital,
 *		and words are connected by underscores, 
 *		e.g. void NiceFunction(char Dummy_Var);
 *	Local variables:  all lower cases, words are connected 
 *		by underscores,
 *		e.g. short local_var;
 *	Function names or data types:  same as Globals.
 *
 ****
 *	Dimension of length: cm.
 *
 ****/

#ifndef __MCML_MODEL_H__
#define __MCML_MODEL_H__

#include <math.h>
#include <new>
#include <string>
#include <algorithm>
#include <random>

#define PI 3.1415926
#define WEIGHT 1e-4		/* Critical weight for roulette. */
#define CHANCE 0.1		/* Chance of roulette survival. */

#define SIGN(x) ((x)>=0 ? 1:-1)

/****************** Classes *****************************/

/****
Class used to describe an optical medium.
Medium class - optical medium class defining the optical properties
Class instance variables:
  n - refractive index
  mua - absorption coefficient. [1/cm]
  mus - scattering coefficient. [1/cm]
  g - anisotropy
Methods:
 ****/

class Medium {
  friend class LayerStruct;
  friend class Photon;
  private:
    double n;			/* refractive index of a layer. */
    double mua;	    /* absorption coefficient. [1/cm] */
    double mus;	    /* scattering coefficient. [1/cm] */
    double g;		    /* anisotropy. */
  public:
    enum MediumName {
        AIR,
        DERMIS,
        TYPE_II_EPIDERMIS
    };
    void SelectMedium (Medium::MediumName mediumName = Medium::DERMIS);    
};


/****
 *	Class used to describe the geometry and optical
 *	properties of a layer.
 *	z0 and z1 are the z coordinates for the upper boundary
 *	and lower boundary respectively.
 *
 *	cos_crit0 and cos_crit1 are the cosines of the 
 *	critical angle of total internal reflection for the
 *	upper boundary and lower boundary respectively.
 *	They are set to zero if no total internal reflection
 *	exists.
 *	They are used for computation speed.

LayerStruct class - multi-layered structure
Class instance variables:
  nIn - refractive index of the incidence medium
  nOut - refractive index of the exit medium
  numLayers - number of layers
  layer - list of layer objects
  layerThickness - layer thickness in [cm]
  layerZ - layer depth z coordinates, top and bottom [cm]
  cosCrit - ciritical angle cosines of each layer, top and bottom
Methods:
 ****/

class LayerStruct {
  friend class Photon;
  friend class ModelInput;
  friend class MCMLModel;
  private:
    short	numLayers;			/* number of layers */
    Medium * layer;     /* layer medium list, 1st layer is incidence medium
                          last layer is the exit medium */
    double * layerThickness;  /* layer thickness array [cm], not include
                          the incidence and exit media */
    double * layerZ[2];	  /* layer z coordinates, top and bottom [cm] */
    double * cosCrit[2];  /* layer ciritical angle cosines, top and bottom */

  public:  
    enum LayerName {
        BARE_DERMIS,
        TYPE_II_SKIN
    };
    LayerStruct() : layer (nullptr), layerThickness (nullptr) {
        layerZ[0] = nullptr;
        layerZ[1] = nullptr;
        cosCrit[0] = nullptr;
        cosCrit[1] = nullptr;
    }
    void SelectLayerStruct (LayerStruct::LayerName layerName =
      LayerStruct::BARE_DERMIS);
    void FreeLayerStruct ();
    double CalcRSpecular ();
};




/****
 *	Model input parameters for each independent run.
 *
 *	z and r are for the cylindrical coordinate system. [cm]
 *	a is for the angle alpha between the photon exiting 
 *	direction and the surface normal. [radian]
 *
 *	The grid line separations in z, r, and alpha
 *	directions are dz, dr, and da respectively.  The numbers 
 *	of grid lines in z, r, and alpha directions are
 *	nz, nr, and na respectively.
 *
 *	The member layerspecs will point to an array of 
 *	structures which store parameters of each layer. 
 *	This array has (number_layers + 2) elements. One
 *	element is for a layer.
 *	The layers 0 and (num_layers + 1) are for top ambient 
 *	medium and the bottom ambient medium respectively.

ModelInput class - multi-layered photon scattering model input
Class instance variables:
  numPhotons - number of photons to be traced
  Wth - play roulette if photon weight < Wth
  dz - z grid separation [cm]
  dr - r grid separation [cm]
  da - alpha grid separation [radian]
  nz - array range 0..nz-1
  nr - array range 0..nr-1
  na - array range 0..na-1
  layerObj - medium layer structure class instance
Methods:
 ****/

class ModelInput {
  friend class Photon;
  public:
    double Wth; 				/* play roulette if photon */
							      /* weight < Wth.*/
  
    double dz;				/* z grid separation.[cm] */ 
    double dr;				/* r grid separation.[cm] */
    double da;				/* alpha grid separation. */
							/* [radian] */
    short nz;					/* array range 0..nz-1. */
    short nr;					/* array range 0..nr-1. */
    short na;					/* array range 0..na-1. */
    LayerStruct layerObj;	  /* layer class object with parameters. */

    enum ModelInputName {
        BARE_DERMIS,
        TYPE_II_SKIN
    };
    void SelectModelInput (ModelInput::ModelInputName modelInputName =
      ModelInput::BARE_DERMIS);
    void FreeModelInput ();  	
};


/****
 *	Class for scoring physical quantities. 
 *	z and r represent z and r coordinates of the 
 *	cylindrical coordinate system. [cm]
 *	a is the angle alpha between the photon exiting 
 *	direction and the normal to the surfaces. [radian]
 *	See comments of the InputStruct.
 *	See manual for the physcial quantities.

ModelOutput class - multi-layered photon scattering model output
Class instance variables:
  Rsp - specular reflectance [-]
  Rd - total diffuse reflectance [-]
  A - total absorption probability [-]
  Tt - total transmittance [-]
  Rd_ra - 2D distribution of diffuse reflectance [1/(cm2 sr)]
  Rd_r - 1D radial distribution of diffuse reflectance [1/cm2]
  Rd_a - 1D angular distribution of diffuse reflectance [1/sr]
  A_rz - 2D probability density in turbid media over r & z [1/cm3]
  A_z - 1D probability density over z [1/cm]
  A_l - each layer's absorption probability [-]
  Tt_ra - 2D distribution of total transmittance [1/(cm2 sr)]
  Tt_r - 1D radial distribution of transmittance [1/cm2]
  Tt_a - 1D angular distribution of transmittance [1/sr]
Methods:
 ****/

class MCMLModel : public ModelInput {
  friend class Photon;
  private:
    void Sum2DRd();
    short IzToLayer(short iz);
    void Sum2DA();
    void Sum2DTt();
    void ScaleRdTt();
    void ScaleA(); 
  public:
    long numPhotons;    // number of photons traced
    double Rsp;	// specular reflectance. [-]
    double ** Rd_ra;	  // 2D distribution of diffusereflectance. [1/(cm2 sr)]
    double * Rd_r;	  // 1D radial distribution of diffuse reflectance. [1/cm2]
    double * Rd_a;	  // 1D angular distribution of diffuse reflectance. [1/sr]
    double Rd;		  // total diffuse reflectance. [-]  
    double ** A_rz; // 2D probability density in turbid media over r & z.
                    // [1/cm3]
    double * A_z;	// 1D probability density over z. [1/cm]
    double * A_l;	// each layer's absorption probability. [-]
    double A;		  // total absorption probability. [-]  
    double ** Tt_ra;	// 2D distribution of total transmittance. [1/(cm2 sr)]
    double * Tt_r;	  // 1D radial distribution of transmittance. [1/cm2]
    double * Tt_a;	  // 1D angular distribution of transmittance. [1/sr]
    double Tt;		// total transmittance. [-]


    MCMLModel () : Rd_ra (nullptr), Rd_r (nullptr), Rd_a (nullptr),
      A_rz (nullptr), A_z (nullptr), A_l (nullptr),
      Tt_ra (nullptr), Tt_r (nullptr), Tt_a (nullptr) {};
    void SelectMCMLModel (std::string modelName);
    void FreeMCMLModel ();
    void DoOneRun (long numPhotons);
    void SumScaleResult();
};



/****
 *	Structure used to describe a photon packet.

Photon class - MCML photon class for Monte Carlo scattering model in
multilayered turbid media. 
Class instance variables:
  x = Cartesian coordinate x [cm]
  y = Cartesian coordinate y [cm]
  z = Cartesian coordinate z [cm]
  ux = directional cosine x of a photon
  uy = directional cosine y of a photon
  uz = directional cosine z of a photon
  w - weight
  dead - true if photon is terminated
  layer - index to layer where the photon packet resides
  s - current step size [cm]
  sleft - step size left, dimensionless [-]           
Methods:
 ****/

class Photon {
  friend class MCMLModel;
  private:
    double x, y ,z;	/* Cartesian coordinates.[cm] */
    double ux, uy, uz;/* directional cosines of a photon. */
    double w;			/* weight. */
    bool dead;		/* true if photon is terminated. */
    short layer;		/* index to layer where the photon */
					/* packet resides. */
    double s;			/* current step size. [cm]. */
    double sleft;		/* step size left. dimensionless [-]. */

    void HopDropSpin(MCMLModel * model);
    void HopInGlass(MCMLModel * model);
    void HopDropSpinInTissue(MCMLModel * model);
    void StepSizeInGlass(MCMLModel * model);
    void StepSizeInTissue(MCMLModel * model);
    void Hop();
    void CrossOrNot(MCMLModel * model);
    void CrossUpOrNot(MCMLModel * model);
    void CrossDnOrNot(MCMLModel * model);
    bool HitBoundary(MCMLModel * model);
    void Drop(MCMLModel * model);
    void Spin(double g);
    void RecordR(MCMLModel * model, double refl);
    void RecordT(MCMLModel * model, double refl);
    void Roulette();    

  public:
    void Reset(MCMLModel * model);
    void RunOnePhoton(MCMLModel * model);   
};




#endif    //__MCML_MODEL_H__

