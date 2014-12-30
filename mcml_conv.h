/***********************************************************
 *	Program Name: conv.
 *
 *	A program used to process the data from mcml -
 *	A Monte Carlo simulation of photon distribution in
 *	multi-layered turbid media in ANSI Standard C.
 ****
 *	Creation Date:	11/1991.
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

#ifndef __MCML_CONV_H__
#define __MCML_CONV_H__

#include <math.h>
#include <new>
#include <string>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include "mcml_model.h"

#define GAUSSLIMIT 4
#define PI 3.1415926

#define SIGN(x) ((x)>=0 ? 1:-1)

/****************** Classes *****************************/

class Beam {
/* Beam class - incident light beam class
        Parameters to describe a photon beam.
        Pencil: infinitely narrow beam. This is default for the
            beam from the mcml output.
        Flat:	Flat beam with radius R.
        Gaussian:	Gaussian with 1/e2 radius R.
        Others: general beam described by points with interpolation.
        Class instance variables:
            type - incident beam type, FLAT or GAUSSIAN
            P - total beam power/energy [W or J]
            R - beam radius, defined as 1/e^2 for Gaussian beam [cm]
        Methods:
            
*/
  public:
    enum BeamType {
        FLAT,
        GAUSSIAN
    };
    BeamType type;			// beam type
    double P;       // total power. [J or W]
    double R;       // radius. [cm]

    enum BeamName {
        TRIA_HRL,
        TRIA_FAN
    };
    
    void SelectBeam (Beam::BeamName beamName = Beam::TRIA_HRL);
};




class Node {
/* Node class - node link list binary tree class
        Data structures for the binary tree used to store part of
        the integrand evaluation.
        
        Class instance variables:
            x - x grid node position
            y - y grid node position
            left - left node pointer
            right - right node pointer          
        Methods:            
*/

  public:
    double x;
    double y;
    Node * left;
    Node * right;

    Node () : left (nullptr), right (nullptr) {};
    static Node * FillNode(double x, double y);
    static Node * SearchNode(Node * tree, double x);
    static void InsertNode(Node * tree, double x, double y);
    static void FreeTree(Node * tree);
};



class ConvVar {
/* ConvVar class - convoluation variables class
        A global structure to pass the current coordinate of the
        physical quantities being evaluated and the pointers of the
        input and output parameters to the integration function.
        
        Class instance variables:
            r - r position
            iz - iz index
            ia - ia index
            tree - A tree to store ITheta() & ExpBessI0().           
        Methods:
            
*/
    
  public:
    double r;
    short iz;
    short ia;
    Node * tree;      // A tree to store ITheta() & ExpBessI0().

    ConvVar () : tree (nullptr) {};
};



class ConvInput {
/* ConvInput class - beam convolution input class
        Input parameters for each independent run.
        z and r are for the cylindrical coordinate system. [cm]
        a is for the angle alpha between the photon exiting
        direction and the surface normal. [radian]
        The grid line separations in z, r, and alpha
        directions are dz, dr, and da respectively.  The numbers
        of grid lines in z, r, and alpha directions are
        nz, nr, and na respectively.
        The member layerspecs will point to an array of
        structures which store parameters of each layer.
        This array has (number_layers + 2) elements. One
        element is for a layer.
        The layers 0 and (num_layers + 1) are for top ambient
        medium and the bottom ambient medium respectively.
        For convolution, the grid line separations in z, and alpha
        directions are still dz, and da respectively.  The numbers
        of grid lines in z, and alpha directions are still
        nz, and na respectively. However, the grid line separation
        and the number of grid lines in r direction are drc and
        nrc respectively.
        Class instance variables:
            beam - incident beam class instance object
            drc - convolution r grid separation.[cm]
            nrc - convolution array range 0..nrc-1.
            eps - relative error in convolution           
        Methods:
            
*/

  public:
    Beam beam;    // incident beam of finite size
    double drc;
    short nrc;
    MCMLModel mcmlModel;
    ConvVar convVar;
    
    enum ConvName {
        TRIA_HRL,
        TRIA_FAN
    };
    
    void SelectConvInput (MCMLModel mcmlModelSet, 
          ConvInput::ConvName convName = ConvInput::TRIA_HRL);
    void FreeConvInput();       
};     



class MCMLConv : public ConvInput {
/* MCMLConv class - multi-layered photon scattering model beam convolution
        inherits from ConvInput beam setup
        Structures for scored physical quantities
        from mcml and to be convolved for photon
        beams of finite size.  Therefore, "Out"
        here means the output of both mcml and conv.
        The member allocated is used to keep the status
        of the arrays.  It is set to 1 if all the arrays
        are allocated and assigned values.  It is set to
        0 otherwise.
        z and r represent z and r coordinates of the
        cylindrical coordinate system. [cm]
        a is the angle alpha between the photon exiting
        direction and the normal to the surfaces. [radian]
        See comments of the InputStruct.
        See manual for the physcial quantities.
        Class instance variables:
            Rd_rac - convolved data. [J/(cm2 sr)]
            Rd_rc - 1D radial distribution of diffuse reflectance [J/cm2]
            A_rzc - 2D probability density in turbid media over r & z [J/cm3]
            Tt_rac - 2D distribution of total transmittance [J/(cm2 sr)]
            Tt_rc - 1D radial distribution of transmittance [J/cm2]
        Methods:
            
*/
  private:
    void ConvRd_ra ();
    void ConvRd_r ();
    void ConvA_rz ();
    void ConvTt_ra ();
    void ConvTt_r ();    
  
  public:
    double ** Rd_rac;
    double * Rd_rc;
    double ** A_rzc;
    double ** Tt_rac;
    double * Tt_rc;

    MCMLConv () : Rd_rac (nullptr), Rd_rc (nullptr), A_rzc (nullptr), 
      Tt_rac (nullptr), Tt_rc (nullptr) {};
    void SelectMCMLConv (MCMLModel mcmlModelSet, std::string convName);
    void FreeMCMLConv ();
    void RunConv ();      
};



#endif    //__MCML_CONV_H__

