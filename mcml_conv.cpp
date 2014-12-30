/***********************************************************
 *	class setup program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/

#include "mcml_conv.h"

double ITheta(double r, double r2, double R);
double ExpBessI0(double r, double r2, double R);
double RT_raInterp(double r2, double ** RT_ra, MCMLConv * mcmlConv);
double Rd_raFGIntegrand(double r2, void * params);
double IntegrateQuad(double (*func) (double x, void * params), 
          double a, double b, MCMLConv * mcmlConv);
double FlatIntegration(double (*func) (double x, void * params), 
          MCMLConv * mcmlConv);
double GaussIntegration(double (*func) (double x, void * params), 
          MCMLConv * mcmlConv);


void Beam::SelectBeam (Beam::BeamName beamName) {
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
  switch (beamName) {
    case Beam::TRIA_HRL:
      type = Beam::FLAT;
      P = 20;     // total power. [J or W]
      R = 0.5;    // radius. [cm]
      break;
    case Beam::TRIA_FAN:
      type = Beam::GAUSSIAN;
      P = 0.012;
      R = 0.025;
      break;
    default:
      type = Beam::FLAT;
      P = 20;
      R = 0.5;
  }
}



Node * Node::FillNode(double x, double y) {
  Node * l = new Node;
  l->x = x;
  l->y = y;
  l->left = nullptr;
  l->right = nullptr;
  return l;
}



Node * Node::SearchNode(Node * tree, double x) {
  Node * l = tree;
  bool found = false;

  while (l != nullptr && !found) {
    if (x < l->x)
      l = l->left;
    else if (x > l->x)
      l = l->right;
    else
      found = true;
  }
  return l;
}



void Node::InsertNode(Node * tree, double x, double y) {
  Node * l1 = nullptr;
  Node * l2 = tree;
  while (l2 != nullptr) {
    l1 = l2;
    if (x < l2->x)
      l2 = l2->left;
    else
      l2 = l2->right;
  }

  if (l1 == nullptr)		// Empty tree.
    tree = Node::FillNode(x, y);
  else if (x < l1->x)
    l1->left = Node::FillNode(x, y);
  else
    l1->right = Node::FillNode(x, y);
}


void Node::FreeTree(Node * tree) {
  if (tree != nullptr) {
    FreeTree(tree->left);
    FreeTree(tree->right);
    delete tree;    
  }
}



void ConvInput::SelectConvInput (MCMLModel mcmlModelSet, 
      ConvInput::ConvName convName) {
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
  this->FreeConvInput();    // free previously set array members if needed
  switch (convName) {
    case ConvInput::TRIA_HRL:    
      beam.SelectBeam(Beam::TRIA_HRL);    // incident beam of finite size
      drc = 0.005;        // convolution r grid separation.[cm]
      nrc = 150;          // convolution array range 0..nrc-1.  
      break; 
    case ConvInput::TRIA_FAN:
      beam.SelectBeam(Beam::TRIA_FAN);
      drc = 0.001;
      nrc = 100;
      break;
    default:
      beam.SelectBeam(Beam::TRIA_HRL);
      drc = 0.005;
      nrc = 150;      
  }
  mcmlModel = mcmlModelSet;
}


void ConvInput::FreeConvInput () {
  mcmlModel.FreeMCMLModel();
  if (convVar.tree != nullptr) {
    Node::FreeTree(convVar.tree);
    convVar.tree = nullptr;
  }    
}





void MCMLConv::SelectMCMLConv (MCMLModel mcmlModelSet, std::string convName) {
  
  this->FreeMCMLConv();

  if (convName.compare("TRIA_HRL") == 0)
    this->SelectConvInput (mcmlModelSet, ConvInput::TRIA_HRL);
  else if (convName.compare("TRIA_FAN") == 0)
    this->SelectConvInput (mcmlModelSet, ConvInput::TRIA_FAN);
  else
    this->SelectConvInput (mcmlModelSet, ConvInput::TRIA_HRL);
   
  short na = mcmlModel.na;
  short nz = mcmlModel.nz;
  Rd_rac = new double * [nrc];
  for (int i = 0; i < nrc; i++) {
    Rd_rac[i] = new double [na];
    // initialize array with 0
    std::fill (Rd_rac[i], Rd_rac[i]+na-1, 0);  
  }

  Rd_rc = new double [nrc];
  std::fill (Rd_rc, Rd_rc+nrc-1, 0);

  A_rzc = new double * [nrc];
  for (int i = 0; i < nrc; i++) {
    A_rzc[i] = new double [nz];
    std::fill (A_rzc[i], A_rzc[i]+nz-1, 0);  
  }

  Tt_rac = new double * [nrc];
  for (int i = 0; i < nrc; i++) {
    Tt_rac[i] = new double [na];
    std::fill (Tt_rac[i], Tt_rac[i]+na-1, 0);    
  }

  Tt_rc = new double [nrc];
  std::fill (Tt_rc, Tt_rc+nrc-1, 0);

}


void MCMLConv::FreeMCMLConv () {
  if (Rd_rac != nullptr) {
    for (int i = 0; i < nrc; i++) 
      delete[] Rd_rac[i];
  }
  Rd_rac = nullptr;
  if (Rd_rc != nullptr) 
    delete[] Rd_rc;
  Rd_rc = nullptr;
  if (A_rzc != nullptr) {
    for (int i = 0; i < nrc; i++) 
      delete[] A_rzc[i];
  }
  A_rzc = nullptr;
  if (Tt_rac != nullptr) {
    for (int i = 0; i < nrc; i++) 
      delete[] Tt_rac[i];
  }
  Tt_rac = nullptr;
  if (Tt_rc != nullptr) 
    delete[] Tt_rc;
  Tt_rc = nullptr;  
  this->FreeConvInput();  
}




void MCMLConv::RunConv () {
  ConvRd_ra();
  ConvRd_r();
  ConvA_rz();
  ConvTt_ra();
  ConvTt_r();
}



void MCMLConv::ConvRd_ra () {
  short irc, ia;
  double rc;
  double P = beam.P;
  double R = beam.R;
  for (irc = 0; irc < nrc; irc++) {
    rc = (irc + 0.5)*drc;
    convVar.r = rc;
    convVar.tree = nullptr;    	// init the tree
    for (ia = 0; ia < mcmlModel.na; ia++) {
      convVar.ia = ia;
      if (beam.type == Beam::FLAT)
        Rd_rac[irc][ia] = 2*P/(R*R)*FlatIntegration(&Rd_raFGIntegrand, this);
      else      // Gaussian
        Rd_rac[irc][ia] = 4*P/(R*R)*GaussIntegration(&Rd_raFGIntegrand, this);
    }
  }
}


void MCMLConv::ConvRd_r () {

}


void MCMLConv::ConvA_rz () {

}


void MCMLConv::ConvTt_ra () {

}


void MCMLConv::ConvTt_r () {

}



double RT_raInterp(double r2, double ** RT_ra, MCMLConv * mcmlConv) {
// Interpolate for the arrays Rd_ra[] or Tt_ra[].
  short nr = mcmlConv->mcmlModel.nr;
  short ia = mcmlConv->convVar.ia;
  short ir2lo;
  double RT_lo, RT_hi, RT_at_r2;
  double ir2 = r2/mcmlConv->mcmlModel.dr;

  if (nr < 3)
    RT_at_r2 = RT_ra[0][ia];
  else if (ir2 < (nr - 1.5)) {      	// interpolation
    ir2lo = std::max<short>(0, (short) (ir2 - 0.5));	    // truncation
    RT_lo = RT_ra[ir2lo][ia];
    RT_hi = RT_ra[ir2lo + 1][ia];
    RT_at_r2 = RT_lo + (RT_hi - RT_lo)*(ir2 - ir2lo - 0.5);
  }
  else {			// extrapolation
    ir2lo = nr - 3;
    RT_lo = RT_ra[ir2lo][ia];
    RT_hi = RT_ra[ir2lo + 1][ia];
    if (RT_lo >= RT_hi)		// Noise test
      RT_at_r2 = RT_lo + (RT_hi - RT_lo)*(ir2 - ir2lo - 0.5);
    else
      RT_at_r2 = 0.0;
  }
  return std::max<double>(0, RT_at_r2);
}




double ITheta(double r, double r2, double R) {
  double temp;

  if (R >= (r + r2))
    temp = 1;
  else if (fabs(r - r2) <= R) {
    temp = (r*r + r2*r2 - R*R)/(2*r*r2);
    if (fabs(temp) > 1)
      temp = SIGN(temp);
    temp = acos(temp)/PI;
  }
  else			// R < fabs(r-r2)
    temp = 0;
  return temp;
}




double ExpBessI0(double r, double r2, double R) {
    double _RR = 1/(R*R);
    double x = 4*r*r2*_RR;
    double y = 2*(r2*r2 + r*r)*_RR;
    double expbess = exp(-y + x)*gsl_sf_bessel_j0 (x);
    return expbess;
}



double Rd_raFGIntegrand(double r2, void * params) {
/* Convolution integrand for either flat or gaussian beams.
    See comments for A_rzFGIntegrand().
    r" in the integration. */
  
  MCMLConv * mcmlConv = (MCMLConv *) params;

  double f;
  double ** RT_ra = mcmlConv->mcmlModel.Rd_ra;
  double Rd_at_r2 = RT_raInterp(r2, RT_ra, mcmlConv);
  double R = mcmlConv->beam.R;
  double r = mcmlConv->convVar.r;
  Node * tree = mcmlConv->convVar.tree;
  Node * link = Node::SearchNode(tree, r2);
  if (link != nullptr)	    // f in tree.
    f = link->y;
  else {
    if (mcmlConv->beam.type == Beam::FLAT)
      f = ITheta(r, r2, R);
    else			// Gaussian
      f = ExpBessI0(r, r2, R);
    Node::InsertNode(tree, r2, f);
  }
  f *= Rd_at_r2*r2;
  return f;
}




double IntegrateQuad(double (*func) (double x, void * params), 
          double a, double b, MCMLConv * mcmlConv) {

  gsl_integration_workspace * workPtr 
      = gsl_integration_workspace_alloc (1000);

  double absError = 1.0e-5;   // to avoid round-off problems
  double relError = 1.0e-5;   // the result will usually be much better
  double result;              // the result from the integration
  double error;               // the estimated error from the integration

  gsl_function integralFunction;
  void * paramsPtr = mcmlConv;

  integralFunction.function = func;
  integralFunction.params = paramsPtr;

  gsl_integration_qags (&integralFunction, a, b, absError, relError, 1000,
        workPtr, &result, &error);

  gsl_integration_workspace_free (workPtr);

  return result;

}




double FlatIntegration(double (*func) (double x, void * params), 
          MCMLConv * mcmlConv) {
    double rc = mcmlConv->convVar.r;
    double R = mcmlConv->beam.R;
    double b_max = (mcmlConv->mcmlModel.nr - 0.5)*mcmlConv->mcmlModel.dr;
    double a = std::max<double>(0, rc - R);
    double b = std::min<double>(b_max, rc + R);

    if (a >= b)
      return 0;
    else
      return IntegrateQuad(func, a, b, mcmlConv);    
}



double GaussIntegration(double (*func) (double x, void * params), 
          MCMLConv * mcmlConv) {
// Used by convolution over Gaussian beam.  Ignore the value
// beyond GAUSSLIMIT radius.

  double rc = mcmlConv->convVar.r;
  double R = mcmlConv->beam.R;
  double b_max = (mcmlConv->mcmlModel.nr - 0.5)*mcmlConv->mcmlModel.dr;
  double a = std::max<double>(0, (rc - GAUSSLIMIT * R));
  double b = std::min<double>(b_max, (rc + GAUSSLIMIT * R));

  if (a >= b)
    return 0;
  else
    return IntegrateQuad(func, a, b, mcmlConv);
}

