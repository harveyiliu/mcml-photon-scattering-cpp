/***********************************************************
 *	class setup program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/

#include "mcml_conv.h"

double ITheta(double r, double r2, double R);
double ModifiedBessI0(double x);
double ExpBessI0(double r, double r2, double R);
double RT_raInterp(double r2, double ** RT_ra, MCMLConv * mcmlConv);
double RT_rInterp(double r2, double * RT_r, MCMLConv * mcmlConv);
double A_rzInterp(double r2, MCMLConv * mcmlConv);
double Rd_raFGIntegrand(double r2, void * params);
double Rd_rFGIntegrand(double r2, void * params);
double A_rzFGIntegrand(double r2, void * params);
double Tt_raFGIntegrand(double r2, void * params);
double Tt_rFGIntegrand(double r2, void * params);
double IntegrateQuad(double (*func) (double x, void * params), 
          double a, double b, MCMLConv * mcmlConv);
double FlatIntegration(double (*func) (double x, void * params), 
          MCMLConv * mcmlConv);
double GaussIntegration(double (*func) (double x, void * params), 
          MCMLConv * mcmlConv);





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
      Beam::BeamType beamType, double P, double R) {
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
  switch (beamType) {
    case Beam::FLAT:    
      beam.type = Beam::FLAT;    // incident beam of finite size  
      break; 
    case Beam::GAUSSIAN:
      beam.type = Beam::GAUSSIAN;
      break;
    default:
      beam.type = Beam::FLAT;      
  }
  mcmlModel = mcmlModelSet;
  beam.P = P;
  beam.R = std::max<double>(R, 1e-5);       // minimum radius required
  double minR = 1.2*(R + mcmlModel.nr*mcmlModel.dr);
  drc = std::max<double>(mcmlModel.dr, minR/200);
  nrc = (short) (minR/drc);
}


void ConvInput::FreeConvInput () {
  mcmlModel.FreeMCMLModel();
  if (convVar.tree != nullptr) {
    Node::FreeTree(convVar.tree);
    convVar.tree = nullptr;
  }    
}





void MCMLConv::SelectMCMLConv (MCMLModel mcmlModelSet, std::string beamType,
        double P, double R) {
  
  this->FreeMCMLConv();

  if (beamType.compare("FLAT") == 0)
    this->SelectConvInput (mcmlModelSet, Beam::FLAT, P, R);
  else if (beamType.compare("GAUSSIAN") == 0)
    this->SelectConvInput (mcmlModelSet, Beam::GAUSSIAN, P, R);
  else
    this->SelectConvInput (mcmlModelSet, Beam::FLAT, P, R);
   
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

  F_rzc = new double * [nrc];
  for (int i = 0; i < nrc; i++) {
    F_rzc[i] = new double [nz];
    std::fill (F_rzc[i], F_rzc[i]+nz-1, 0);  
  }

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
  if (F_rzc != nullptr) {
    for (int i = 0; i < nrc; i++) 
      delete[] F_rzc[i];
  }
  F_rzc = nullptr;  
  this->FreeConvInput();  
}




void MCMLConv::RunConv () {
  ConvRd_ra();
  ConvRd_r();
  ConvA_rz();
  ConvTt_ra();
  ConvTt_r();
  ConvA2F();
}


void MCMLConv::ConvA2F () {
  short irc, iz;
  short nz = mcmlModel.nz;
  double mua;
  
  for (irc = 0; irc < nrc; irc++) {
    for (iz = 0; iz < nz; iz++) {
      mua = mcmlModel.GetMuaAtIz(iz);
      if (mua > 0.0)
        F_rzc[irc][iz] = A_rzc[irc][iz]/mua;     // F in J/cm2
    }
  }
}





double MCMLConv::CenterHalfMaxDepth () {
  short iz;
  short nz = mcmlModel.nz;
  double dz = mcmlModel.dz;
  double depth = 0;

  for (iz = 0; iz < nz; iz++) {
    if (F_rzc[0][iz] <= 0.5*F_rzc[0][0]) {
      depth = (iz + 0.5)*dz;
      break;
    }
  }
  return depth;
}



double MCMLConv::SurfaceHalfMaxWidth () {
  short irc;
  double width = 0;

  for (irc = 0; irc < nrc; irc++) {
    if (F_rzc[irc][0] <= 0.5*F_rzc[0][0]) {
      width = 2*(irc + 0.5)*drc;
      break;
    }
  }
  return width;
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
    Node::FreeTree(convVar.tree);
  }
}


void MCMLConv::ConvRd_r () {
  short irc;
  double rc;
  double P = beam.P;
  double R = beam.R;

  for (irc = 0; irc < nrc; irc++) {
    rc = (irc + 0.5)*drc;
    convVar.r = rc;
    if (beam.type == Beam::FLAT)
      Rd_rc[irc] = 2*P/(R*R)*FlatIntegration(&Rd_rFGIntegrand, this);
    else       // Gaussian
      Rd_rc[irc] = 4*P/(R*R)*GaussIntegration(&Rd_rFGIntegrand, this);
  }
}


void MCMLConv::ConvA_rz () {
  short irc, iz;
  double rc;
  double P = beam.P;
  double R = beam.R;

  for (irc = 0; irc < nrc; irc++) {
    rc = (irc + 0.5)*drc;
    convVar.r = rc;
    convVar.tree = nullptr;    	// init the tree
    for (iz = 0; iz < mcmlModel.nz; iz++) {
      convVar.iz = iz;
      if (beam.type == Beam::FLAT)
        A_rzc[irc][iz] = 2*P/(R*R)*FlatIntegration(&A_rzFGIntegrand, this);
      else       // Gaussian
        A_rzc[irc][iz] = 4*P/(R*R)*GaussIntegration(&A_rzFGIntegrand, this);
    }
    Node::FreeTree(convVar.tree);
  }
}


void MCMLConv::ConvTt_ra () {
  short irc, ia;
  double rc;
  double P = beam.P;
  double R = beam.R;

  for (irc = 0; irc < nrc; irc++) {
    rc = (irc + 0.5)*drc;
    convVar.r = rc;
    convVar.tree = nullptr;        // init the tree
    for (ia = 0; ia < mcmlModel.na; ia++) {
      convVar.ia = ia;
      if (beam.type == Beam::FLAT)
        Tt_rac[irc][ia] = 2*P/(R*R)*FlatIntegration(&Tt_raFGIntegrand, this);
      else       // Gaussian
        Tt_rac[irc][ia] = 4*P/(R*R)*GaussIntegration(&Tt_raFGIntegrand, this);
    }
    Node::FreeTree(convVar.tree);
  }
}


void MCMLConv::ConvTt_r () {
  short irc;
  double rc;
  double P = beam.P;
  double R = beam.R;

  for (irc = 0; irc < nrc; irc++) {
    rc = (irc + 0.5)*drc;
    convVar.r = rc;
    if (beam.type == Beam::FLAT)
      Tt_rc[irc] = 2*P/(R*R)*FlatIntegration(&Tt_rFGIntegrand, this);
    else	       // Gaussian
      Tt_rc[irc] = 4*P/(R*R)*GaussIntegration(&Tt_rFGIntegrand, this);
  }
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



double RT_rInterp(double r2, double * RT_r, MCMLConv * mcmlConv) {
// Interpolate for the arrays Rd_r[] or Tt_r[].
  short nr = mcmlConv->mcmlModel.nr;
  short ir2lo;
  double RT_lo, RT_hi, RT_at_r2;
  double ir2 = r2/mcmlConv->mcmlModel.dr;
  
  if (nr < 3)
    RT_at_r2 = RT_r[0];
  else if (ir2 < (nr - 1.5)) {      // interpolation
    ir2lo = std::max<short>(0, (short) (ir2 - 0.5));       // truncation
    RT_lo = RT_r[ir2lo];
    RT_hi = RT_r[ir2lo + 1];
    RT_at_r2 = RT_lo + (RT_hi - RT_lo)*(ir2 - ir2lo - 0.5);
  }
  else {           // extrapolation
    ir2lo = nr - 3;
    RT_lo = RT_r[ir2lo];
    RT_hi = RT_r[ir2lo + 1];
    if (RT_lo >= RT_hi)      // Noise test
      RT_at_r2 = RT_lo + (RT_hi - RT_lo)*(ir2 - ir2lo - 0.5);
    else
      RT_at_r2 = 0.0;
  }
    return std::max<double>(0, RT_at_r2);
}



double A_rzInterp(double r2, MCMLConv * mcmlConv) {
// Interpolate for the arrays A_rz[]
  double ** A_rz = mcmlConv->mcmlModel.A_rz;
  short nr = mcmlConv->mcmlModel.nr;
  short iz = mcmlConv->convVar.iz;
  short ir2lo;
  double A_lo, A_hi, A_at_r2;
  double ir2 = r2/mcmlConv->mcmlModel.dr;
  
  if (nr < 3)
    A_at_r2 = A_rz[0][iz];
  else if (ir2 < (nr - 1.5)) {       // interpolation
    ir2lo = std::max<short>(0, (short) (ir2 - 0.5));   // truncation
    A_lo = A_rz[ir2lo][iz];
    A_hi = A_rz[ir2lo + 1][iz];
    A_at_r2 = A_lo + (A_hi - A_lo)*(ir2 - ir2lo - 0.5);
  }
  else {       // extrapolation
    ir2lo = nr - 3;
    A_lo = A_rz[ir2lo][iz];
    A_hi = A_rz[ir2lo + 1][iz];
    if (A_lo >= A_hi)       // Noise test
      A_at_r2 = A_lo + (A_hi - A_lo)*(ir2 - ir2lo - 0.5);
    else
      A_at_r2 = 0.0;
  }
  return std::max<double>(0, A_at_r2);
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



double ModifiedBessI0(double x) {
/* Modified Bessel function exp(-x) I0(x), for x >=0.
    We modified from the original bessi0(). Instead of
    I0(x) itself, it returns I0(x) exp(-x).
*/
  double ax, y, ans;

  ax = fabs(x);
  if (ax < 3.75) {
    y = x/3.75;
    y *= y;
    ans = exp(-ax)*(1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492
            + y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2))))));
  } else {
    y = 3.75/ax;
    ans = (1/sqrt(ax))*(0.39894228 + y*(0.1328592e-1
            + y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2
            + y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1
            + y * 0.392377e-2))))))));
  }
  return ans;
}




double ExpBessI0(double r, double r2, double R) {
    double _RR = 1/(R*R);
    double x = 4*r*r2*_RR;
    double y = 2*(r2*r2 + r*r)*_RR;
    double expbess = exp(-y + x)*ModifiedBessI0(x);
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



double Rd_rFGIntegrand(double r2, void * params) {
/* Convolution integrand for either flat or gaussian beams.
    See comments for A_rzFGIntegrand().
    r" in the integration
*/

  MCMLConv * mcmlConv = (MCMLConv *) params;

  double f;
  double * RT_r = mcmlConv->mcmlModel.Rd_r;
  double Rd_at_r2 = RT_rInterp(r2, RT_r, mcmlConv);
  double R = mcmlConv->beam.R;
  double r = mcmlConv->convVar.r;

  if (mcmlConv->beam.type == Beam::FLAT)
    f = Rd_at_r2*ITheta(r, r2, R)*r2;
  else       // Gaussian
    f = Rd_at_r2*ExpBessI0(r, r2, R)*r2;
  return f;
}




double A_rzFGIntegrand(double r2, void * params) {
/* Convolution integrand for either flat or gaussian beams.
    Return the integrand for the convolution integral.
    r2 is the r" in the formula shown in the manual.
    When r2 is in the range of recorded array, interpolation
    is used to evaluate the diffuse reflectance at r2.
    Note that since the last grid elements collect all the
    photon weight that falls beyond the grid system, we should
    avoid using them in the convolution.
    r" in the integration
*/

  MCMLConv * mcmlConv = (MCMLConv *) params;

  double f;
  double A_at_r2 = A_rzInterp(r2, mcmlConv);
  double R = mcmlConv->beam.R;
  double r = mcmlConv->convVar.r;
  Node * tree = mcmlConv->convVar.tree;
  Node * link = Node::SearchNode(tree, r2);

  if (link != nullptr)        // f in tree
    f = link->y;
  else {
    if (mcmlConv->beam.type == Beam::FLAT)
      f = ITheta(r, r2, R);
    else       // Gaussian
      f = ExpBessI0(r, r2, R);
      Node::InsertNode(tree, r2, f);
  }
  f *= A_at_r2*r2;
  return f;
}



double Tt_raFGIntegrand(double r2, void * params) {
/* Convolution integrand for either flat or gaussian beams.
    See comments for A_rzFGIntegrand().
    r" in the integration.
*/

  MCMLConv * mcmlConv = (MCMLConv *) params;

  double f;
  double ** TT_ra = mcmlConv->mcmlModel.Tt_ra;
  double Tt_at_r2 = RT_raInterp(r2, TT_ra, mcmlConv);
  double R = mcmlConv->beam.R;
  double r = mcmlConv->convVar.r;
  Node * tree = mcmlConv->convVar.tree;
  Node * link = Node::SearchNode(tree, r2);
  
  if (link != nullptr)        // f in tree
    f = link->y;
  else {
    if (mcmlConv->beam.type == Beam::FLAT)
      f = ITheta(r, r2, R);
    else	       // Gaussian
      f = ExpBessI0(r, r2, R);
    Node::InsertNode(tree, r2, f);
  }
  f *= Tt_at_r2*r2;
  return f;
}



double Tt_rFGIntegrand(double r2, void * params) {
/* Convolution integrand for either flat or gaussian beams.
    See comments for A_rzFGIntegrand().
    r" in the integration
*/
  MCMLConv * mcmlConv = (MCMLConv *) params;

  double f;
  double * TT_r = mcmlConv->mcmlModel.Tt_r;
  double Tt_at_r2 = RT_rInterp(r2, TT_r, mcmlConv);
  double R = mcmlConv->beam.R;
  double r = mcmlConv->convVar.r;

  if (mcmlConv->beam.type == Beam::FLAT)
    f = Tt_at_r2*ITheta(r, r2, R)*r2;
  else	       // Gaussian
    f = Tt_at_r2*ExpBessI0(r, r2, R)*r2;
  return f;
}



double IntegrateQuad(double (*func) (double x, void * params), 
          double a, double b, MCMLConv * mcmlConv) {

  gsl_integration_workspace * workPtr 
      = gsl_integration_workspace_alloc (2500);

  double absError = 1.0e-5;   // to avoid round-off problems
  double relError = 1.0e-5;   // the result will usually be much better
  double result;              // the result from the integration
  double error;               // the estimated error from the integration

  gsl_function integralFunction;
  void * paramsPtr = mcmlConv;

  integralFunction.function = func;
  integralFunction.params = paramsPtr;

  gsl_integration_qags (&integralFunction, a, b, absError, relError, 2500,
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

