///////////////////////////////////////////////////////////////////////////////
/// \file EulerState1D.cc
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 24 2021
///
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <math.h>
#include <iostream>
#include <stdexcept>
using namespace std;

// Include
#include "EulerState.h"


///////////////////////////////////////////////////////////////////////////////
// EULER 1D 
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// PrimtoCons
// **************************************************************************//
void EulerState1D::PrimtoCons(const double& rho, const double& u, const double& p, double* U){

  // DESCRIPTION
  // ----------------------
  // Convert primitive variable rho, u, p to conservative solution state.


  // INPUTS
  // ----------------------
  // rho    - Density
  // u      - Velocity
  // p      - Pressure


  // OUTPUTS
  // ----------------------
  // U      - Conservative solution state


  // CONSERVATIVE VARIABLE
  // ----------------------
  U[0] = rho;
  U[1] = rho*u;
  U[2] = p/(g-1.0) + 0.5*rho*u*u;
}


// **************************************************************************//
// PrimtoCons
// **************************************************************************//
void EulerState1D::PrimtoCons(const double* W, double* U){

  // DESCRIPTION
  // ----------------------
  // Convert primitive variable W to conservative solution state U.


  // INPUTS
  // ----------------------
  // W      - Primitive solution state


  // OUTPUTS
  // ----------------------
  // U      - Conservative solution state


  // CONSERVATIVE VARIABLE
  // ----------------------
  U[0] = W[0];
  U[1] = W[0]*W[1];
  U[2] = W[2]/(g-1.0) + 0.5*W[0]*W[1]*W[1];
}


// **************************************************************************//
// ConstoPrim
// **************************************************************************//
void EulerState1D::ConstoPrim(const double* U, double& rho, double& u, double& p){

  // DESCRIPTION
  // ----------------------
  // Convert conservative solution state to primitive variable rho, u, p.


  // INPUTS
  // ----------------------
  // U      - Conservative solution state


  // OUTPUTS
  // ----------------------
  // rho    - Density
  // u      - Velocity
  // p      - Pressure


  // PRIMITIVE VARIABLE
  // ----------------------
  // Density
  rho = U[0];
  // Velocity
  u = U[1]/U[0];
  // Pressure
  p = (g-1.0)*(U[2] - 0.5*U[1]*U[1]/U[0]);
}


// **************************************************************************//
// ConstoPrim
// **************************************************************************//
void EulerState1D::ConstoPrim(const double* U, double* W){

  // DESCRIPTION
  // ----------------------
  // Convert conservative solution state to primitive variable rho, u, p.


  // INPUTS
  // ----------------------
  // U      - Conservative solution state


  // OUTPUTS
  // ----------------------
  // W      - Conservative solution state


  // PRIMITIVE VARIABLE
  // ----------------------
  // Density
  W[0] = U[0];
  // Velocity
  W[1] = U[1]/U[0];
  // Pressure
  W[2] = (g-1.0)*(U[2] - 0.5*U[1]*U[1]/U[0]);
}


// **************************************************************************//
// InviscidFlux
// **************************************************************************//
void EulerState1D::InviscidFlux(const double* U, double* F){

  // DESCRIPTION
  // ----------------------
  // Compute the 1D inviscid flux from the conservative solution state.


  // INPUTS
  // ----------------------
  // U      - Conservative solution state


  // OUTPUTS
  // ----------------------
  // F      - Inviscid flux


  // VARIABLE DECLARATION
  // ----------------------
  double p;


  // 1D INVISCID FLUX
  // ----------------------
  // Compute pressure
  p = (g-1.0)*(U[2] - 0.5*U[1]*U[1]/U[0]);

  // Compute flux
  F[0] = U[1];
  F[1] = U[1]*U[1]/U[0] + p;
  F[2] = U[1]/U[0]*(U[2] + p);
}


// **************************************************************************//
// InviscidFlux
// **************************************************************************//
void EulerState1D::InviscidFlux(const double& rho, const double& u, const double& p, double* F){

  // DESCRIPTION
  // ----------------------
  // Compute the 1D inviscid flux from primitive variables rho, u, p.


  // INPUTS
  // ----------------------
  // rho    - Density
  // u      - Velocity
  // p      - Pressure


  // OUTPUTS
  // ----------------------
  // F      - Inviscid flux


  // 1D INVISCID FLUX
  // ----------------------
  F[0] = rho*u;
  F[1] = rho*u*u + p;
  F[2] = u*(g*p/(g-1.0) + 0.5*rho*u*u);
}


// **************************************************************************//
// RoeSolver
// **************************************************************************//
void EulerState1D::RoeSolver(const double* Wl, const double* Wr, double* Fupw){

  // DESCRIPTION
  // ----------------------
  // Roe approximate Riemann solver for 1D Euler.


  // INPUTS
  // ----------------------
  // Wl      - Left solution state
  // Wr      - Right solution state


  // OUTPUTS
  // ----------------------
  // Fupw    - Upwind flux


  // VARIABLE DECLARATION
  // ----------------------
  double rhol, ul, pl, Hl, cl;
  double rhor, ur, pr, Hr, cr;
  double rhoa, ua, pa, Ha, ca;
  double* Fl = new double[3];
  double* Fr = new double[3];
  double* Rc = new double[3];
  double* Lp = new double[3];
  double* alpha = new double[3];
  double* Lambda = new double[3];


  // 1D ROE RIEMANN SOLVER
  // ----------------------
  // Get primitive variables
  rhol = Wl[0]; ul = Wl[1]; pl = Wl[2];
  rhor = Wr[0]; ur = Wr[1]; pr = Wr[2];
  Hl = g/(g-1.0)*pl/rhol + ul*ul/2.0;
  Hr = g/(g-1.0)*pr/rhor + ur*ur/2.0;
  cl = sqrt(g*pl/rhol);
  cr = sqrt(g*pr/rhor);

  // Get Roe average
  rhoa = sqrt(rhol)*sqrt(rhor);
  ua = (sqrt(rhol)*ul+sqrt(rhor)*ur)/(sqrt(rhol)+sqrt(rhor));
  Ha = (sqrt(rhol)*Hl+sqrt(rhor)*Hr)/(sqrt(rhol)+sqrt(rhor));
  pa = (Ha-ua*ua/2.0)*(g-1.0)/g*rhoa;
  ca = sqrt(g*pa/rhoa);
        
  // Compute flux        
  InviscidFlux(rhol,ul,pl,Fl);
  InviscidFlux(rhor,ur,pr,Fr); 
  EigenValue(ua,ca,Lambda);
  EntropyCorrection(ul,cl,ur,cr,Lambda);
  for (int i=0; i<3; i++){
    Fupw[i] = 0.5*(Fl[i]+Fr[i]);
    alpha[i] = 0.0;
    LeftPrimEigenVector(rhoa,ca,i,Lp);
    for (int k=0; k<3; k++){
      alpha[i] += Lp[k]*(Wr[k]-Wl[k]);
    }
  }
  for (int j=0; j<3; j++){
    RightConsEigenVector(ua,ca,j,Rc);
    for (int i=0; i<3; i++){
      Fupw[i] -= 0.5*(fabs(Lambda[j])*alpha[j])*Rc[i];
    }
  }
  
  // Free memory
  delete[] Fl; Fl=NULL;
  delete[] Fr; Fr=NULL;
  delete[] Rc; Rc=NULL;
  delete[] Lp; Lp=NULL;
  delete[] alpha; alpha=NULL;
  delete[] Lambda; Lambda=NULL;
}


// **************************************************************************//
// HLLESolver
// **************************************************************************//
void EulerState1D::HLLESolver(const double* Wl, const double* Wr, double* Fupw){

  // DESCRIPTION
  // ----------------------
  // HLLE approximate Riemann solver for Euler 1D.


  // INPUTS
  // ----------------------
  // Wl      - Left solution state
  // Wr      - Right solution state


  // OUTPUTS
  // ----------------------
  // Fupw    - Upwind flux


  // VARIABLE DECLARATION
  // ----------------------
  double rhol, ul, pl, Hl, cl;
  double rhor, ur, pr, Hr, cr;
  double rhoa, ua, pa, Ha, ca;
  double Lambda_m, Lambda_p;
  double* Fl = new double[3];
  double* Fr = new double[3];
  double* Ul = new double[3];
  double* Ur = new double[3];


  // 1D HLLE RIEMANN SOLVER
  // ----------------------
  // Get primitive variables
  rhol = Wl[0]; ul = Wl[1]; pl = Wl[2];
  rhor = Wr[0]; ur = Wr[1]; pr = Wr[2];
  Hl = g/(g-1.0)*pl/rhol + ul*ul/2.0;
  Hr = g/(g-1.0)*pr/rhor + ur*ur/2.0;
  cl = sqrt(g*pl/rhol);
  cr = sqrt(g*pr/rhor);
        
  // Get Roe average
  rhoa = sqrt(rhol)*sqrt(rhor);
  ua = (sqrt(rhol)*ul+sqrt(rhor)*ur)/(sqrt(rhol)+sqrt(rhor));
  Ha = (sqrt(rhol)*Hl+sqrt(rhor)*Hr)/(sqrt(rhol)+sqrt(rhor));
  pa = (Ha-ua*ua/2.0)*(g-1.0)/g*rhoa;
  ca = sqrt(g*pa/rhoa);

  // Compute flux        
  PrimtoCons(rhol,ul,pl,Ul);
  PrimtoCons(rhor,ur,pr,Ur);
  InviscidFlux(rhol,ul,pl,Fl);
  InviscidFlux(rhor,ur,pr,Fr); 
  Lambda_m = fmin(ul-cl, ua-ca);
  Lambda_p = fmax(ur+cr, ua+ca);
  for (int i=0; i<3; i++){
    if (Lambda_m > 0.0){
      Fupw[i] = Fl[i];
    }else if (Lambda_p < 0.0){
      Fupw[i] = Fr[i];
    }else{
      Fupw[i] = (Lambda_p*Fl[i]-Lambda_m*Fr[i] + 
                Lambda_m*Lambda_p*(Ur[i]-Ul[i])) / (Lambda_p-Lambda_m);
    }
  }

  // Free memory
  delete[] Fl; Fl=NULL;
  delete[] Fr; Fr=NULL;
  delete[] Ul; Ul=NULL;
  delete[] Ur; Ur=NULL;
}


// **************************************************************************//
// EigenValue
// **************************************************************************//
void EulerState1D::EigenValue(const double& u, const double& c, double* Lambda){

  // DESCRIPTION
  // ----------------------
  // Compute the eigenvalue for the 1D euler system.


  // INPUTS
  // ----------------------
  // u       - Velocity
  // c       - Sound speed


  // OUTPUTS
  // ----------------------
  // Lambda  - Eigenvalue


  // EIGEN VALUE
  // ----------------------
  Lambda[0] = u-c;
  Lambda[1] = u;
  Lambda[2] = u+c;
}


// **************************************************************************//
// EntropyCorrection
// **************************************************************************//
void EulerState1D::EntropyCorrection(const double& ul, const double& cl, 
                                     const double& ur, const double& cr,
                                     double* Lambda){

  // DESCRIPTION
  // ----------------------


  // INPUTS
  // ----------------------
  // ul      - Left solution state velocity
  // cl      - Left solution state sound speed
  // ur      - Right solution state velocity
  // cr      - Right solution state sound speed
  // Lambda  - Uncorrected eigenvalues


  // OUTPUTS
  // ----------------------
  // Lambda  - Corrected eigenvalues


  // VARIABLE DECLARATION
  // ----------------------
  double theta, eps;


  // ENTROPY CORRECTION
  // ----------------------
  theta = fmax(0.0, 2.0*((ur-cr) - (ul-cl)));
  if ( fabs(Lambda[0]) < theta ){
    eps = 0.5*(theta + Lambda[0]*Lambda[0]/theta) - abs(Lambda[0]);
  }else{
    eps = 0.0;
  }
  Lambda[0] = Lambda[0] - eps;

  theta = fmax(0.0, 2.0*((ur+cr) - (ul+cl)));
  if ( fabs(Lambda[2]) < theta ){
    eps = 0.5*(theta + Lambda[2]*Lambda[2]/theta) - abs(Lambda[2]);
  }else{
    eps = 0.0;
  }
  Lambda[2] = Lambda[2] + eps;
}


// **************************************************************************//
// LeftPrimEigenVector
// **************************************************************************//
void EulerState1D::LeftPrimEigenVector(const double& rho, const double& c, 
                                       const int& i, double* Lp){

  // DESCRIPTION
  // ----------------------
  // Compute the left eigenvector for primitive variables.


  // INPUTS
  // ----------------------
  // rho     - Density
  // c       - Sound speed
  // i       - Eigenvector index


  // OUTPUTS
  // ----------------------
  // Lp      - Left eigenvector


  // LEFT EIGENVECTOR
  // ----------------------
  switch (i){
    case 0:{
      Lp[0] = 0.0; 
      Lp[1] =-0.5*rho/c;
      Lp[2] = 0.5/(c*c);
      break;}
    case 1:{
      Lp[0] = 1.0;
      Lp[1] = 0.0;
      Lp[2] =-1.0/(c*c);
      break;}
    case 2:{
      Lp[0] = 0.0;
      Lp[1] = 0.5*rho/c;
      Lp[2] = 0.5/(c*c);
      break;}
    default:{
      throw invalid_argument("index out of range"); 
    }
  }
}


// **************************************************************************//
// LeftConsEigenVector
// **************************************************************************//
void EulerState1D::LeftConsEigenVector(const double& u, const double& c, 
                                       const int& i, double* Lc){

  // DESCRIPTION
  // ----------------------
  // Compute the left eigenvector for conservative variables.


  // INPUTS
  // ----------------------
  // u       - Velocity
  // c       - Sound speed
  // i       - Eigenvector index


  // OUTPUTS
  // ----------------------
  // Lc      - Left eigenvector


  // LEFT EIGENVECTOR
  // ----------------------
  switch (i){
    case 0:{
      Lc[0] = ((g-1.0)*u*u/2.0+c*u)/2.0/(c*c);
      Lc[1] = ((1.0-g)*u-c)/2.0/(c*c); 
      Lc[2] = (g-1.0)/2.0/(c*c);
      break;}
    case 1:{
      Lc[0] = (c*c-(g-1.0)*u*u/2.0)/(c*c);
      Lc[1] = (g-1.0)*u/(c*c);
      Lc[2] = (1.0-g)/(c*c);
      break;}
    case 2:{
      Lc[0] = ((g-1.0)*u*u/2.0-c*u)/2.0/(c*c);
      Lc[1] = ((1.0-g)*u+c)/2.0/(c*c);
      Lc[2] = (g-1.0)/2.0/(c*c);
      break;}
    default:{
      throw invalid_argument("index out of range"); 
    }
  }
}


// **************************************************************************//
// RightPrimEigenVector
// **************************************************************************//
void EulerState1D::RightPrimEigenVector(const double& rho, const double& c, 
                                        const int& i, double* Rp){

  // DESCRIPTION
  // ----------------------
  // Compute the right eigenvector for primitive variables.


  // INPUTS
  // ----------------------
  // rho     - Density
  // c       - Sound speed
  // i       - Eigenvector index


  // OUTPUTS
  // ----------------------
  // Rp      - Right eigenvector


  // LEFT EIGENVECTOR
  // ----------------------
  switch (i){
    case 0:{
      Rp[0] = 1.0; 
      Rp[1] =-c/rho;
      Rp[2] = c*c;
      break;}
    case 1:{
      Rp[0] = 1.0;
      Rp[1] = 0.0;
      Rp[2] = 0.0;
      break;}
    case 2:{
      Rp[0] = 1.0;
      Rp[1] = c/rho;
      Rp[2] = c*c;
      break;}
    default:{
      throw invalid_argument("index out of range"); 
    }
  }
}


// **************************************************************************//
// RightConsEigenVector1D
// **************************************************************************//
void EulerState1D::RightConsEigenVector(const double& u, const double& c, 
                                        const int& i, double* Rc){

  // DESCRIPTION
  // ----------------------
  // Compute the right eigenvector for conservative variables.


  // INPUTS
  // ----------------------
  // u       - Velocity
  // c       - Sound speed
  // rho     - Density
  // p       - Pressure
  // i       - Eigenvector index


  // OUTPUTS
  // ----------------------
  // Rc      - Right eigenvector


  // RIGHT EIGENVECTOR
  // ----------------------
  switch (i){
    case 0:{
      Rc[0] = 1.0;
      Rc[1] = u-c;
      Rc[2] = c*c/(g-1.0)+u*u/2.0-u*c;
      break;}
    case 1:{
      Rc[0] = 1.0;
      Rc[1] = u;
      Rc[2] = u*u/2.0;
      break;}
    case 2:{
      Rc[0] = 1.0; 
      Rc[1] = u+c;
      Rc[2] = c*c/(g-1.0)+u*u/2.0+u*c;
      break;}
    default:{
      throw invalid_argument("index out of range"); 
    }
  }
}


// **************************************************************************//
// FluxJacobian
// **************************************************************************//
void EulerState1D::FluxJacobian(const double* W, double** dFdU){
    
  // DESCRIPTION
  // ----------------------
  // Compute the flux Jacobian


  // INPUTS
  // ----------------------
  // W      - Primitive solution state


  // OUTPUTS
  // ----------------------
  // dFdU   - Flux Jacobian

    
  // FLUX JACOBIAN
  // ----------------------
  dFdU[0][0] = 0.0;
  dFdU[0][1] = 1.0;
  dFdU[0][2] = 0.0;

  dFdU[1][0] = (g-3.0)/2.0*W[1]*W[1];
  dFdU[1][1] = (3.0-g)*W[1];
  dFdU[1][2] = (g-1.0);
  
  dFdU[2][0] = (g-2.0)/2.0*W[1]*W[1]*W[1] - g/(g-1.0)*W[1]*W[2]/W[0];
  dFdU[2][1] = g/(g-1.0)*W[2]/W[0] + (3.0-2.0*g)/2.0*W[1]*W[1];
  dFdU[2][2] = g*W[1];  
}
