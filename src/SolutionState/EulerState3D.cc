///////////////////////////////////////////////////////////////////////////////
/// \file EulerState3D.cc
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
// EULER 3D
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// PrimtoCons3D
// **************************************************************************//
void EulerState3D::PrimtoCons(const double& rho, const double& u, const double& v, const double& w, const double& p, double* U){

  // DESCRIPTION
  // ----------------------
  // Convert primitive variable rho, u, v, w, p to conservative solution state.


  // INPUTS
  // ----------------------
  // rho    - Density
  // u      - Velocity
  // v      - Velocity
  // w      - Velocity
  // p      - Pressure


  // OUTPUTS
  // ----------------------
  // U      - Conservative solution state


  // CONSERVATIVE VARIABLE
  // ----------------------
  U[0] = rho;
  U[1] = rho*u;
  U[2] = rho*v;
  U[3] = rho*w;
  U[4] = p/(g-1.0) + 0.5*rho*(u*u + v*v + w*w);
}


// **************************************************************************//
// PrimtoCons3D
// **************************************************************************//
void EulerState3D::PrimtoCons(const double* W, double* U){

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
  U[2] = W[0]*W[2];
  U[3] = W[0]*W[3];
  U[4] = W[4]/(g-1.0) + 0.5*W[0]*(W[1]*W[1] + W[2]*W[2] + W[3]*W[3]);
}


// **************************************************************************//
// ConstoPrim3D
// **************************************************************************//
void EulerState3D::ConstoPrim(const double* U, double& rho, double& u, double& v, double& w, double& p){

  // DESCRIPTION
  // ----------------------
  // Convert conservative solution state to primitive variable rho, u, v, w, p.


  // INPUTS
  // ----------------------
  // U      - Conservative solution state


  // OUTPUTS
  // ----------------------
  // rho    - Density
  // u      - Velocity
  // v      - Velocity
  // w      - Velocity
  // p      - Pressure


  // PRIMITIVE VARIABLE
  // ----------------------
  // Density
  rho = U[0];
  // Velocity
  u = U[1]/U[0];
  // Velocity
  v = U[2]/U[0];
  // Velocity
  w = U[3]/U[0];
  // Pressure
  p = (g-1.0)*(U[4] - 0.5*(U[1]*U[1]+U[2]*U[2]+U[3]*U[3])/U[0]);
}


// **************************************************************************//
// ConstoPrim3D
// **************************************************************************//
void EulerState3D::ConstoPrim(const double* U, double* W){

  // DESCRIPTION
  // ----------------------
  // Convert conservative solution state to primitive variable rho, u, v, w, p.


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
  // Velocity
  W[2] = U[2]/U[0];
  // Velocity
  W[3] = U[3]/U[0];
  // Pressure
  W[4] = (g-1.0)*(U[4] - 0.5*(U[1]*U[1]+U[2]*U[2]+U[3]*U[3])/U[0]);
}


// **************************************************************************//
// InviscidFlux3D
// **************************************************************************//
void EulerState3D::InviscidFlux(const double* U, double* F){

  // DESCRIPTION
  // ----------------------
  // Compute the 3D inviscid flux from the conservative solution state.


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
  p = (g-1.0)*(U[4] - 0.5*(U[1]*U[1]+U[2]*U[2]+U[3]*U[3])/U[0]);

  // Compute flux
  F[0] = U[1];
  F[1] = U[1]*U[1]/U[0] + p;
  F[2] = U[1]*U[2]/U[0];
  F[3] = U[1]*U[3]/U[0];
  F[4] = U[1]/U[0]*(U[4] + p);
}


// **************************************************************************//
// InviscidFlux3D
// **************************************************************************//
void EulerState3D::InviscidFlux(const double& rho, const double& u, const double& v, const double& w, const double& p, double* F){

  // DESCRIPTION
  // ----------------------
  // Compute the 3D inviscid flux from primitive variables rho, u, v, w, p.


  // INPUTS
  // ----------------------
  // rho    - Density
  // u      - Velocity
  // v      - Velocity
  // w      - Velocity
  // p      - Pressure


  // OUTPUTS
  // ----------------------
  // F      - Inviscid flux


  // 1D INVISCID FLUX
  // ----------------------
  F[0] = rho*u;
  F[1] = rho*u*u + p;
  F[2] = rho*u*v;
  F[3] = rho*u*w;
  F[4] = u*(g*p/(g-1.0) + 0.5*rho*(u*u+v*v+w*w));
}


// **************************************************************************//
// RoeSolver3D
// **************************************************************************//
void EulerState3D::RoeSolver(const double* Wl, const double* Wr, double* Fupw){

  // DESCRIPTION
  // ----------------------
  // Roe approximate Riemann solver for 3D Euler.


  // INPUTS
  // ----------------------
  // Wl      - Left solution state
  // Wr      - Right solution state


  // OUTPUTS
  // ----------------------
  // Fupw    - Upwind flux


  // VARIABLE DECLARATION
  // ----------------------
  double rhol, ul, vl, wl, pl, Hl, cl;
  double rhor, ur, vr, wr, pr, Hr, cr;
  double rhoa, ua, va, wa, pa, Ha, ca;


  // 1D ROE RIEMANN SOLVER
  // ----------------------
  // Get primitive variables
  rhol = Wl[0]; ul = Wl[1]; vl = Wl[2]; wl = Wl[3]; pl = Wl[4];
  rhor = Wr[0]; ur = Wr[1]; vr = Wr[2]; wr = Wr[3]; pr = Wr[4];
  Hl = g/(g-1.0)*pl/rhol + (ul*ul+vl*vl+wl*wl)/2.0;
  Hr = g/(g-1.0)*pr/rhor + (ur*ur+vr*vr+wr*wr)/2.0;
  cl = sqrt(g*pl/rhol);
  cr = sqrt(g*pr/rhor);

  // Get Roe average
  rhoa = sqrt(rhol)*sqrt(rhor);
  ua = (sqrt(rhol)*ul+sqrt(rhor)*ur)/(sqrt(rhol)+sqrt(rhor));
  va = (sqrt(rhol)*vl+sqrt(rhor)*vr)/(sqrt(rhol)+sqrt(rhor));
  wa = (sqrt(rhol)*wl+sqrt(rhor)*wr)/(sqrt(rhol)+sqrt(rhor));
  Ha = (sqrt(rhol)*Hl+sqrt(rhor)*Hr)/(sqrt(rhol)+sqrt(rhor));
  pa = (Ha-(ua*ua+va*va+wa*wa)/2.0)*(g-1.0)/g*rhoa;
  ca = sqrt(g*pa/rhoa);
        
  // Compute flux        
  InviscidFlux(rhol,ul,vl,wl,pl,Fl);
  InviscidFlux(rhor,ur,vr,wr,pr,Fr); 
  EigenValue(ua,ca,Lambda);
  EntropyCorrection(ul,cl,ur,cr,Lambda);
  for (int i=0; i<5; i++){
    Fupw[i] = 0.5*(Fl[i]+Fr[i]);
    alpha[i] = 0.0;
    LeftPrimEigenVector(rhoa,ca,i,Lp);
    for (int k=0; k<5; k++){
      alpha[i] += Lp[k]*(Wr[k]-Wl[k]);
    }
  }
  for (int j=0; j<5; j++){
    RightConsEigenVector(ua,va,wa,ca,j,Rc);
    for (int i=0; i<5; i++){
      Fupw[i] -= 0.5*(fabs(Lambda[j])*alpha[j])*Rc[i];
    }
  }
}


// **************************************************************************//
// HLLESolver3D
// **************************************************************************//
void EulerState3D::HLLESolver(const double* Wl, const double* Wr, double* Fupw){

  // DESCRIPTION
  // ----------------------
  // HLLE approximate Riemann solver for Euler 3D.


  // INPUTS
  // ----------------------
  // Wl      - Left solution state
  // Wr      - Right solution state


  // OUTPUTS
  // ----------------------
  // Fupw    - Upwind flux


  // VARIABLE DECLARATION
  // ----------------------
  double rhol, ul, vl, wl, pl, Hl, cl;
  double rhor, ur, vr, wr, pr, Hr, cr;
  double rhoa, ua, va, wa, pa, Ha, ca;
  double Lambda_m, Lambda_p;


  // 1D HLLE RIEMANN SOLVER
  // ----------------------
  // Get primitive variables
  rhol = Wl[0]; ul = Wl[1]; vl = Wl[2]; wl = Wl[3]; pl = Wl[4];
  rhor = Wr[0]; ur = Wr[1]; vr = Wr[2]; wr = Wr[3]; pr = Wr[4];
  Hl = g/(g-1.0)*pl/rhol + (ul*ul+vl*vl+wl*wl)/2.0;
  Hr = g/(g-1.0)*pr/rhor + (ur*ur+vr*vr+wr*wr)/2.0;
  cl = sqrt(g*pl/rhol);
  cr = sqrt(g*pr/rhor);
        
  // Get Roe average
  rhoa = sqrt(rhol)*sqrt(rhor);
  ua = (sqrt(rhol)*ul+sqrt(rhor)*ur)/(sqrt(rhol)+sqrt(rhor));
  va = (sqrt(rhol)*vl+sqrt(rhor)*vr)/(sqrt(rhol)+sqrt(rhor));
  wa = (sqrt(rhol)*wl+sqrt(rhor)*wr)/(sqrt(rhol)+sqrt(rhor));
  Ha = (sqrt(rhol)*Hl+sqrt(rhor)*Hr)/(sqrt(rhol)+sqrt(rhor));
  pa = (Ha-(ua*ua+va*va+wa*wa)/2.0)*(g-1.0)/g*rhoa;
  ca = sqrt(g*pa/rhoa);

  // Compute flux        
  PrimtoCons(rhol,ul,vl,wl,pl,Ul);
  PrimtoCons(rhor,ur,vr,wr,pr,Ur);
  InviscidFlux(rhol,ul,vl,wl,pl,Fl);
  InviscidFlux(rhor,ur,vr,wr,pr,Fr); 
  Lambda_m = min(ul-cl, ua-ca);
  Lambda_p = max(ur+cr, ua+ca);
  for (int i=0; i<5; i++){
    Fupw[i] = (Lambda_p*Fl[i]-Lambda_m*Fr[i] + 
                Lambda_m*Lambda_p*(Ur[i]-Ul[i])) / (Lambda_p-Lambda_m);
  }
}


// **************************************************************************//
// EigenValue3D
// **************************************************************************//
void EulerState3D::EigenValue(const double& u, const double& c, double* Lambda){

  // DESCRIPTION
  // ----------------------
  // Compute the eigenvalue for the 3D euler system.


  // INPUTS
  // ----------------------
  // u       - Velocity
  // c       - Sound speed


  // OUTPUTS
  // ----------------------
  // Lambda  - Eigenvalue


  // VARIABLE DECLARATION
  // ----------------------


  // EIGEN VALUE
  // ----------------------
  Lambda[0] = u-c;
  Lambda[1] = u;
  Lambda[2] = u;
  Lambda[3] = u;
  Lambda[4] = u+c;
}


// **************************************************************************//
// EntropyCorrection3D
// **************************************************************************//
void EulerState3D::EntropyCorrection(const double& ul, const double& cl, 
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
  if ( fabs(Lambda[4]) < theta ){
    eps = 0.5*(theta + Lambda[4]*Lambda[4]/theta) - abs(Lambda[4]);
  }else{
    eps = 0.0;
  }
  Lambda[4] = Lambda[4] + eps;
}


// **************************************************************************//
// LeftPrimEigenVector3D
// **************************************************************************//
void EulerState3D::LeftPrimEigenVector(const double& rho, const double& c, 
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
      Lp[2] = 0.0; 
      Lp[3] = 0.0; 
      Lp[4] = 0.5/(c*c);
      break;}
    case 1:{
      Lp[0] = 1.0;
      Lp[1] = 0.0;
      Lp[2] = 0.0;
      Lp[3] = 0.0;
      Lp[4] =-1.0/(c*c);
      break;}
    case 2:{
      Lp[0] = 0.0;
      Lp[1] = 0.0;
      Lp[2] = 1.0;
      Lp[3] = 0.0;
      Lp[4] = 0.0;
      break;}
    case 3:{
      Lp[0] = 0.0;
      Lp[1] = 0.0;
      Lp[2] = 0.0;
      Lp[3] = 1.0;
      Lp[4] = 0.0;
      break;}
    case 4:{
      Lp[0] = 0.0;
      Lp[1] = 0.5*rho/c;
      Lp[2] = 0.0;
      Lp[3] = 0.0;
      Lp[4] = 0.5/(c*c);
      break;}
    default:{
      throw invalid_argument("index out of range"); 
    }
  }
}


// **************************************************************************//
// LeftConsEigenVector3D
// **************************************************************************//
void EulerState3D::LeftConsEigenVector(const double& u, const double& v, const double& w, const double& c, 
                                       const int& i, double* Lc){

  // DESCRIPTION
  // ----------------------
  // Compute the left eigenvector for conservative variables.


  // INPUTS
  // ----------------------
  // u       - Velocity
  // v       - Velocity
  // w       - Velocity
  // c       - Sound speed
  // i       - Eigenvector index


  // OUTPUTS
  // ----------------------
  // Lc      - Left eigenvector


  // LEFT EIGENVECTOR
  // ----------------------
  switch (i){
    case 0:{
      Lc[0] = ((g-1.0)*(u*u+v*v+w*w)/2.0+c*u)/2.0/(c*c);
      Lc[1] = ((1.0-g)*u-c)/2.0/(c*c); 
      Lc[2] = ((1.0-g)*v)/2.0/(c*c); 
      Lc[3] = ((1.0-g)*w)/2.0/(c*c); 
      Lc[4] = (g-1.0)/2.0/(c*c);
      break;}
    case 1:{
      Lc[0] = (c*c-(g-1.0)*(u*u+v*v+w*w)/2.0)/(c*c);
      Lc[1] = (g-1.0)*u/(c*c);
      Lc[2] = (g-1.0)*v/(c*c);
      Lc[3] = (g-1.0)*w/(c*c);
      Lc[4] = (1.0-g)/(c*c);
      break;}
    case 2:{
      Lc[0] = -v;
      Lc[1] = 0.0;
      Lc[2] = 1.0;
      Lc[3] = 0.0;
      Lc[4] = 0.0;
      break;}
    case 3:{
      Lc[0] = -w;
      Lc[1] = 0.0;
      Lc[2] = 0.0;
      Lc[3] = 1.0;
      Lc[4] = 0.0;
      break;}
    case 4:{
      Lc[0] = ((g-1.0)*(u*u+v*v+w*w)/2.0-c*u)/2.0/(c*c);
      Lc[1] = ((1.0-g)*u+c)/2.0/(c*c);
      Lc[2] = ((1.0-g)*v)/2.0/(c*c);
      Lc[3] = ((1.0-g)*w)/2.0/(c*c);
      Lc[4] = (g-1.0)/2.0/(c*c);
      break;}
    default:{
      throw invalid_argument("index out of range"); 
    }
  }
}


// **************************************************************************//
// RightPrimEigenVector3D
// **************************************************************************//
void EulerState3D::RightPrimEigenVector(const double& rho, const double& c, 
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
      Rp[2] = 0.0; 
      Rp[3] = 0.0; 
      Rp[4] = c*c;
      break;}
    case 1:{
      Rp[0] = 1.0;
      Rp[1] = 0.0;
      Rp[2] = 0.0;
      Rp[3] = 0.0;
      Rp[4] = 0.0;
      break;}
    case 2:{
      Rp[0] = 0.0;
      Rp[1] = 0.0;
      Rp[2] = 1.0;
      Rp[3] = 0.0;
      Rp[4] = 0.0;
      break;}
    case 3:{
      Rp[0] = 0.0;
      Rp[1] = 0.0;
      Rp[2] = 0.0;
      Rp[3] = 1.0;
      Rp[4] = 0.0;
      break;}
    case 4:{
      Rp[0] = 1.0;
      Rp[1] = c/rho;
      Rp[2] = 0.0;
      Rp[3] = 0.0;
      Rp[4] = c*c;
      break;}
    default:{
      throw invalid_argument("index out of range"); 
    }
  }
}


// **************************************************************************//
// RightConsEigenVector3D
// **************************************************************************//
void EulerState3D::RightConsEigenVector(const double& u, const double& v, const double& w, const double& c, 
                                        const int& i, double* Rc){

  // DESCRIPTION
  // ----------------------
  // Compute the right eigenvector for conservative variables.


  // INPUTS
  // ----------------------
  // u       - Velocity
  // v       - Velocity
  // w       - Velocity
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
      Rc[2] = v;
      Rc[3] = w;
      Rc[4] = c*c/(g-1.0)+(u*u+v*v+w*w)/2.0-u*c;
      break;}
    case 1:{
      Rc[0] = 1.0;
      Rc[1] = u;
      Rc[2] = v;
      Rc[3] = w;
      Rc[4] = (u*u+v*v+w*w)/2.0;
      break;}
    case 2:{
      Rc[0] = 0.0;
      Rc[1] = 0.0;
      Rc[2] = 1.0;
      Rc[3] = 0.0;
      Rc[4] = v;
      break;}
    case 3:{
      Rc[0] = 0.0;
      Rc[1] = 0.0;
      Rc[2] = 0.0;
      Rc[3] = 1.0;
      Rc[4] = w;
      break;}
    case 4:{
      Rc[0] = 1.0; 
      Rc[1] = u+c;
      Rc[2] = v;
      Rc[3] = w;
      Rc[4] = c*c/(g-1.0)+(u*u+v*v+w*w)/2.0+u*c;
      break;}
    default:{
      throw invalid_argument("index out of range"); 
    }
  }
}


// **************************************************************************//
// FluxJacobian
// **************************************************************************//
void EulerState3D::FluxJacobian(const double* W, double** dFdU){

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
  dFdU[0][3] = 0.0;
  dFdU[0][4] = 0.0;
  
  dFdU[1][0] = (g-3.0)/2.0*W[1]*W[1] + (g-1.0)/2.0*(W[2]*W[2]+W[3]*W[3]);
  dFdU[1][1] = (3.0-g)*W[1];
  dFdU[1][2] = (1.0-g)*W[2];
  dFdU[1][3] = (1.0-g)*W[3];
  dFdU[1][4] = (g-1.0);

  dFdU[2][0] =-W[1]*W[2];
  dFdU[2][1] = W[2];
  dFdU[2][2] = W[1];
  dFdU[2][3] = 0.0;
  dFdU[2][4] = 0.0;
  
  dFdU[3][0] =-W[1]*W[3];
  dFdU[3][1] = W[3];
  dFdU[3][2] = 0.0;
  dFdU[3][3] = W[1];
  dFdU[3][4] = 0.0;
  
  dFdU[4][0] = (g-2.0)/2.0*W[1]*(W[1]*W[1]+W[2]*W[2]+W[3]*W[3]) - g/(g-1.0)*W[1]*W[4]/W[0];
  dFdU[4][1] = g/(g-1.0)*W[4]/W[0] + ((3.0-2.0*g)*W[1]*W[1]+W[2]*W[2]+W[3]*W[3])/2.0;
  dFdU[4][2] = (1.0-g)*W[1]*W[2];
  dFdU[4][3] = (1.0-g)*W[1]*W[3];
  dFdU[4][4] = g*W[1];  
}
