///////////////////////////////////////////////////////////////////////////////
/// \file FRelement.cc
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 21 2021
///
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <stdexcept>
using namespace std;

// Includes
#include "FRelement.h"


///////////////////////////////////////////////////////////////////////////////
// FR ELEMENT
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// setElement
// **************************************************************************//
void FRelement1D::setElement(){

  // DESCRIPTION
  // ----------------------
  // Compute the solution points, eta, the left solution propagation vector, 
  // Lvec, the right solution propagation vector, Rvec, the differential matrix,
  // Dmat, the left flux correction derivatives, gl_prime, and the right flux
  // correction derivatives, gr_prime.


  // SET ELEMENTS
  // ----------------------
  // Solution points
  setSolutionPoints();

  // Left and right solution propagation vector
  LagrangeBasis1D(K, eta, -1.0, Lvec);
  LagrangeBasis1D(K, eta,  1.0, Rvec);

  // Differential matrix
  for (int i=0; i<K; i++){
    LagrangeBasisDerivative1D(K, eta, eta[i], &Dmat[i][0]);
  }

  // Flux correction function derivative
  setFluxCorrectionDerivative();
}


// **************************************************************************//
// setSolutionPoints
// **************************************************************************//
void FRelement1D::setSolutionPoints(){

  // DESCRIPTION
  // ----------------------
  // Compute the solution points specified by the solution type, Ktype, and the
  // number of solution point, K.
  

  // SOLUTION POINTS
  // ----------------------
  switch (Ktype){
    case 1:{ // Equidistant without boundary
      for (int i=0; i<K; i++){
        eta[i] = -1.0 + ((i+1)-0.5)*2.0/K;
      }
      break;}
    case 2:{ // Equidistant with boundary
      for (int i=0; i<K; i++){
        eta[i] = -1.0 + i*2.0/(K-1.0);
      }
      break;}
    case 3:{ // Gauss-Legendre points
      double* dump=new double[K];
      GaussPoint(K, eta, dump);
      delete[] dump;
      break;}
    case 4:{ // Gauss-Lobatto points
      double* dump=new double[K];
      GaussLobattoPoint(K, eta, dump);
      delete[] dump;
      break;}
    case 5:{ // Chebyshev-Gauss points
      for (int i=0; i<K; i++){
        eta[i] = -cos(((i+1)-0.5)*M_PI/K);
      }
      break;}
    case 6:{ // Chebyshev-Lobatto points
      for (int i=0; i<K; i++){
        eta[i] = -cos(i*M_PI/(K-1));
      }
      break;}
    default:{
      throw invalid_argument("Unsupported solution points"); 
    }
  }
}


// **************************************************************************//
// setFluxCorrectionDerivative
// **************************************************************************//
void FRelement1D::setFluxCorrectionDerivative(){

  // DESCRIPTION
  // ----------------------
  // Compute the left and right flux correction function derivatives. The left 
  // boundary flux correction function derivatives are computed directly for the 
  // flux correction type, FRtype, and for each solution point, eta[i]. The right 
  // boundary flux correction derivatives are computed by reflecting the values
  // of the left boundary flux correction derivatives.
  

  // FLUX CORR. DERIVATIVES
  // ----------------------
  // Left boundary flux correction function derivatives
  switch (FRtype){
    case 1:{ // DG (Discontinuous Galerkin)
      for (int i=0; i<K; i++){
        gl_prime[i] = gDG_prime(eta[i], K);
      }
      strcpy(FRname, "DG");
      break;}
    case 2:{ // g2
      for (int i=0; i<K; i++){
        gl_prime[i] = g2_prime(eta[i], K);
      }
      strcpy(FRname, "g2");
      break;}
    case 3:{ // Ga
      for (int i=0; i<K; i++){
        gl_prime[i] = gGa_prime(eta[i], K);
      }
      strcpy(FRname, "Ga");
      break;}
    case 4:{ // SG (Staggered Grid)
      for (int i=0; i<K; i++){
        gl_prime[i] = gSG_prime(eta[i], K);
      }
      strcpy(FRname, "SG");
      break;}
    case 5:{ // Lo
      for (int i=0; i<K; i++){
        gl_prime[i] = gLo_prime(eta[i], K);
      }
      strcpy(FRname, "Lo");
      break;}
    case 6:{ // K-1
      for (int i=0; i<K; i++){
        gl_prime[i] = gKm1_prime(eta[i], K);
      }
      break;}
    case 7:{ // K
      for (int i=0; i<K; i++){
        gl_prime[i] = gK_prime(eta[i], K);
      }
      break;}
    default:{
      throw invalid_argument("Unsupported flux correction function"); 
    }
  }

  // Right boundary flux correction function derivatives
  for (int i=0; i<K; i++){
    gr_prime[i] = -gl_prime[K-(i+1)];
  }
}
