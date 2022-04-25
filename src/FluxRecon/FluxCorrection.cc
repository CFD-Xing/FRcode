///////////////////////////////////////////////////////////////////////////////
/// \file FluxCorrection.cc
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 20 2021
///
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <math.h>
#include <iostream>
using namespace std;

// Includes
#include "FluxCorrection.h"


///////////////////////////////////////////////////////////////////////////////
// FLUX CORRECTION
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// gDG flux correction function
// **************************************************************************//
double gDG(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute DG (Discontinuous Galerkin, aka g1) flux correction function. The 
  // scheme is based on right Radau polynomial of order K, is stable, and has 
  // an accuracy of order 2K-1.
  //
  // NOTE: Only the left boundary flux correction function value is returned. 
  // To obtain the right boundary flux correction function value, the left 
  // boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g       - DG flux correction function


  // VARIABLE DECLARATION
  // ----------------------
  double g;


  // FLUX CORRECTION FUN.
  // ----------------------
  g = RightRadauPolynomial(x, K);


  // RETURN VARIABLE
  // ----------------------
  return g;
}


// **************************************************************************//
// gDG_prime flux correction function derivative
// **************************************************************************//
double gDG_prime(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute DG (Discontinuous Galerkin, aka g1) flux correction function 
  // derivative. The scheme is based on right Radau polynomial of order K, is 
  // stable, and has an accuracy of order 2K-1.
  //
  // NOTE: Only the left boundary flux correction function derivative value is 
  // returned. To obtain the right boundary flux correction function derivative
  // value, the left boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g_prime - DG flux correction function derivative


  // VARIABLE DECLARATION
  // ----------------------
  double g_prime, PolK, PolKm1, PolKm2;


  // FLUX CORRECTION FUN.
  // ----------------------
  if (x==-1.0){
    g_prime = -K*K/2.0;
  }else if (x==1.0){
    g_prime = pow(-1,K)*K/2.0;
  }else{
    PolK   = LegendrePolynomial(x,K);
    PolKm1 = LegendrePolynomial(x,K-1);
    PolKm2 = LegendrePolynomial(x,K-2);
    g_prime = pow(-1,K-1)/(1.0-x*x)/2.0*(K*x*PolK - (K+(K-1.0)*x)*PolKm1 + (K-1.0)*PolKm2);
  }


  // RETURN VARIABLE
  // ----------------------
  return g_prime;
}


// **************************************************************************//
// gSG flux correction function
// **************************************************************************//
double gSG(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute SG (Staggered-Grid) flux correction function. The scheme is based
  // on the K+1 Chebyshev-Lobatto points, is mildly unstable, and has an 
  // accuracy of order K.
  //
  // NOTE: Only the left boundary flux correction function value is returned. 
  // To obtain the right boundary flux correction function value, the left 
  // boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g       - SG flux correction function


  // VARIABLE DECLARATION
  // ----------------------
  double* xi = new double[K+1];
  double g = 1.0;


  // CHEBYSHEV-LOBATTO QUAD.
  // ----------------------
  for (int i=0; i<K+1; i++){
    xi[i] = -cos(i*M_PI/K);  
  }


  // FLUX CORRECTION FUN.
  // ----------------------
  for (int m=0; m<K; m++){
    g*=(x-xi[m+1])/(-1.0-xi[m+1]);
  }
  delete[] xi;


  // RETURN VARIABLE
  // ----------------------
  return g;
}


// **************************************************************************//
// gSG_prime flux correction function derivative
// **************************************************************************//
double gSG_prime(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute SG (Staggered-Grid) flux correction function derivative. The 
  // scheme is based on the K+1 Chebyshev-Lobatto points, is mildly unstable, 
  // and has an accuracy of order K.
  //
  // NOTE: Only the left boundary flux correction function derivative value is 
  // returned. To obtain the right boundary flux correction function derivative
  // value, the left boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g_prime - SG flux correction function derivative


  // VARIABLE DECLARATION
  // ----------------------
  double* xi = new double[K+1];
  double g_prime = 0.0, tmp;


  // CHEBYSHEV-LOBATTO QUAD.
  // ----------------------
  for (int i=0; i<K+1; i++){
    xi[i] = -cos(i*M_PI/K);  
  }


  // FLUX CORRECTION FUN.
  // ----------------------
  for (int m=0; m<K; m++){
    tmp = 1.0;
    for (int n=0; n<K; n++){
      if (n!=m) tmp*=(x-xi[n+1])/(-1.0-xi[n+1]);
    }
    g_prime += tmp/(-1.0-xi[m+1]);
  }
  delete[] xi;


  // RETURN VARIABLE
  // ----------------------
  return g_prime;
}


// **************************************************************************//
// gLo flux correction function
// **************************************************************************//
double gLo(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute Lo (Lobatto) flux correction function. The scheme is based on the
  // K+1 Gauss-Lobatto points, is mildly unstable, and has an accuracy of order
  // K.
  //
  // NOTE: Only the left boundary flux correction function value is returned. 
  // To obtain the right boundary flux correction function value, the left 
  // boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g       - Lo flux correction function


  // VARIABLE DECLARATION
  // ----------------------
  double* xi = new double[K+1];
  double* wi = new double[K+1];
  double g = 1.0;


  // GAUSS-LOBATTO QUAD.
  // ----------------------
  GaussLobattoPoint(K+1, xi, wi);


  // FLUX CORRECTION FUN.
  // ----------------------
  for (int m=0; m<K; m++){
    g*=(x-xi[m+1])/(-1.0-xi[m+1]);
  }
  delete[] xi;
  delete[] wi;


  // RETURN VARIABLE
  // ----------------------
  return g;
}


// **************************************************************************//
// gLo_prime flux correction function derivative
// **************************************************************************//
double gLo_prime(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute Lo (Lobatto) flux correction function derivative. The scheme is 
  // based on the K+1 Gauss-Lobatto points, is mildly unstable, and has an 
  // accuracy of order K.
  //
  // NOTE: Only the left boundary flux correction function derivative value is 
  // returned. To obtain the right boundary flux correction function derivative
  // value, the left boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g_prime - Lo flux correction function derivative


  // VARIABLE DECLARATION
  // ----------------------
  double* xi = new double[K+1];
  double* wi = new double[K+1];
  double g_prime = 0.0, tmp;


  // GAUSS-LOBATTO QUAD.
  // ----------------------
  GaussLobattoPoint(K+1, xi, wi);


  // FLUX CORRECTION FUN.
  // ----------------------
  for (int m=0; m<K; m++){
    tmp = 1.0;
    for (int n=0; n<K; n++){
      if (n!=m) tmp*=(x-xi[n+1])/(-1.0-xi[n+1]);
    }
    g_prime += tmp/(-1.0-xi[m+1]);
  }
  delete[] xi;
  delete[] wi;


  // RETURN VARIABLE
  // ----------------------
  return g_prime;
}


// **************************************************************************//
// gGa flux correction function
// **************************************************************************//
double gGa(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute Ga (Gauss) flux correction function. The scheme is based on the
  // K-1 Gauss-Legendre points, is stable, and has an accuracy of order 2K-2.
  //
  // NOTE: Only the left boundary flux correction function value is returned. 
  // To obtain the right boundary flux correction function value, the left 
  // boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g       - Ga flux correction function


  // VARIABLE DECLARATION
  // ----------------------
  double* xi = new double[K+1];
  double* wi = new double[K+1];
  double g = 1.0;


  // GAUSS-LEGENDRE QUAD.
  // ----------------------
  GaussPoint(K-1, &xi[1], &wi[1]); xi[0] = -1.0; xi[K] = 1.0;


  // FLUX CORRECTION FUN.
  // ----------------------
  for (int m=0; m<K; m++){
    g*=(x-xi[m+1])/(-1.0-xi[m+1]);
  }
  delete[] xi;
  delete[] wi;


  // RETURN VARIABLE
  // ----------------------
  return g;
}


// **************************************************************************//
// gGa_prime flux correction function derivative
// **************************************************************************//
double gGa_prime(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute Ga (Gauss) flux correction function derivative. The scheme is 
  // based on the K-1 Gauss-Legendre points, is stable, and has an accuracy of 
  // order 2K-2.
  //
  // NOTE: Only the left boundary flux correction function derivative value is 
  // returned. To obtain the right boundary flux correction function derivative
  // value, the left boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g_prime - Ga flux correction function


  // VARIABLE DECLARATION
  // ----------------------
  double* xi = new double[K+1];
  double* wi = new double[K+1];
  double g_prime = 0.0, tmp;


  // GAUSS-LEGENDRE QUAD.
  // ----------------------
  GaussPoint(K-1, &xi[1], &wi[1]); xi[0] = -1.0; xi[K] = 1.0;


  // FLUX CORRECTION FUN.
  // ----------------------
  for (int m=0; m<K; m++){
    tmp = 1.0;
    for (int n=0; n<K; n++){
      if (n!=m) tmp*=(x-xi[n+1])/(-1.0-xi[n+1]);
    }
    g_prime += tmp/(-1.0-xi[m+1]);
  }
  delete[] xi;
  delete[] wi;


  // RETURN VARIABLE
  // ----------------------
  return g_prime;
}


// **************************************************************************//
// g2 flux correction function
// **************************************************************************//
double g2(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute g2 flux correction function. The scheme is based on the K Gauss-
  // Lobatto points, is stable, and has an accuracy of order 2K-2. The scheme
  // should be used with K Gauss-Lobatto solution points.
  //
  // NOTE: Only the left boundary flux correction function value is returned. 
  // To obtain the right boundary flux correction function value, the left 
  // boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g       - g2 flux correction function


  // VARIABLE DECLARATION
  // ----------------------
  double g, RadK, RadKm1;


  // FLUX CORRECTION FUNCTION
  // ----------------------
  if (K<2){cout <<"Error, for g2 flux correction function, K>=2!"; exit(-1);}
  RadK   = RightRadauPolynomial(x, K);
  RadKm1 = RightRadauPolynomial(x, K-1);
  g = (K-1.0)/(2.0*K-1.0)*RadK + K/(2.0*K-1.0)*RadKm1;


  // RETURN VARIABLE
  // ----------------------
  return g;
}


// **************************************************************************//
// g2_prime flux correction function derivative
// **************************************************************************//
double g2_prime(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute g2 flux correction function derivative. The scheme is based on the 
  // K Gauss-Lobatto points, is stable, and has an accuracy of order 2K-2. The 
  // scheme should be used with K Gauss-Lobatto solution points.
  //
  // NOTE: Only the left boundary flux correction function derivative value is 
  // returned. To obtain the right boundary flux correction function derivative
  // value, the left boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g_prime - g2 flux correction function derivative


  // VARIABLE DECLARATION
  // ----------------------
  double g_prime, PolK, PolKm2;


  // FLUX CORRECTION FUN.
  // ----------------------
  if (K<2){cout <<"Error, for g2 flux correction function, K>=2!"; exit(-1);}
  if (x==-1.0){
    g_prime = -K*(K-1.0)/2.0;
  }else if (x==1.0){
    g_prime = 0.0;
  }else{
    PolK   = LegendrePolynomial(-x, K);
    PolKm2 = LegendrePolynomial(-x, K-2);
    g_prime = K*(K-1.0)/(2.0*K-1.0)/2.0 * 1.0/(1.0+x) * (PolK - PolKm2);
  }


  // RETURN VARIABLE
  // ----------------------
  return g_prime;
}


// **************************************************************************//
// gKm1 flux correction function
// **************************************************************************//
double gKm1(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute gKm1 flux correction function. The scheme is based on zero z=1 of
  // multiplicity K-1.
  //
  // NOTE: Only the left boundary flux correction function value is returned. 
  // To obtain the right boundary flux correction function value, the left 
  // boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g       - gKm1 flux correction function


  // VARIABLE DECLARATION
  // ----------------------
  double g;


  // FLUX CORRECTION FUN.
  // ----------------------
  g = (K+1.0)*pow((1.0-x)/2.0, K) - K*pow((1.0-x)/2.0, K-1);


  // RETURN VARIABLE
  // ----------------------
  return g;
}


// **************************************************************************//
// gKm1_prime flux correction function derivative
// **************************************************************************//
double gKm1_prime(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute gKm1 flux correction function derivative. The scheme is based on 
  // zero z=1 of multiplicity K-1.
  //
  // NOTE: Only the left boundary flux correction function derivative value is 
  // returned. To obtain the right boundary flux correction function derivative
  // value, the left boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g_prime - gKm1 flux correction function derivative


  // VARIABLE DECLARATION
  // ----------------------
  double g_prime;


  // FLUX CORRECTION FUN.
  // ----------------------
  g_prime = -K/2.0*(K+1.0)*pow((1.0-x)/2.0, K-1) + K/2.0*(K-1.0)*pow((1.0-x)/2.0, K-2);


  // RETURN VARIABLE
  // ----------------------
  return g_prime;
}


// **************************************************************************//
// gK flux correction function
// **************************************************************************//
double gK(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute gK flux correction function. The scheme is based on zero z=1 of
  // multiplicity K.
  //
  // NOTE: Only the left boundary flux correction function value is returned. 
  // To obtain the right boundary flux correction function value, the left 
  // boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g       - gK flux correction function


  // VARIABLE DECLARATION
  // ----------------------
  double g;


  // FLUX CORRECTION FUN.
  // ----------------------
  g = pow((1.0-x)/2.0, K);


  // RETURN VARIABLE
  // ----------------------
  return g;
}


// **************************************************************************//
// gK_prime flux correction function derivative
// **************************************************************************//
double gK_prime(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute gK flux correction function derivative. The scheme is based on zero 
  // z=1 of multiplicity K.
  //
  // NOTE: Only the left boundary flux correction function derivative value is 
  // returned. To obtain the right boundary flux correction function derivative
  // value, the left boundary value need to be reflected at x=0. 


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // g_prime - gK flux correction function derivative


  // VARIABLE DECLARATION
  // ----------------------
  double g_prime;


  // FLUX CORRECTION FUN.
  // ----------------------
  g_prime = -K/2.0*pow((1.0-x)/2.0, K-1);


  // RETURN VARIABLE
  // ----------------------
  return g_prime;
}
