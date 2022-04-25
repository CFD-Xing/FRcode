///////////////////////////////////////////////////////////////////////////////
/// \file Polynomial.h
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 20 2021
///
///////////////////////////////////////////////////////////////////////////////
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <math.h>
#include <iostream>
using namespace std;


///////////////////////////////////////////////////////////////////////////////
// ORTHOGONAL POLYNOMIAL
///////////////////////////////////////////////////////////////////////////////
double LegendrePolynomial(const double& x, const int& K);
double LeftRadauPolynomial(const double& x, const int& K);
double RightRadauPolynomial(const double& x, const int& K);
double LobattoPolynomial(const double& x, const int& K);


// **************************************************************************//
// LegendrePolynomial
// **************************************************************************//
inline double LegendrePolynomial(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute Legendre orthogonal polynomial of order K at solution point x.


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // PolK    - Legendre polynomial of order K at solution point x


  // VARIABLE DECLARATION
  // ----------------------
  double PolK, PolKm1, PolKm2;


  // ORTHOGONAL POLYNOMIALS
  // ----------------------
  if (K==-1){
    PolK = 0.0;
  }else if (K==0){
    PolK = 1.0;
  }else{
    PolKm2 = 0.0;
    PolKm1 = 1.0;
    for (int k=1; k<=K; k++){
      PolK   = (2.0*k-1.0)/k*x*PolKm1 - (k-1.0)/k*PolKm2;
      PolKm2 = PolKm1;
      PolKm1 = PolK;
    }
  }


  // RETURN VARIABLE
  // ----------------------
  return PolK;
}


// **************************************************************************//
// LeftRadauPolynomial
// **************************************************************************//
inline double LeftRadauPolynomial(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute left Radau orthogonal polynomial of order K at solution point x.


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // RadK    - (Left) Radau polynomial of order K at solution point x


  // VARIABLE DECLARATION
  // ----------------------
  double RadK;


  // ORTHOGONAL POLYNOMIALS
  // ----------------------
  if (K<1){cout <<"Error, for Radau polynomial, K>=1!"; exit(-1);}
  RadK = 1./2.*(LegendrePolynomial(x, K) + LegendrePolynomial(x, K-1));


  // RETURN VARIABLE
  // ----------------------
  return RadK;
}


// **************************************************************************//
// RightRadauPolynomial
// **************************************************************************//
inline double RightRadauPolynomial(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute right Radau orthogonal polynomial of order K at solution point x.


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // RadK    - (Right) Radau polynomial of order K at solution point x


  // VARIABLE DECLARATION
  // ----------------------
  double RadK;


  // ORTHOGONAL POLYNOMIALS
  // ----------------------
  if (K<1){cout <<"Error, for Radau polynomial, K>=1!"; exit(-1);}
  RadK = pow(-1.,K)/2.*(LegendrePolynomial(x, K) - LegendrePolynomial(x, K-1));


  // RETURN VARIABLE
  // ----------------------
  return RadK;
}


// **************************************************************************//
// LobattoPolynomial
// **************************************************************************//
inline double LobattoPolynomial(const double& x, const int& K){

  // DESCRIPTION
  // ----------------------
  // Compute Lobatto orthogonal polynomial of order K at solution point x.


  // INPUTS
  // ----------------------
  // x       - Solution point
  // K       - Polynomial order


  // OUTPUTS
  // ----------------------
  // LoK    - Lobatto polynomial of order K at solution point x


  // VARIABLE DECLARATION
  // ----------------------
  double LoK;


  // ORTHOGONAL POLYNOMIALS
  // ----------------------
  if (K<2){cout <<"Error, for Lobatto polynomial, K>=2!"; exit(-1);}
  LoK = 2.*(LegendrePolynomial(x, K) - LegendrePolynomial(x, K-2));


  // RETURN VARIABLE
  // ----------------------
  return LoK;
}
#endif
