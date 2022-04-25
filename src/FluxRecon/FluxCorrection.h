///////////////////////////////////////////////////////////////////////////////
/// \file FluxCorrection.h
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 20 2021
///
///////////////////////////////////////////////////////////////////////////////
#ifndef FR_CORRECTION_H
#define FR_CORRECTION_H


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <math.h>
#include <iostream>
using namespace std;

// Includes
#include "Quadrature.h"
#include "Polynomial.h"
#define M_PI       3.14159265358979323846


///////////////////////////////////////////////////////////////////////////////
// FLUX CORRECTION
///////////////////////////////////////////////////////////////////////////////
double gDG(const double& x, const int& K);
double gDG_prime(const double& x, const int& K);
double gSG(const double& x, const int& K);
double gSG_prime(const double& x, const int& K);
double gLo(const double& x, const int& K);
double gLo_prime(const double& x, const int& K);
double gGa(const double& x, const int& K);
double gGa_prime(const double& x, const int& K);
double g2(const double& x, const int& K);
double g2_prime(const double& x, const int& K);
double gKm1(const double& x, const int& K);
double gKm1_prime(const double& x, const int& K);
double gK(const double& x, const int& K);
double gK_prime(const double& x, const int& K);
#endif
