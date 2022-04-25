///////////////////////////////////////////////////////////////////////////////
/// \file Lagrange.h
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 20 2021
///
///////////////////////////////////////////////////////////////////////////////
#ifndef LAGRANGE_H
#define LAGRANGE_H


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <iostream>
using namespace std;


///////////////////////////////////////////////////////////////////////////////
// LAGRANGE POLYNOMIAL INTERPOLATION
///////////////////////////////////////////////////////////////////////////////
void LagrangeBasis1D(const int& nx,
                     const double* xi,
                     const double& xp,
                     double* li);

void LagrangeBasisDerivative1D(const int& nx,
                               const double* xi,
                               const double& xp,
                               double* dlidx);

double LagrangeInterpol1D(const int& nx,
                          const double* xi,
                          const double* f, 
                          const double& xp);

void LagrangeBasis2D(const int& nx,
                     const int& ny,
                     const double* xi, 
                     const double* yj,
                     const double& xp,
                     const double& yp,
                     double** lij);

void LagrangeBasisDerivative2D(const int& nx,
                               const int& ny,
                               const double* xi, 
                               const double* yj,
                               const double& xp,
                               const double& yp,
                               double** dlijdx,
                               double** dlijdy);

double LagrangeInterpol2D(const int& nx,
                          const int& ny,
                          const double* xi, 
                          const double* yj,
                          const double** f,
                          const double& xp,
                          const double& yp);

void LagrangeBasis3D(const int& nx,
                     const int& ny,
                     const int& nz,
                     const double* xi, 
                     const double* yj,
                     const double* zk,
                     const double& xp,
                     const double& yp,
                     const double& zp,
                     double*** lijk);

void LagrangeBasisDerivative3D(const int& nx,
                               const int& ny,
                               const int& nz,
                               const double* xi, 
                               const double* yj,
                               const double* zk,
                               const double& xp,
                               const double& yp,
                               const double& zp,
                               double*** dlijkdx,
                               double*** dlijkdy,
                               double*** dlijkdz);

double LagrangeInterpol3D(const int& nx,
                          const int& ny,
                          const int& nz,
                          const double* xi, 
                          const double* yj,
                          const double* zk,
                          const double*** f,
                          const double& xp,
                          const double& yp,
                          const double& zp);

#endif
