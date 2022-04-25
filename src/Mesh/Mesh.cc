///////////////////////////////////////////////////////////////////////////////
/// \file Mesh.cc
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 30 2021
///
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <iostream>
using namespace std;

// Include
#include "Mesh.h"


///////////////////////////////////////////////////////////////////////////////
// CLASS DEFINITIONS
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// setUniformMesh (1D)
// **************************************************************************//
void Mesh1D::setUniformMesh(const double& x0, const double& x1){

  // DESCRIPTION
  // ----------------------
  // Generate a 1D uniform mesh.

    
  // VARIABLE DECLARATION
  // ----------------------
  double xi, xf;
    
    
  // 1D UNIFORM MESH
  // ----------------------
  // Cell-center coordinates
  for (int i=0; i<nx+2*ng; i++){
    xc[i] = x0 + (i-ng+0.5)*(x1-x0)/nx;
  }
  
  // Cell size
  for (int i=1; i<nx+2*ng-1; i++){
    xi = (xc[i]   + xc[i-1])/2.0;
    xf = (xc[i+1] + xc[i])/2.0;
    dx[i] = xf-xi;
  }
  dx[0]         = 2.0*(xc[0]-x0);
  dx[nx+2*ng-1] = 2.0*(xf-xc[nx+2*ng-1]);
}


// **************************************************************************//
// setUniformMesh (2D)
// **************************************************************************//
void Mesh2D::setUniformMesh(const double& x0, const double& x1,
                            const double& y0, const double& y1){

  // DESCRIPTION
  // ----------------------
  // Generate a 2D uniform mesh.

    
  // VARIABLE DECLARATION
  // ----------------------
  double xi, xf, yi, yf;
    

  // 2D UNIFORM MESH
  // ----------------------
  // Cell-center coordinates
  for (int i=0; i<nx+2*ng; i++){
    for (int j=0; j<ny+2*ng; j++){
      xc[i][j] = x0 + (i-ng+0.5)*(x1-x0)/nx;
      yc[i][j] = y0 + (j-ng+0.5)*(y1-y0)/ny;
    }
  }
  
  // Cell size
  for (int j=0; j<ny+2*ng; j++){
    for (int i=1; i<nx+2*ng-1; i++){
      xi = (xc[i][j]   + xc[i-1][j])/2.0;
      xf = (xc[i+1][j] + xc[i][j])/2.0;
      dx[i][j] = xf-xi;
    }
    dx[0][j]         = 2.0*(xc[0][j]-x0);
    dx[nx+2*ng-1][j] = 2.0*(xf-xc[nx+2*ng-1][j]);
  }
  
  for (int i=0; i<nx+2*ng; i++){
    for (int j=1; j<ny+2*ng-1; j++){
      yi = (yc[i][j]   + yc[i][j-1])/2.0;
      yf = (yc[i][j+1] + yc[i][j])/2.0;
      dy[i][j] = yf-yi;
    }
    dy[i][0]         = 2.0*(yc[i][0]-y0);
    dy[i][ny+2*ng-1] = 2.0*(yf-yc[i][ny+2*ng-1]);
  }
}


// **************************************************************************//
// setUniformMesh (3D)
// **************************************************************************//
void Mesh3D::setUniformMesh(const double& x0, const double& x1,
                            const double& y0, const double& y1,
                            const double& z0, const double& z1){

  // DESCRIPTION
  // ----------------------
  // Generate a 3D uniform mesh.

    
  // VARIABLE DECLARATION
  // ----------------------
  double xi, xf, yi, yf, zi, zf;
    
  
  
  // 3D UNIFORM MESH
  // ----------------------
  for (int i=0; i<nx+2*ng; i++){
    for (int j=0; j<ny+2*ng; j++){
      for (int k=0; k<nz+2*ng; k++){
        xc[i][j][k] = x0 + (i-ng+0.5)*(x1-x0)/nx;
        yc[i][j][k] = y0 + (j-ng+0.5)*(y1-y0)/ny;
        zc[i][j][k] = z0 + (k-ng+0.5)*(z1-z0)/nz;
      }
    }
  }
  
  // Cell size
  for (int k=0; k<nz+2*ng; k++){
    for (int j=0; j<ny+2*ng; j++){
      for (int i=1; i<nx+2*ng-1; i++){
        xi = (xc[i][j][k]   + xc[i-1][j][k])/2.0;
        xf = (xc[i+1][j][k] + xc[i][j][k])/2.0;
        dx[i][j][k] = xf-xi;
      }
      dx[0][j][k]         = 2.0*(xc[0][j][k]-x0);
      dx[nx+2*ng-1][j][k] = 2.0*(xf-xc[nx+2*ng-1][j][k]);
    } 
  }
  
  for (int k=0; k<nz+2*ng; k++){
    for (int i=0; i<nx+2*ng; i++){
      for (int j=1; j<ny+2*ng-1; j++){
        yi = (yc[i][j][k]   + yc[i][j-1][k])/2.0;
        yf = (yc[i][j+1][k] + yc[i][j][k])/2.0;
        dy[i][j][k] = yf-yi;
      }
      dy[i][0][k]         = 2.0*(yc[i][0][k]-y0);
      dy[i][ny+2*ng-1][k] = 2.0*(yf-yc[i][ny+2*ng-1][k]);
    } 
  }
  
  for (int i=0; i<nx+2*ng; i++){
    for (int j=0; j<ny+2*ng; j++){
      for (int k=1; k<nz+2*ng-1; k++){
        zi = (zc[i][j][k]   + zc[i][j][k-1])/2.0;
        zf = (zc[i][j][k+1] + zc[i][j][k])/2.0;
        dz[i][j][k] = zf-zi;
      }
      dz[i][j][0]         = 2.0*(zc[i][j][0]-z0);
      dz[i][j][nz+2*ng-1] = 2.0*(zf-zc[i][j][nz+2*ng-1]);
    } 
  }
}
