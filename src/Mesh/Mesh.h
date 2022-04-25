///////////////////////////////////////////////////////////////////////////////
/// \file Mesh.h
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 30 2021
///
///////////////////////////////////////////////////////////////////////////////
#ifndef MESH_H
#define MESH_H


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <iostream>
using namespace std;


///////////////////////////////////////////////////////////////////////////////
// CLASS DEFINITIONS
///////////////////////////////////////////////////////////////////////////////
class Mesh1D{
  public:
  int ng;
  int nx;
  double* xc;
  double* dx;
  // Default constructor
  Mesh1D(){ Nullify(); }
  // Constructor
  Mesh1D(int ng, int nx){
    this->ng = ng;
    this->nx = nx;
    Allocate();
  }
  // Default destructor
  ~Mesh1D(){ Deallocate(); Nullify(); }
  // Some functions
  void setMesh1D(int ng, int nx){
    Deallocate();
    this->ng = ng;
    this->nx = nx;
    Allocate();
  }
  void setUniformMesh(const double& x0, const double& x1);
  // Allocate
  void Allocate(){
    xc = new double[nx+2*ng];
    dx = new double[nx+2*ng];
  }
  // Deallocate
  void Deallocate(){
    if (xc!=NULL) delete[] xc;
    if (dx!=NULL) delete[] dx;
  }
  // Nullify
  void Nullify(){
    this->ng = -1;
    this->nx = -1;
    xc = NULL;
    dx = NULL;
  }
};


///////////////////////////////////////////////////////////////////////////////
// CLASS DEFINITIONS
///////////////////////////////////////////////////////////////////////////////
class Mesh2D{
  public:
  int ng;
  int nx;
  int ny;
  double** xc;
  double** yc;
  double** dx;
  double** dy;
  // Default constructor
  Mesh2D(){ Nullify(); }
  // Constructor
  Mesh2D(int ng, int nx, int ny){
    this->ng = ng;
    this->nx = nx;
    this->ny = ny;
    Allocate();
  }
  // Default constructor
  ~Mesh2D(){ Deallocate(); Nullify(); }
  // Some functions
  void setMesh2D(int ng, int nx, int ny){
    Deallocate();
    this->ng = ng;
    this->nx = nx;
    this->ny = ny;
    Allocate();
  }
  void setUniformMesh(const double& x0, const double& x1,
                      const double& y0, const double& y1);
  // Allocate
  void Allocate(){
    xc = new double*[nx+2*ng];
    yc = new double*[nx+2*ng];
    dx = new double*[nx+2*ng];
    dy = new double*[nx+2*ng];
    for (int i=0; i<nx+2*ng; i++){
      xc[i] = new double[ny+2*ng];
      yc[i] = new double[ny+2*ng];
      dx[i] = new double[ny+2*ng];
      dy[i] = new double[ny+2*ng];
    }
  }
  // Deallocate
  void Deallocate(){
    for (int i=0; i<nx+2*ng; i++){
      if (xc!=NULL) delete[] xc[i];
      if (yc!=NULL) delete[] yc[i];
      if (dx!=NULL) delete[] dx[i];
      if (dy!=NULL) delete[] dy[i];
    }
    if (xc!=NULL) delete[] xc;
    if (yc!=NULL) delete[] yc;
    if (dx!=NULL) delete[] dx;
    if (dy!=NULL) delete[] dy;
  }
  // Nullify
  void Nullify(){
    this->ng = -1;
    this->nx = -1;
    this->ny = -1;
    xc = NULL;
    yc = NULL;
    dx = NULL;
    dy = NULL;
  }
};


///////////////////////////////////////////////////////////////////////////////
// CLASS DEFINITIONS
///////////////////////////////////////////////////////////////////////////////
class Mesh3D{
  public:
  int ng;
  int nx;
  int ny;
  int nz;
  double*** xc;
  double*** yc;
  double*** zc;
  double*** dx;
  double*** dy;
  double*** dz;
  // Default constructor
  Mesh3D(){ Nullify(); }
  // Constructor
  Mesh3D(int ng, int nx, int ny, int nz){
    this->ng = ng;
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
    Allocate();
  }
  // Default destructor
  ~Mesh3D(){ Deallocate(); Nullify(); }
  // Some functions
  void setMesh3D(int ng, int nx, int ny, int nz){
    Deallocate();
    this->ng = ng;
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
    Allocate();
  }
  void setUniformMesh(const double& x0, const double& x1,
                      const double& y0, const double& y1,
                      const double& z0, const double& z1);
  // Allocate
  void Allocate(){
    xc = new double**[nx+2*ng];
    yc = new double**[nx+2*ng];
    zc = new double**[nx+2*ng];
    dx = new double**[nx+2*ng];
    dy = new double**[nx+2*ng];
    dz = new double**[nx+2*ng];
    for (int i=0; i<nx+2*ng; i++){
      xc[i] = new double*[ny+2*ng];
      yc[i] = new double*[ny+2*ng];
      zc[i] = new double*[ny+2*ng];
      dx[i] = new double*[ny+2*ng];
      dy[i] = new double*[ny+2*ng];
      dz[i] = new double*[ny+2*ng];
      for (int j=0; j<ny+2*ng; j++){
        xc[i][j] = new double[nz+2*ng];
        yc[i][j] = new double[nz+2*ng];
        zc[i][j] = new double[nz+2*ng];
        dx[i][j] = new double[nz+2*ng];
        dy[i][j] = new double[nz+2*ng];
        dz[i][j] = new double[nz+2*ng];
      }
    }
  }
  // Deallocate
  void Deallocate(){
    for (int i=0; i<nx+2*ng; i++){
      for (int j=0; j<ny+2*ng; j++){
        if (xc!=NULL) delete[] xc[i][j];
        if (yc!=NULL) delete[] yc[i][j];
        if (zc!=NULL) delete[] zc[i][j];
        if (dx!=NULL) delete[] dx[i][j];
        if (dy!=NULL) delete[] dy[i][j];
        if (dz!=NULL) delete[] dz[i][j];
      }
      if (xc!=NULL) delete[] xc[i];
      if (yc!=NULL) delete[] yc[i];
      if (zc!=NULL) delete[] zc[i];
      if (dx!=NULL) delete[] dx[i];
      if (dy!=NULL) delete[] dy[i];
      if (dz!=NULL) delete[] dz[i];
    }
    if (xc!=NULL) delete[] xc;
    if (yc!=NULL) delete[] yc;
    if (zc!=NULL) delete[] zc;
    if (dx!=NULL) delete[] dx;
    if (dy!=NULL) delete[] dy;
    if (dz!=NULL) delete[] dz;
  }
  // Nullify
  void Nullify(){
    this->ng = -1;
    this->nx = -1;
    this->ny = -1;
    this->nz = -1;
    xc = NULL;
    yc = NULL;
    zc = NULL;
    dx = NULL;
    dy = NULL;
    dz = NULL;
  }
};
#endif
