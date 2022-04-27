///////////////////////////////////////////////////////////////////////////////
/// \file FlowSolver.h
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date January 14, 2022
///
///////////////////////////////////////////////////////////////////////////////
#ifndef FLOWSOLVER_H
#define FLOWSOLVER_H


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <math.h>
#include <iostream>
#include <stdexcept>
using namespace std;

// Include
#include "../Mesh/Mesh.h"
#include "../SolutionState/EulerState.h"
#include "../FluxRecon/FRelement.h"


///////////////////////////////////////////////////////////////////////////////
// CLASS DEFINITIONS
///////////////////////////////////////////////////////////////////////////////
class Solver{
  public:
  int nx, ny, nz;
  int nvar;
  int K;
  int SolverType;
  int EquationType;
  int LimiterType;
  int BCtype;
  int Ktype;
  int FRtype;
  int nResidual=6;
  char LimiterName[256];
  FRelement1D FR;
  double wavespeed=20.0;
};


///////////////////////////////////////////////////////////////////////////////
// CLASS DEFINITIONS
///////////////////////////////////////////////////////////////////////////////
class Solver1D: public Solver{
  public:
  Mesh1D mesh;
  EulerState1D Euler1D;
  double*  dt;
  double*** dUdt;
  double** phi;
  double** dWdx;
  double** W;
  double** U;
  double** Wn;
  double** Un;
  double** W0;
  double** U0;
  double** F;
  double** Fupw;
  double** Fl;
  double** Fr;
  double** Ul;
  double** Ur;
  double** Wl;
  double** Wr;
  Solver1D(){ Nullify(); }
  Solver1D(int K, int nx, int SolverType, int EquationType, int LimiterType, int Ktype, int FRtype){
    this->K = (SolverType == 2) ? K : 1;
    this->nx = nx/this->K;
    this->SolverType = SolverType;
    this->EquationType = EquationType;
    this->LimiterType = LimiterType;
    this->Ktype = Ktype;
    this->FRtype = FRtype;
    Allocate();
  }
  ~Solver1D(){ Deallocate(); Nullify(); }
  void setSolver1D(int K, int nx, int SolverType, int EquationType, int LimiterType, int Ktype, int FRtype){
    Deallocate();
    this->K = (SolverType == 2) ? K : 1;
    this->nx = nx/this->K;;
    this->SolverType = SolverType;
    this->EquationType = EquationType;
    this->LimiterType = LimiterType;
    this->Ktype = Ktype;
    this->FRtype = FRtype;
    Allocate();
  }
  // Allocate
  void Allocate(){
    if (EquationType==0)      nvar = 1;
    else if (EquationType==1) nvar = 1;
    else if (EquationType==2) nvar = 3;
    mesh.setMesh1D(1,nx);
    FR.setFRelement1D(K, Ktype, FRtype);
    dt = new double[(nx+2)*K];
    dUdt = new double**[(nx+2)*K];
    phi = new double*[(nx+2)*K];
    dWdx = new double*[(nx+2)*K];
    W = new double*[(nx+2)*K];
    U = new double*[(nx+2)*K];
    Wn = new double*[(nx+2)*K];
    Un = new double*[(nx+2)*K];
    W0 = new double*[(nx+2)*K];
    U0 = new double*[(nx+2)*K];
    F = new double*[(nx+2)*K];
    Fupw = new double*[nx+1];
    Fl = new double*[nx+2];
    Fr = new double*[nx+2];
    Ul = new double*[nx+2];
    Ur = new double*[nx+2];
    Wl = new double*[nx+2];
    Wr = new double*[nx+2];
    for (int i=0; i<(nx+2)*K; i++){
      dUdt[i] = new double*[nResidual];
      phi[i] = new double[3];
      dWdx[i] = new double[3];
      W[i] = new double[3];
      U[i] = new double[3];
      Wn[i] = new double[3];
      Un[i] = new double[3];
      W0[i] = new double[3];
      U0[i] = new double[3];
      F[i] = new double[3];
      for (int n=0; n<nResidual; n++){
        dUdt[i][n] = new double[3];
      }
    }
    for (int i=0; i<nx+1; i++){
      Fupw[i] = new double[3];
    }
    for (int i=0; i<nx+2; i++){
      Fl[i] = new double[3];
      Fr[i] = new double[3];
      Ul[i] = new double[3];
      Ur[i] = new double[3];
      Wl[i] = new double[3];
      Wr[i] = new double[3];
    }
    switch (LimiterType){
      case 0:
        strcpy(LimiterName, "FirstOrder");
        break;
      case 1:
        strcpy(LimiterName, "NoLimiter");
        break;
      case 2:
        strcpy(LimiterName, "Barth-Jespersen");
        break;
      case 3:
        strcpy(LimiterName, "Venkatakrishnan");
        break;
    }
  }
  // Deallocate
  void Deallocate(){
    mesh.Deallocate();
    FR.Deallocate();
    for (int i=0; i<(nx+2)*K; i++){
      for (int n=0; n<nResidual; n++){
        if (dUdt!=NULL) delete[] dUdt[i][n];
      }
      if (dUdt!=NULL) delete[] dUdt[i];
      if (phi!=NULL) delete[] phi[i];
      if (dWdx!=NULL) delete[] dWdx[i];
      if (W!=NULL) delete[] W[i];
      if (U!=NULL) delete[] U[i];
      if (Wn!=NULL) delete[] Wn[i];
      if (Un!=NULL) delete[] Un[i];
      if (W0!=NULL) delete[] W0[i];
      if (U0!=NULL) delete[] U0[i];
      if (F!=NULL) delete[] F[i];
    }
    for (int i=0; i<nx+1; i++){
      if (Fupw!=NULL) delete[] Fupw[i];
    }
    for (int i=0; i<nx+2; i++){
      if (Fl!=NULL) delete[] Fl[i];
      if (Fr!=NULL) delete[] Fr[i];
      if (Ul!=NULL) delete[] Ul[i];
      if (Ur!=NULL) delete[] Ur[i];
      if (Wl!=NULL) delete[] Wl[i];
      if (Wr!=NULL) delete[] Wr[i];
    }
    if (dt!=NULL) delete[] dt;
    if (dUdt!=NULL) delete[] dUdt;
    if (phi!=NULL) delete[] phi;
    if (dWdx!=NULL) delete[] dWdx;
    if (W!=NULL) delete[] W;
    if (U!=NULL) delete[] U;
    if (Wn!=NULL) delete[] Wn;
    if (Un!=NULL) delete[] Un;
    if (W0!=NULL) delete[] W0;
    if (U0!=NULL) delete[] U0;
    if (F!=NULL) delete[] F;
    if (Fupw!=NULL) delete[] Fupw;
    if (Fl!=NULL) delete[] Fl;
    if (Fr!=NULL) delete[] Fr;
    if (Ul!=NULL) delete[] Ul;
    if (Ur!=NULL) delete[] Ur;
    if (Wl!=NULL) delete[] Wl;
    if (Wr!=NULL) delete[] Wr;
  }
  // Nullify
  void Nullify(){
    this->K = -1;
    this->nx = -1;
    this->nvar = -1;
    this->SolverType = -1;
    this->EquationType = -1;
    this->LimiterType = -1;
    this->BCtype = -1;
    this->Ktype = -1;
    this->FRtype = -1;
    mesh.Nullify();
    FR.Nullify();
    dt = NULL;
    dUdt = NULL;
    phi = NULL;
    dWdx = NULL;
    W = NULL;
    U = NULL;
    Wn = NULL;
    Un = NULL;
    W0 = NULL;
    U0 = NULL;
    F = NULL;
    Fupw = NULL;
    Fl = NULL;
    Fr = NULL;
    Ul = NULL;
    Ur = NULL;
    Wl = NULL;
    Wr = NULL;
  }
  void setMesh(const int& option, 
               const double& x0, const double& x1);
  double TimeStep(const double& CFL);
  void ApplyBCs(void);
  void UpdatePrimitive(void);
  void EnforcePositivity(void);
  void Gradient(void);
  void Limiter(const int& i);
  void Limiter(void);
  void ResetResidual(const int& index);
  void Residual(const int& index);
  void FVMResidual(const int& index);
  void FRResidual(const int& index);
  void CommonFluxXdir(double* commonflux, int i);
  void LinearAdvectionInitialCondition(const int& option, const double& ncycle, double& t);
  void LinearAdvectionExactSolution(const int& option, const double& t);
  void LinearAdvectionOutput(const int& option, const double& CFL, const int& nStage, const double& t);
  void EulerInitialCondition(const int& option, double& t);
  void EulerExactSolution(const int& option);
  void EulerOutput(const int& option, const double& CFL, const int& nStage);
  void RungeKutta(const int& nstage);
};


///////////////////////////////////////////////////////////////////////////////
// CLASS DEFINITIONS
///////////////////////////////////////////////////////////////////////////////
class Solver2D: public Solver{
  public:
  Mesh2D mesh;
  EulerState2D Euler2D;
  double**  dt;
  double**** dUdt;
  double*** phi;
  double*** dWdx;
  double*** dWdy;
  double*** W;
  double*** U;
  double*** Wn;
  double*** Un;
  double*** W0;
  double*** U0;
  double*** F;
  double*** Fupw;
  double*** Fl;
  double*** Fr;
  double*** Ul;
  double*** Ur;
  double*** Wl;
  double*** Wr;
  Solver2D(){ Nullify(); }
  Solver2D(int K, int nx, int ny, int SolverType, int EquationType, int LimiterType, int Ktype, int FRtype){
    this->K = (SolverType == 2) ? K : 1;
    this->nx = nx/this->K;
    this->ny = ny/this->K;
    this->SolverType = SolverType;
    this->EquationType = EquationType;
    this->LimiterType = LimiterType;
    this->Ktype = Ktype;
    this->FRtype = FRtype;
    Allocate();
  }
  ~Solver2D(){ Deallocate(); Nullify(); }
  void setSolver2D(int K, int nx, int ny, int SolverType, int EquationType, int LimiterType, int Ktype, int FRtype){
    Deallocate();
    this->K = (SolverType == 2) ? K : 1;
    this->nx = nx/this->K;
    this->ny = ny/this->K;
    this->SolverType = SolverType;
    this->EquationType = EquationType;
    this->LimiterType = LimiterType;
    this->Ktype = Ktype;
    this->FRtype = FRtype;
    Allocate();
  }
  // Allocate
  void Allocate(){
    if (EquationType==0)      nvar = 1;
    else if (EquationType==1) nvar = 1;
    else if (EquationType==2) nvar = 4;
    mesh.setMesh2D(1,nx,ny);
    FR.setFRelement1D(K, Ktype, FRtype);
    dt = new double*[(nx+2)*K];
    dUdt = new double***[(nx+2)*K];
    phi = new double**[(nx+2)*K];
    dWdx = new double**[(nx+2)*K];
    dWdy = new double**[(nx+2)*K];
    W = new double**[(nx+2)*K];
    U = new double**[(nx+2)*K];
    Wn = new double**[(nx+2)*K];
    Un = new double**[(nx+2)*K];
    W0 = new double**[(nx+2)*K];
    U0 = new double**[(nx+2)*K];
    F = new double**[(nx+2)*K];
    Fupw = new double**[nx+1];
    Fl = new double**[nx+2];
    Fr = new double**[nx+2];
    Ul = new double**[nx+2];
    Ur = new double**[nx+2];
    Wl = new double**[nx+2];
    Wr = new double**[nx+2];
    for (int i=0; i<(nx+2)*K; i++){
      dt[i] = new double[(ny+2)*K];
      dUdt[i] = new double**[(ny+2)*K];
      phi[i] = new double*[(ny+2)*K];
      dWdx[i] = new double*[(ny+2)*K];
      dWdy[i] = new double*[(ny+2)*K];
      W[i] = new double*[(ny+2)*K];
      U[i] = new double*[(ny+2)*K];
      Wn[i] = new double*[(ny+2)*K];
      Un[i] = new double*[(ny+2)*K];
      W0[i] = new double*[(ny+2)*K];
      U0[i] = new double*[(ny+2)*K];
      F[i] = new double*[(ny+2)*K];
      for (int j=0; j<(ny+2)*K; j++){
        dUdt[i][j] = new double*[nResidual];
        phi[i][j] = new double[4];
        dWdx[i][j] = new double[4];
        dWdy[i][j] = new double[4];
        W[i][j] = new double[4];
        U[i][j] = new double[4];
        Wn[i][j] = new double[4];
        Un[i][j] = new double[4];
        W0[i][j] = new double[4];
        U0[i][j] = new double[4];
        F[i][j] = new double[4];
        for (int n=0; n<nResidual; n++){
          dUdt[i][j][n] = new double[4];
        }
      }
    }
    for (int i=0; i<nx+1; i++){
      Fupw[i] = new double*[ny+1];
      for (int j=0; j<ny+1; j++){
        Fupw[i][j] = new double[4];
      }
    }
    for (int i=0; i<nx+2; i++){
      Fl[i] = new double*[ny+2];
      Fr[i] = new double*[ny+2];
      Ul[i] = new double*[ny+2];
      Ur[i] = new double*[ny+2];
      Wl[i] = new double*[ny+2];
      Wr[i] = new double*[ny+2];
      for (int j=0; j<ny+2; j++){
        Fl[i][j] = new double[4];
        Fr[i][j] = new double[4];
        Ul[i][j] = new double[4];
        Ur[i][j] = new double[4];
        Wl[i][j] = new double[4];
        Wr[i][j] = new double[4];
      }
    }
    switch (LimiterType){
      case 0:
        strcpy(LimiterName, "FirstOrder");
        break;
      case 1:
        strcpy(LimiterName, "NoLimiter");
        break;
      case 2:
        strcpy(LimiterName, "Barth-Jespersen");
        break;
      case 3:
        strcpy(LimiterName, "Venkatakrishnan");
        break;
    }
  }
  // Deallocate
  void Deallocate(){
    mesh.Deallocate();
    FR.Deallocate();
    for (int i=0; i<(nx+2)*K; i++){
      for (int j=0; j<(ny+2)*K; j++){
        for (int n=0; n<nResidual; n++){
          if (dUdt!=NULL) delete[] dUdt[i][j][n];
        }
        if (dUdt!=NULL) delete[] dUdt[i][j];
        if (phi!=NULL) delete[] phi[i][j];
        if (dWdx!=NULL) delete[] dWdx[i][j];
        if (dWdy!=NULL) delete[] dWdy[i][j];
        if (W!=NULL) delete[] W[i][j];
        if (U!=NULL) delete[] U[i][j];
        if (Wn!=NULL) delete[] Wn[i][j];
        if (Un!=NULL) delete[] Un[i][j];
        if (W0!=NULL) delete[] W0[i][j];
        if (U0!=NULL) delete[] U0[i][j];
        if (F!=NULL) delete[] F[i][j];
      }
      if (dt!=NULL) delete[] dt[i];
      if (dUdt!=NULL) delete[] dUdt[i];
      if (phi!=NULL) delete[] phi[i];
      if (dWdx!=NULL) delete[] dWdx[i];
      if (dWdy!=NULL) delete[] dWdy[i];
      if (W!=NULL) delete[] W[i];
      if (U!=NULL) delete[] U[i];
      if (Wn!=NULL) delete[] Wn[i];
      if (Un!=NULL) delete[] Un[i];
      if (W0!=NULL) delete[] W0[i];
      if (U0!=NULL) delete[] U0[i];
      if (F!=NULL) delete[] F[i];
    }
    for (int i=0; i<nx+1; i++){
      for (int j=0; j<ny+1; j++){
        if (Fupw!=NULL) delete[] Fupw[i][j];
      }
      if (Fupw!=NULL) delete[] Fupw[i];
    }
    for (int i=0; i<nx+2; i++){
      for (int j=0; j<ny+2; j++){
        if (Fl!=NULL) delete[] Fl[i][j];
        if (Fr!=NULL) delete[] Fr[i][j];
        if (Ul!=NULL) delete[] Ul[i][j];
        if (Ur!=NULL) delete[] Ur[i][j];
        if (Wl!=NULL) delete[] Wl[i][j];
        if (Wr!=NULL) delete[] Wr[i][j];
      }
      if (Fl!=NULL) delete[] Fl[i];
      if (Fr!=NULL) delete[] Fr[i];
      if (Ul!=NULL) delete[] Ul[i];
      if (Ur!=NULL) delete[] Ur[i];
      if (Wl!=NULL) delete[] Wl[i];
      if (Wr!=NULL) delete[] Wr[i];
    }
    if (dt!=NULL) delete[] dt;
    if (dUdt!=NULL) delete[] dUdt;
    if (phi!=NULL) delete[] phi;
    if (dWdx!=NULL) delete[] dWdx;
    if (dWdy!=NULL) delete[] dWdy;
    if (W!=NULL) delete[] W;
    if (U!=NULL) delete[] U;
    if (Wn!=NULL) delete[] Wn;
    if (Un!=NULL) delete[] Un;
    if (W0!=NULL) delete[] W0;
    if (U0!=NULL) delete[] U0;
    if (F!=NULL) delete[] F;
    if (Fupw!=NULL) delete[] Fupw;
    if (Fl!=NULL) delete[] Fl;
    if (Fr!=NULL) delete[] Fr;
    if (Ul!=NULL) delete[] Ul;
    if (Ur!=NULL) delete[] Ur;
    if (Wl!=NULL) delete[] Wl;
    if (Wr!=NULL) delete[] Wr;
  }
  // Nullify
  void Nullify(){
    this->K = -1;
    this->nx = -1;
    this->ny = -1;
    this->nvar = -1;
    this->SolverType = -1;
    this->EquationType = -1;
    this->LimiterType = -1;
    this->BCtype = -1;
    this->Ktype = -1;
    this->FRtype = -1;
    mesh.Nullify();
    FR.Nullify();
    dt = NULL;
    dUdt = NULL;
    phi = NULL;
    dWdx = NULL;
    dWdy = NULL;
    W = NULL;
    U = NULL;
    Wn = NULL;
    Un = NULL;
    W0 = NULL;
    U0 = NULL;
    F = NULL;
    Fupw = NULL;
    Fl = NULL;
    Fr = NULL;
    Ul = NULL;
    Ur = NULL;
    Wl = NULL;
    Wr = NULL;
  }
  void setMesh(const int& option, 
               const double& x0, const double& x1, 
               const double& y0, const double& y1);
  double TimeStep(const double& CFL);
  void ApplyBCs(void);
  void UpdatePrimitive(void);
  void EnforcePositivity(void);
  void Gradient(void);
  void Limiter(const int& i, const int& j);
  void Limiter(void);
  void ResetResidual(const int& index);
  void Residual(const int& index);
  void FVMResidual(const int& index);
  void FRResidual(const int& index);
  void CommonFluxXdir(double* commonflux, int i, int j, int jp);
  void CommonFluxYdir(double* commonflux, int i, int j, int ip);
  void LinearAdvectionInitialCondition(const int& index, const double& ncycle, double& t);
  void LinearAdvectionExactSolution(const int& index, const double& t);
  void LinearAdvectionOutput(const int& option, const double& CFL, const int& nStage, const double& t);
  void EulerInitialCondition(const int& option, double& t);
  void EulerExactSolution(const int& option);
  void EulerOutput(const int& option, const double& CFL, const int& nStage);
  void RungeKutta(const int& nstage);
};


///////////////////////////////////////////////////////////////////////////////
// CLASS DEFINITIONS
///////////////////////////////////////////////////////////////////////////////
class Solver3D: public Solver{
  public:
  Mesh3D mesh;
  EulerState3D Euler3D;
  double***  dt;
  double***** dUdt;
  double**** phi;
  double**** dWdx;
  double**** dWdy;
  double**** dWdz;
  double**** W;
  double**** U;
  double**** Wn;
  double**** Un;
  double**** W0;
  double**** U0;
  double**** F;
  double**** Fupw;
  double**** Fl;
  double**** Fr;
  double**** Ul;
  double**** Ur;
  double**** Wl;
  double**** Wr;
  Solver3D(){ Nullify(); }
  Solver3D(int K, int nx, int ny, int nz, int SolverType, int EquationType, int LimiterType, int Ktype, int FRtype){
    this->K = (SolverType == 2) ? K : 1;
    this->nx = nx/this->K;
    this->ny = ny/this->K;
    this->nz = nz/this->K;
    this->nvar = nvar;
    this->SolverType = SolverType;
    this->EquationType = EquationType;
    this->LimiterType = LimiterType;
    this->Ktype = Ktype;
    this->FRtype = FRtype;
    Allocate();
  }
  ~Solver3D(){ Deallocate(); Nullify(); }
  void setSolver3D(int K, int nx, int ny, int nz, int SolverType, int EquationType, int LimiterType, int Ktype, int FRtype){
    Deallocate();
    this->K = (SolverType == 2) ? K : 1;
    this->nx = nx/this->K;
    this->ny = ny/this->K;
    this->nz = nz/this->K;
    this->SolverType = SolverType;
    this->EquationType = EquationType;
    this->LimiterType = LimiterType;
    this->Ktype = Ktype;
    this->FRtype = FRtype;
    Allocate();
  }
  // Allocate
  void Allocate(){
    if (EquationType==0)      nvar = 1;
    else if (EquationType==1) nvar = 1;
    else if (EquationType==2) nvar = 5;
    mesh.setMesh3D(1,nx,ny,nz);
    FR.setFRelement1D(K, Ktype, FRtype);
    dt = new double**[(nx+2)*K];
    dUdt = new double****[(nx+2)*K];
    phi = new double***[(nx+2)*K];
    dWdx = new double***[(nx+2)*K];
    dWdy = new double***[(nx+2)*K];
    dWdz = new double***[(nx+2)*K];
    W = new double***[(nx+2)*K];
    U = new double***[(nx+2)*K];
    Wn = new double***[(nx+2)*K];
    Un = new double***[(nx+2)*K];
    W0 = new double***[(nx+2)*K];
    U0 = new double***[(nx+2)*K];
    F = new double***[(nx+2)*K];
    Fupw = new double***[nx+1];
    Fl = new double***[nx+2];
    Fr = new double***[nx+2];
    Ul = new double***[nx+2];
    Ur = new double***[nx+2];
    Wl = new double***[nx+2];
    Wr = new double***[nx+2];
    for (int i=0; i<(nx+2)*K; i++){
      dt[i] = new double*[(ny+2)*K];
      dUdt[i] = new double***[(ny+2)*K];
      phi[i] = new double**[(ny+2)*K];
      dWdx[i] = new double**[(ny+2)*K];
      dWdy[i] = new double**[(ny+2)*K];
      dWdz[i] = new double**[(ny+2)*K];
      W[i] = new double**[(ny+2)*K];
      U[i] = new double**[(ny+2)*K];
      Wn[i] = new double**[(ny+2)*K];
      Un[i] = new double**[(ny+2)*K];
      W0[i] = new double**[(ny+2)*K];
      U0[i] = new double**[(ny+2)*K];
      F[i] = new double**[(ny+2)*K];
      for (int j=0; j<(ny+2)*K; j++){
        dt[i][j] = new double[(nz+2)*K];
        dUdt[i][j] = new double**[(nz+2)*K];
        phi[i][j] = new double*[(nz+2)*K];
        dWdx[i][j] = new double*[(nz+2)*K];
        dWdy[i][j] = new double*[(nz+2)*K];
        dWdz[i][j] = new double*[(nz+2)*K];
        W[i][j] = new double*[(nz+2)*K];
        U[i][j] = new double*[(nz+2)*K];
        Wn[i][j] = new double*[(nz+2)*K];
        Un[i][j] = new double*[(nz+2)*K];
        W0[i][j] = new double*[(nz+2)*K];
        U0[i][j] = new double*[(nz+2)*K];
        F[i][j] = new double*[(nz+2)*K];
        for (int k=0; k<(nz+2)*K; k++){
          dUdt[i][j][k] = new double*[nResidual];
          phi[i][j][k] = new double[5];
          dWdx[i][j][k] = new double[5];
          dWdy[i][j][k] = new double[5];
          dWdz[i][j][k] = new double[5];
          W[i][j][k] = new double[5];
          U[i][j][k] = new double[5];
          Wn[i][j][k] = new double[5];
          Un[i][j][k] = new double[5];
          W0[i][j][k] = new double[5];
          U0[i][j][k] = new double[5];
          F[i][j][k] = new double[5];
          for (int n=0; n<nResidual; n++){
            dUdt[i][j][k][n] = new double[5];
          }
        }
      }
    }
    for (int i=0; i<nx+1; i++){
      Fupw[i] = new double**[ny+1];
      for (int j=0; j<ny+1; j++){
        Fupw[i][j] = new double*[nz+1];
        for (int k=0; k<nz+1; k++){
          Fupw[i][j][k] = new double[5];
        }
      }
    }
    for (int i=0; i<nx+2; i++){
      Fl[i] = new double**[ny+2];
      Fr[i] = new double**[ny+2];
      Ul[i] = new double**[ny+2];
      Ur[i] = new double**[ny+2];
      Wl[i] = new double**[ny+2];
      Wr[i] = new double**[ny+2];
      for (int j=0; j<ny+2; j++){
        Fl[i][j] = new double*[nz+2];
        Fr[i][j] = new double*[nz+2];
        Ul[i][j] = new double*[nz+2];
        Ur[i][j] = new double*[nz+2];
        Wl[i][j] = new double*[nz+2];
        Wr[i][j] = new double*[nz+2];
        for (int k=0; k<nz+2; k++){
          Fl[i][j][k] = new double[5];
          Fr[i][j][k] = new double[5];
          Ul[i][j][k] = new double[5];
          Ur[i][j][k] = new double[5];
          Wl[i][j][k] = new double[5];
          Wr[i][j][k] = new double[5];
        }
      }
    }
    switch (LimiterType){
      case 0:
        strcpy(LimiterName, "FirstOrder");
        break;
      case 1:
        strcpy(LimiterName, "NoLimiter");
        break;
      case 2:
        strcpy(LimiterName, "Barth-Jespersen");
        break;
      case 3:
        strcpy(LimiterName, "Venkatakrishnan");
        break;
    }
  }
  // Deallocate
  void Deallocate(){
    mesh.Deallocate();
    FR.Deallocate();
    for (int i=0; i<(nx+2)*K; i++){
      for (int j=0; j<(ny+2)*K; j++){
        for (int k=0; k<(nz+2)*K; k++){
          for (int n=0; n<nResidual; n++){
            if (dUdt!=NULL) delete[] dUdt[i][j][k][n];
          }
          if (dUdt!=NULL) delete[] dUdt[i][j][k];
          if (phi!=NULL) delete[] phi[i][j][k];
          if (dWdx!=NULL) delete[] dWdx[i][j][k];
          if (dWdy!=NULL) delete[] dWdy[i][j][k];
          if (dWdz!=NULL) delete[] dWdz[i][j][k];
          if (W!=NULL) delete[] W[i][j][k];
          if (U!=NULL) delete[] U[i][j][k];
          if (Wn!=NULL) delete[] Wn[i][j][k];
          if (Un!=NULL) delete[] Un[i][j][k];
          if (W0!=NULL) delete[] W0[i][j][k];
          if (U0!=NULL) delete[] U0[i][j][k];
          if (F!=NULL) delete[] F[i][j][k];
        }
        if (dt!=NULL) delete[] dt[i][j];
        if (dUdt!=NULL) delete[] dUdt[i][j];
        if (phi!=NULL) delete[] phi[i][j];
        if (dWdx!=NULL) delete[] dWdx[i][j];
        if (dWdy!=NULL) delete[] dWdy[i][j];
        if (dWdz!=NULL) delete[] dWdz[i][j];
        if (W!=NULL) delete[] W[i][j];
        if (U!=NULL) delete[] U[i][j];
        if (Wn!=NULL) delete[] Wn[i][j];
        if (Un!=NULL) delete[] Un[i][j];
        if (W0!=NULL) delete[] W0[i][j];
        if (U0!=NULL) delete[] U0[i][j];
        if (F!=NULL) delete[] F[i][j];
      }
      if (dt!=NULL) delete[] dt[i];
      if (dUdt!=NULL) delete[] dUdt[i];
      if (phi!=NULL) delete[] phi[i];
      if (dWdx!=NULL) delete[] dWdx[i];
      if (dWdy!=NULL) delete[] dWdy[i];
      if (dWdz!=NULL) delete[] dWdz[i];
      if (W!=NULL) delete[] W[i];
      if (U!=NULL) delete[] U[i];
      if (Wn!=NULL) delete[] Wn[i];
      if (Un!=NULL) delete[] Un[i];
      if (W0!=NULL) delete[] W0[i];
      if (U0!=NULL) delete[] U0[i];
      if (F!=NULL) delete[] F[i];
    }
    for (int i=0; i<nx+1; i++){
      for (int j=0; j<ny+1; j++){
        for (int k=0; k<nz+1; k++){
          if (Fupw!=NULL) delete[] Fupw[i][j][k];
        }
        if (Fupw!=NULL) delete[] Fupw[i][j];
      }
      if (Fupw!=NULL) delete[] Fupw[i];
    }
    for (int i=0; i<nx+2; i++){
      for (int j=0; j<ny+2; j++){
        for (int k=0; k<nz+2; k++){
          if (Fl!=NULL) delete[] Fl[i][j][k];
          if (Fr!=NULL) delete[] Fr[i][j][k];
          if (Ul!=NULL) delete[] Ul[i][j][k];
	  if (Ur!=NULL) delete[] Ur[i][j][k];
	  if (Wl!=NULL) delete[] Wl[i][j][k];
          if (Wr!=NULL) delete[] Wr[i][j][k];
        }
        if (Fl!=NULL) delete[] Fl[i][j];
        if (Fr!=NULL) delete[] Fr[i][j];
        if (Ul!=NULL) delete[] Ul[i][j];
        if (Ur!=NULL) delete[] Ur[i][j];
        if (Wl!=NULL) delete[] Wl[i][j];
        if (Wr!=NULL) delete[] Wr[i][j];
      }
      if (Fl!=NULL) delete[] Fl[i];
      if (Fr!=NULL) delete[] Fr[i];
      if (Ul!=NULL) delete[] Ul[i];
      if (Ur!=NULL) delete[] Ur[i];
      if (Wl!=NULL) delete[] Wl[i];
      if (Wr!=NULL) delete[] Wr[i];
    }
    if (dt!=NULL) delete[] dt;
    if (dUdt!=NULL) delete[] dUdt;
    if (phi!=NULL) delete[] phi;
    if (dWdx!=NULL) delete[] dWdx;
    if (dWdy!=NULL) delete[] dWdy;
    if (dWdz!=NULL) delete[] dWdz;
    if (W!=NULL) delete[] W;
    if (U!=NULL) delete[] U;
    if (Wn!=NULL) delete[] Wn;
    if (Un!=NULL) delete[] Un;
    if (W0!=NULL) delete[] W0;
    if (U0!=NULL) delete[] U0;
    if (F!=NULL) delete[] F;
    if (Fupw!=NULL) delete[] Fupw;
    if (Fl!=NULL) delete[] Fl;
    if (Fr!=NULL) delete[] Fr;
    if (Ul!=NULL) delete[] Ul;
    if (Ur!=NULL) delete[] Ur;
    if (Wl!=NULL) delete[] Wl;
    if (Wr!=NULL) delete[] Wr;
  }
  // Nullify
  void Nullify(){
    this->K = -1;
    this->nx = -1;
    this->ny = -1;
    this->nz = -1;
    this->nvar = -1;
    this->SolverType = -1;
    this->EquationType = -1;
    this->LimiterType = -1;
    this->BCtype = -1;
    this->Ktype = -1;
    this->FRtype = -1;
    mesh.Nullify();
    FR.Nullify();
    dt = NULL;
    dUdt = NULL;
    phi = NULL;
    dWdx = NULL;
    dWdy = NULL;
    dWdz = NULL;
    W = NULL;
    U = NULL;
    Wn = NULL;
    Un = NULL;
    W0 = NULL;
    U0 = NULL;
    F = NULL;
    Fupw = NULL;
    Fl = NULL;
    Fr = NULL;
    Ul = NULL;
    Ur = NULL;
    Wl = NULL;
    Wr = NULL;
  }
  void setMesh(const int& option, 
               const double& x0, const double& x1, 
               const double& y0, const double& y1, 
               const double& z0, const double& z1);
  double TimeStep(const double& CFL);
  void ApplyBCs(void);
  void UpdatePrimitive(void);
  void EnforcePositivity(void);
  void Gradient(void);
  void Limiter(const int& i, const int& j, const int& k);
  void Limiter(void);
  void ResetResidual(const int& index);
  void Residual(const int& index);
  void FVMResidual(const int& index);
  void FRResidual(const int& index);
  void CommonFluxXdir(double* commonflux, int i, int j, int k, int jp, int kp);
  void CommonFluxYdir(double* commonflux, int i, int j, int k, int ip, int kp);
  void CommonFluxZdir(double* commonflux, int i, int j, int k, int ip, int jp);
  void LinearAdvectionInitialCondition(const int& option, const double& ncycle, double& t);
  void LinearAdvectionExactSolution(const int& option, const double& t);
  void LinearAdvectionOutput(const int& option, const double& CFL, const int& nStage, const double& t);
  void EulerInitialCondition(const int& option, double& t);
  void EulerExactSolution(const int& option);
  void EulerOutput(const int& option, const double& CFL, const int& nStage);
  void RungeKutta(const int& nstage);
};


///////////////////////////////////////////////////////////////////////////////
// MESH
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// setMesh (1D)
// **************************************************************************//
inline void Solver1D::setMesh(const int& option, 
                              const double& x0, const double& x1){

  // DESCRIPTION
  // ----------------------
  // Generate 1D structured mesh.


  // INPUTS
  // ----------------------
  // option  - Mesh option
  // x0      - Domain boundary
  // x1      - Domain boundary


  // SET 1D MESH
  // ----------------------
  switch (option){
    case 1: // Uniform mesh
      mesh.setUniformMesh(x0, x1);
      break;
    default:
      throw invalid_argument("Only uniform mesh avalaible"); 
  }
}


// **************************************************************************//
// setMesh (2D)
// **************************************************************************//
inline void Solver2D::setMesh(const int& option, 
                              const double& x0, const double& x1,
                              const double& y0, const double& y1){

  // DESCRIPTION
  // ----------------------
  // Generate 2D structured mesh.


  // INPUTS
  // ----------------------
  // option  - Mesh option
  // x0      - Domain boundary
  // x1      - Domain boundary
  // y0      - Domain boundary
  // y1      - Domain boundary


  // SET 2D MESH
  // ----------------------
  switch (option){
    case 1: // Uniform mesh
      mesh.setUniformMesh(x0, x1, y0, y1);
      break;
    default:
      throw invalid_argument("Only uniform mesh avalaible"); 
  }
}


// **************************************************************************//
// setMesh (3D)
// **************************************************************************//
inline void Solver3D::setMesh(const int& option, 
                              const double& x0, const double& x1,
                              const double& y0, const double& y1,
                              const double& z0, const double& z1){

  // DESCRIPTION
  // ----------------------
  // Generate 3D structured mesh.

    
  // INPUTS
  // ----------------------
  // option  - Mesh option
  // x0      - Domain boundary
  // x1      - Domain boundary
  // y0      - Domain boundary
  // y1      - Domain boundary
  // z0      - Domain boundary
  // z1      - Domain boundary


  // SET 3D MESH
  // ----------------------
  switch (option){
    case 1: // Uniform mesh
      mesh.setUniformMesh(x0, x1, y0, y1, z0, z1);
      break;
    default:
      throw invalid_argument("Only uniform mesh avalaible"); 
  }
}


///////////////////////////////////////////////////////////////////////////////
// TIMESTEP
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// TimeStep (1D)
// **************************************************************************//
inline double Solver1D::TimeStep(const double& CFL){
  
  // DESCRIPTION
  // ----------------------
  // Compute maximum timestep from the CFL number.


  // INPUTS
  // ----------------------
  // CFL     - CFL number


  // OUTPUTS
  // ----------------------
  // dtmin   - (Global) minimum timestep


  // VARIABLE DECLARATION
  // ----------------------
  double dtmin=1e9;
  double dx, rho, u, p, g=Euler1D.getHeatRatio();


  // TIMESTEP
  // ----------------------
  // Initialize timestep
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      dt[K*i+ip] = dtmin;
    }
  }

  // Local timestep
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      dx = mesh.dx[i];
      switch (EquationType){  
        case 1: // Linear Advection
          dt[K*i+ip] = CFL*dx/fabs(wavespeed);
          break;
        case 2: // Euler
          Euler1D.ConstoPrim(&U[K*i+ip][0], rho, u, p);
          dt[K*i+ip] = CFL*dx/(fabs(u) + sqrt(g*p/rho));
          break;
        default:
          throw invalid_argument("Unknown solver type"); 
      }
      dtmin = fmin(dtmin, dt[K*i+ip]);
    }
  }

  // Global minimum
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      dt[K*i+ip] = dtmin;
    }
  }


  // RETURN VARIABLE
  // ----------------------
  return dtmin;
}


// **************************************************************************//
// TimeStep (2D)
// **************************************************************************//
inline double Solver2D::TimeStep(const double& CFL){
  
  // DESCRIPTION
  // ----------------------
  // Compute maximum timestep from the CFL number.


  // INPUTS
  // ----------------------
  // CFL     - CFL number


  // OUTPUTS
  // ----------------------
  // dtmin   - (Global) minimum timestep


  // VARIABLE DECLARATION
  // ----------------------
  double dtmin=1e9;
  double dx, dy, rho, u, v, p, g=Euler2D.getHeatRatio();


  // TIMESTEP
  // ----------------------
  // Initialize timestep
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int ip=0; ip<K; ip++){
        for (int jp=0; jp<K; jp++){
          dt[K*i+ip][K*j+jp] = dtmin;
        }
      }
    }
  }

  // Local timestep
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int ip=0; ip<K; ip++){
        for (int jp=0; jp<K; jp++){
          dx = mesh.dx[i][j];
          dy = mesh.dy[i][j];
          switch (EquationType){  
            case 1: // Linear Advection
              dt[K*i+ip][K*j+jp] = CFL*(dx*dy)/(dx+dy)/fabs(wavespeed);
              break;
            case 2: // Euler
              Euler2D.ConstoPrim(&U[K*i+ip][K*j+jp][0], rho, u, v, p);
              //dt[K*i+ip][K*j+jp] = CFL*fmin(dx/(fabs(u) + sqrt(g*p/rho)),
              //                              dy/(fabs(v) + sqrt(g*p/rho)));
              dt[K*i+ip][K*j+jp] = CFL*(dx*dy)/(dy*(fabs(u) + sqrt(g*p/rho)) + 
                                                dx*(fabs(v) + sqrt(g*p/rho)));
             break;
             default:
              throw invalid_argument("Unknown solver type"); 
          }
          dtmin = fmin(dtmin, dt[K*i+ip][K*j+jp]);
        }
      }
    }
  }

  // Global minimum
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int ip=0; ip<K; ip++){
        for (int jp=0; jp<K; jp++){
          dt[K*i+ip][K*j+jp] = dtmin;
        }
      }
    }
  }


  // RETURN VARIABLE
  // ----------------------
  return dtmin;
}


// **************************************************************************//
// TimeStep (3D)
// **************************************************************************//
inline double Solver3D::TimeStep(const double& CFL){
  
  // DESCRIPTION
  // ----------------------
  // Compute maximum timestep from the CFL number.


  // INPUTS
  // ----------------------
  // CFL     - CFL number


  // OUTPUTS
  // ----------------------
  // dtmin   - (Global) minimum timestep


  // VARIABLE DECLARATION
  // ----------------------
  double dtmin=1e9;
  double dx, dy, dz, rho, u, v, w, p, g=Euler3D.getHeatRatio();


  // TIMESTEP
  // ----------------------
  // Initialize timestep
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
        for (int ip=0; ip<K; ip++){
          for (int jp=0; jp<K; jp++){
            for (int kp=0; kp<K; kp++){
              dt[K*i+ip][K*j+jp][K*k+kp] = dtmin;
            }
          }
        }
      }
    }
  }

  // Local timestep
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
        for (int ip=0; ip<K; ip++){
          for (int jp=0; jp<K; jp++){
            for (int kp=0; kp<K; kp++){
              dx = mesh.dx[i][j][k];
              dy = mesh.dy[i][j][k];
              dz = mesh.dz[i][j][k];
              switch (EquationType){  
                case 1: // Linear Advection
                  dt[K*i+ip][K*j+jp][K*k+kp] = CFL*(dx*dy*dz)/(dx*dy+dy*dz+dx*dz)/fabs(wavespeed);
                  break;
                case 2: // Euler
                  Euler3D.ConstoPrim(&U[K*i+ip][K*j+jp][K*k+kp][0], rho, u, v, w, p);
                  //dt[K*i+ip][K*j+jp][K*k+kp] = CFL*fmin(fmin(dx/(fabs(u) + sqrt(g*p/rho)),
                  //                                           dy/(fabs(v) + sqrt(g*p/rho))),
                  //                                           dz/(fabs(w) + sqrt(g*p/rho)));
                  dt[K*i+ip][K*j+jp][K*k+kp] = CFL*(dx*dy*dx)/(dy*dz*(fabs(u) + sqrt(g*p/rho)) + 
                                                               dx*dz*(fabs(v) + sqrt(g*p/rho)) + 
                                                               dx*dy*(fabs(w) + sqrt(g*p/rho)));
                  break;
                default:
                  throw invalid_argument("Unknown solver type"); 
              }
              dtmin = fmin(dtmin, dt[K*i+ip][K*j+jp][K*k+kp]);
            }
          }
        }
      }
    }
  }

  // Global minimum
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
        for (int ip=0; ip<K; ip++){
          for (int jp=0; jp<K; jp++){
            for (int kp=0; kp<K; kp++){
              dt[K*i+ip][K*j+jp][K*k+kp] = dtmin;
            }
          }
        }
      }
    }
  }


  // RETURN VARIABLE
  // ----------------------
  return dtmin;
}


///////////////////////////////////////////////////////////////////////////////
// BOUNDARY CONDITIONS
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// ApplyBCs (1D)
// **************************************************************************//
inline void Solver1D::ApplyBCs(void){

  // DESCRIPTION
  // ----------------------
  // Apply boundary conditions (BCs).
 
 
  // APPLY BCs
  // ----------------------
  switch (BCtype){
    case 1: // Fixed
      for (int i=0; i<mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          for (int n=0; n<nvar; n++){
            U[K*i+ip][n] = U0[K*i+ip][n]; 
            U[K*(mesh.nx+mesh.ng+i)+ip][n] = U0[K*(mesh.nx+mesh.ng+i)+ip][n];
          }
        }
      }
      break;
    case 2: // Periodic
      for (int i=0; i<mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          for (int n=0; n<nvar; n++){
            U[K*i+ip][n] = U[K*(mesh.nx+i)+ip][n]; 
            U[K*(mesh.nx+mesh.ng+i)+ip][n] = U[K*(mesh.ng+i)+ip][n];
          }
        }
      }
      break;
    case 3: // None
      for (int i=0; i<mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          for (int n=0; n<nvar; n++){
            U[K*i+ip][n] = U[K*mesh.ng][n]; 
            U[K*(mesh.nx+mesh.ng+i)+ip][n] = U[K*(mesh.nx+mesh.ng)-1][n];
          }
        }
      }
      break;
    //default:
    //  throw invalid_argument("Unknown boundary condition type"); 
  }
}


// **************************************************************************//
// ApplyBCs (2D)
// **************************************************************************//
inline void Solver2D::ApplyBCs(void){

  // DESCRIPTION
  // ----------------------
  // Apply boundary conditions (BCs).
  

  // APPLY BCs
  // ----------------------
  switch (BCtype){
    case 1: // Fixed
      for (int j=0; j<mesh.ny+2*mesh.ng; j++){
        for (int jp=0; jp<K; jp++){
          for (int i=0; i<mesh.ng; i++){
            for (int ip=0; ip<K; ip++){
              for (int n=0; n<nvar; n++){
                U[K*i+ip][K*j+jp][n] = U0[K*i+ip][K*j+jp][n]; 
                U[K*(mesh.nx+mesh.ng+i)+ip][K*j+jp][n] = U0[K*(mesh.nx+mesh.ng+i)][K*j+jp][n];
              }
            }
          }
        }
      }
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int j=0; j<mesh.ng; j++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){
              for (int n=0; n<nvar; n++){
                U[K*i+ip][K*j+jp][n] = U0[K*i+ip][K*j+jp][n]; 
                U[K*i+ip][K*(mesh.ny+mesh.ng+j)+jp][n] = U0[K*i+ip][K*(mesh.ny+mesh.ng+j)+jp][n];
              }
            }
          }
        }
      }
      break;
    case 2: // Periodic
      for (int j=0; j<mesh.ny+2*mesh.ng; j++){
        for (int i=0; i<mesh.ng; i++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){
              for (int n=0; n<nvar; n++){
                U[K*i+ip][K*j+jp][n] = U[K*(mesh.nx+i)+ip][K*j+jp][n]; 
                U[K*(mesh.nx+mesh.ng+i)+ip][K*j+jp][n] = U[K*(mesh.ng+i)+ip][K*j+jp][n];
              }
            }
          }
        }
      }
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int j=0; j<mesh.ng; j++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){
              for (int n=0; n<nvar; n++){
                U[K*i+ip][K*j+jp][n] = U[K*i+ip][K*(mesh.ny+j)+jp][n]; 
                U[K*i+ip][K*(mesh.ny+mesh.ng+j)+jp][n] = U[K*i+ip][K*(mesh.ng+j)+jp][n];
              }
            }
          }
        }
      }
      break;
    case 3: // None
      for (int j=0; j<mesh.ny+2*mesh.ng; j++){
        for (int i=0; i<mesh.ng; i++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){
              for (int n=0; n<nvar; n++){
                U[K*i+ip][K*j+jp][n] = U[K*mesh.ng][K*j+jp][n]; 
                U[K*(mesh.nx+mesh.ng+i)+ip][K*j+jp][n] = U[K*(mesh.nx+mesh.ng)-1][K*j+jp][n];
              }
            }
          }
        }
      }
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int j=0; j<mesh.ng; j++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){
              for (int n=0; n<nvar; n++){
                U[K*i+ip][K*j+jp][n] = U[K*i+ip][K*mesh.ng][n]; 
                U[K*i+ip][K*(mesh.ny+mesh.ng+j)+jp][n] = U[K*i+ip][K*(mesh.ny+mesh.ng)-1][n];
              }
            }
          }
        }
      }
      break;
    case 4: // Wedge
      for (int j=0; j<mesh.ny+2*mesh.ng; j++){
        for (int i=0; i<mesh.ng; i++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){
              for (int n=0; n<nvar; n++){
                U[K*i+ip][K*j+jp][n] = U0[K*i+ip][K*j+jp][n]; 
                U[K*(mesh.nx+mesh.ng+i)+ip][K*j+jp][n] = U[K*(mesh.nx+mesh.ng)-1][K*j+jp][n];
              }
            }
          }
        }
      }
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int j=0; j<mesh.ng; j++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){
              for (int n=0; n<nvar; n++){
                U[K*i+ip][K*j+jp][n] = U[K*i+ip][K*mesh.ng][n]; 
                U[K*i+ip][K*(mesh.ny+mesh.ng+j)+jp][n] = U0[K*i+ip][K*(mesh.ny+mesh.ng+j)+jp][n];
              }
              U[K*i+ip][K*j+jp][2] = 0.0;
            }
          }
        }
      }
      break;
    //default:
    //  throw invalid_argument("Unknown boundary condition type"); 
  }
}


// **************************************************************************//
// ApplyBCs (3D)
// **************************************************************************//
inline void Solver3D::ApplyBCs(void){

  // DESCRIPTION
  // ----------------------
  // Apply boundary conditions (BCs).
  

  // APPLY BCs
  // ----------------------
  switch (BCtype){
    case 1: // Fixed
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int j=0; j<mesh.ny+2*mesh.ng; j++){
          for (int i=0; i<mesh.ng; i++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  for (int n=0; n<nvar; n++){
                    U[K*i+ip][K*j+jp][K*k+kp][n] = U0[K*i+ip][K*j+jp][K*k+kp][n]; 
                    U[K*(mesh.nx+mesh.ng+i)+ip][K*j+jp][K*k+kp][n] = U0[K*(mesh.nx+mesh.ng+i)+ip][K*j+jp][K*k+kp][n];
                  }
                }
              } 
            }
          }
        }
      }
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int i=0; i<mesh.nx+2*mesh.ng; i++){
          for (int j=0; j<mesh.ng; j++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  for (int n=0; n<nvar; n++){
                    U[K*i+ip][K*j+jp][K*k+kp][n] = U0[K*i+ip][K*j+jp][K*k+kp][n]; 
                    U[K*i+ip][K*(mesh.ny+mesh.ng+j)+jp][K*k+kp][n] = U0[K*i+ip][K*(mesh.ny+mesh.ng+j)+jp][K*k+kp][n];
                  }
                }
              }
            }
          }
        }
      }
      for (int j=0; j<mesh.ny+2*mesh.ng; j++){
        for (int i=0; i<mesh.nx+2*mesh.ng; i++){
          for (int k=0; k<mesh.ng; k++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  for (int n=0; n<nvar; n++){
                    U[K*i+ip][K*j+jp][K*k+kp][n] = U0[K*i+ip][K*j+jp][K*k+kp][n]; 
                    U[K*i+ip][K*j+jp][K*(mesh.nz+mesh.ng+k)+kp][n] = U0[K*i+ip][K*j+jp][K*(mesh.nz+mesh.ng+k)+kp][n];
                  }
                }
              }
            }
          }
        }
      }
      break;
    case 2: // Periodic
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int j=0; j<mesh.ny+2*mesh.ng; j++){
          for (int i=0; i<mesh.ng; i++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  for (int n=0; n<nvar; n++){
                    U[K*i+ip][K*j+jp][K*k+kp][n] = U[K*(mesh.nx+i)+ip][K*j+jp][K*k+kp][n]; 
                    U[K*(mesh.nx+mesh.ng+i)+ip][K*j+jp][K*k+kp][n] = U[K*(mesh.ng+i)+ip][K*j+jp][K*k+kp][n];
                  }
                }
              }
            }
          }
        }
      }
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int i=0; i<mesh.nx+2*mesh.ng; i++){
          for (int j=0; j<mesh.ng; j++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  for (int n=0; n<nvar; n++){
                    U[K*i+ip][K*j+jp][K*k+kp][n] = U[K*i+ip][K*(mesh.ny+j)+jp][K*k+kp][n]; 
                    U[K*i+ip][K*(mesh.ny+mesh.ng+j)+jp][K*k+kp][n] = U[K*i+ip][K*(mesh.ng+j)+jp][K*k+kp][n];
                  }
                }
              }
            }
          }
        }
      }
      for (int j=0; j<mesh.ny+2*mesh.ng; j++){
        for (int i=0; i<mesh.nx+2*mesh.ng; i++){
          for (int k=0; k<mesh.ng; k++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  for (int n=0; n<nvar; n++){
                    U[K*i+ip][K*j+jp][K*k+kp][n] = U[K*i+ip][K*j+jp][K*(mesh.nz+k)+kp][n]; 
                    U[K*i+ip][K*j+jp][K*(mesh.nz+mesh.ng+k)+kp][n] = U[K*i+ip][K*j+jp][K*(mesh.ng+k)+kp][n];
                  }
                }
              }
            }
          }
        }
      }
      break;
    case 3: // None
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int j=0; j<mesh.ny+2*mesh.ng; j++){
          for (int i=0; i<mesh.ng; i++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  for (int n=0; n<nvar; n++){
                    U[K*i+ip][K*j+jp][K*k+kp][n] = U[K*(mesh.ng+i)+ip][K*j+jp][K*k+kp][n]; 
                    U[K*(mesh.nx+mesh.ng+i)+ip][K*j+jp][K*k+kp][n] = U[K*(mesh.nx+mesh.ng)-1][K*j+jp][K*k+kp][n];
                  }
                }
              }
            }
          }
        }
      }
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int i=0; i<mesh.nx+2*mesh.ng; i++){
          for (int j=0; j<mesh.ng; j++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  for (int n=0; n<nvar; n++){
                    U[K*i+ip][K*j+jp][K*k+kp][n] = U[K*i+ip][K*(mesh.ng+j)+jp][K*k+kp][n]; 
                    U[K*i+ip][K*(mesh.ny+mesh.ng+j)+jp][K*k+kp][n] = U[K*i+ip][K*(mesh.ny+mesh.ng)-1][K*k+kp][n];
                  }
                }
              }
            }
          }
        }
      }
      for (int j=0; j<mesh.ny+2*mesh.ng; j++){
        for (int i=0; i<mesh.nx+2*mesh.ng; i++){
          for (int k=0; k<mesh.ng; k++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  for (int n=0; n<nvar; n++){
                    U[K*i+ip][K*j+jp][K*k+kp][n] = U[K*i+ip][K*j+jp][K*(mesh.ng+k)+kp][n]; 
                    U[K*i+ip][K*j+jp][K*(mesh.nz+mesh.ng+k)+kp][n] = U[K*i+ip][K*j+jp][K*(mesh.nz+mesh.ng)-1][n];
                  }
                }
              }
            }
          }
        }
      }
      break;
    case 4: // Wedge
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int j=0; j<mesh.ny+2*mesh.ng; j++){
          for (int i=0; i<mesh.ng; i++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  for (int n=0; n<nvar; n++){
                    U[K*i+ip][K*j+jp][K*k+kp][n] = U0[K*i+ip][K*j+jp][K*k+kp][n]; 
                    U[K*(mesh.nx+mesh.ng+i)+ip][K*j+jp][K*k+kp][n] = U[K*(mesh.nx+mesh.ng)-1][K*j+jp][K*k+kp][n];
                  }
                }
              }
            }
          }
        }
      }
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int i=0; i<mesh.nx+2*mesh.ng; i++){
          for (int j=0; j<mesh.ng; j++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  for (int n=0; n<nvar; n++){
                    U[K*i+ip][K*j+jp][K*k+kp][n] = U[K*i+ip][K*mesh.ng][K*k+kp][n]; 
                    U[K*i+ip][K*(mesh.ny+mesh.ng+j)+jp][K*k+kp][n] = U[K*i+ip][K*(mesh.ny+mesh.ng)-1][K*k+kp][n];
                  }
                }
              }
            }
          }
        }
      }
      for (int j=0; j<mesh.ny+2*mesh.ng; j++){
        for (int i=0; i<mesh.nx+2*mesh.ng; i++){
          for (int k=0; k<mesh.ng; k++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  for (int n=0; n<nvar; n++){
                    U[K*i+ip][K*j+jp][K*k+kp][n] = U[K*i+ip][K*j+jp][K*mesh.ng][n]; 
                    U[K*i+ip][K*j+jp][K*(mesh.nz+mesh.ng+k)+kp][n] = U0[K*i+ip][K*j+jp][K*(mesh.nz+mesh.ng+k)+kp][n];
                  }
                  U[K*i+ip][K*j+jp][K*k+kp][3] = 0.0;
                }
              }
            }
          }
        }
      }
      break;
    //default:
    //  throw invalid_argument("Unknown boundary condition type"); 
  }
}


///////////////////////////////////////////////////////////////////////////////
// UPDATE PRIMITIVE
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// UpdatePrimitive (1D)
// **************************************************************************//
inline void Solver1D::UpdatePrimitive(void){

  // DESCRIPTION
  // ----------------------
  // Update primitive variable, W, from conservative variable, U.


  // UPDATE PRIMITIVE VAR
  // ----------------------
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      if (EquationType!=2){
        W[K*i+ip][0]=U[K*i+ip][0];
      }else{
        Euler1D.ConstoPrim(&U[K*i+ip][0],&W[K*i+ip][0]); 
      }
    }
  }
}


// **************************************************************************//
// UpdatePrimitive (2D)
// **************************************************************************//
inline void Solver2D::UpdatePrimitive(void){

  // DESCRIPTION
  // ----------------------
  // Update primitive variable, W, from conservative variable, U.


  // UPDATE PRIMITIVE VAR
  // ----------------------
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int ip=0; ip<K; ip++){
        for (int jp=0; jp<K; jp++){
          if (EquationType!=2){
            W[K*i+ip][K*j+jp][0]=U[K*i+ip][K*j+jp][0];
          }else{
            Euler2D.ConstoPrim(&U[K*i+ip][K*j+jp][0],&W[K*i+ip][K*j+jp][0]); 
          }
        }
      }
    }
  }
}


// **************************************************************************//
// UpdatePrimitive (3D)
// **************************************************************************//
inline void Solver3D::UpdatePrimitive(void){

  // DESCRIPTION
  // ----------------------
  // Update primitive variable, W, from conservative variable, U.

  // UPDATE PRIMITIVE VAR
  // ----------------------
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int ip=0; ip<K; ip++){
          for (int jp=0; jp<K; jp++){
            for (int kp=0; kp<K; kp++){
              if (EquationType!=2){
                W[K*i+ip][K*j+jp][K*k+kp][0]=U[K*i+ip][K*j+jp][K*k+kp][0];
              }else{
                Euler3D.ConstoPrim(&U[K*i+ip][K*j+jp][K*k+kp][0],&W[K*i+ip][K*j+jp][K*k+kp][0]);
              }
            } 
          }
        }
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// ENFORCE POSITIVITY
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// EnforcePositivity (1D)
// **************************************************************************//
inline void Solver1D::EnforcePositivity(void){

  // DESCRIPTION
  // ----------------------
  // Enforce positivity of density and pressure for the Euler equations.


  // ENFORCE POSITIVITY
  // ----------------------
  if (EquationType!=2) return;
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      // Enforce positivity of density and pressure
      W[K*i+ip][0] = fmax(1.0E-08, W[K*i+ip][0]);  
      W[K*i+ip][2] = fmax(1.0E-08, W[K*i+ip][2]);  
      // Update conservative variable
      Euler1D.PrimtoCons(&W[K*i+ip][0],&U[K*i+ip][0]); 
    }
  }
}


// **************************************************************************//
// EnforcePositivity (2D)
// **************************************************************************//
inline void Solver2D::EnforcePositivity(void){

  // DESCRIPTION
  // ----------------------
  // Enforce positivity of density and pressure for the Euler equations.


  // ENFORCE POSITIVITY
  // ----------------------
  if (EquationType!=2) return;
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int ip=0; ip<K; ip++){
        for (int jp=0; jp<K; jp++){
          // Enforce positivity of density and pressure
          W[K*i+ip][K*j+jp][0] = fmax(1.0E-08, W[K*i+ip][K*j+jp][0]);
          W[K*i+ip][K*j+jp][3] = fmax(1.0E-08, W[K*i+ip][K*j+jp][3]);
          // Update conservative variable
          Euler2D.PrimtoCons(&W[K*i+ip][K*j+jp][0],&U[K*i+ip][K*j+jp][0]); 
        }
      }
    }
  }
}


// **************************************************************************//
// EnforcePositivity (3D)
// **************************************************************************//
inline void Solver3D::EnforcePositivity(void){

  // DESCRIPTION
  // ----------------------
  // Enforce positivity of density and pressure for the Euler equations.


  // ENFORCE POSITIVITY
  // ----------------------
  if (EquationType!=2) return;
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
        for (int ip=0; ip<K; ip++){
          for (int jp=0; jp<K; jp++){
            for (int kp=0; kp<K; kp++){
              // Enforce positivity of density and pressure
              W[K*i+ip][K*j+jp][K*k+kp][0] = fmax(1.0E-08, W[K*i+ip][K*j+jp][K*k+kp][0]);  
              W[K*i+ip][K*j+jp][K*k+kp][4] = fmax(1.0E-08, W[K*i+ip][K*j+jp][K*k+kp][4]);  
              // Update conservative variable
              Euler3D.PrimtoCons(&W[K*i+ip][K*j+jp][K*k+kp][0],&U[K*i+ip][K*j+jp][K*k+kp][0]); 
            }
          }
        }
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// GRADIENT
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// Gradient (1D)
// **************************************************************************//
inline void Solver1D::Gradient(void){

  // DESCRIPTION
  // ----------------------
  // Compute the 1D solution gradient.
 

  // 1D GRADIENT
  // ----------------------
  if(SolverType == 2) return;
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for(int n=0; n<nvar; n++){
      dWdx[i][n] = 0.25*(W[i+1][n] - W[i-1][n]);
    }
  }
}


// **************************************************************************//
// Gradient (2D)
// **************************************************************************//
inline void Solver2D::Gradient(void){

  // DESCRIPTION
  // ----------------------
  // Compute the 2D solution gradient.
 

  // 2D GRADIENT
  // ----------------------
  if(SolverType == 2) return;
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for(int n=0; n<nvar; n++){
        dWdx[i][j][n] = 0.25*(W[i+1][j][n] - W[i-1][j][n]);
        dWdy[i][j][n] = 0.25*(W[i][j+1][n] - W[i][j-1][n]);
      }
    }
  }
}


// **************************************************************************//
// Gradient (3D)
// **************************************************************************//
inline void Solver3D::Gradient(void){

  // DESCRIPTION
  // ----------------------
  // Compute the 3D solution gradient.
 

  // 3D GRADIENT
  // ----------------------
  if(SolverType == 2) return;
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
        for(int n=0; n<nvar; n++){
          dWdx[i][j][k][n] = 0.25*(W[i+1][j][k][n] - W[i-1][j][k][n]);
          dWdy[i][j][k][n] = 0.25*(W[i][j+1][k][n] - W[i][j-1][k][n]);
          dWdz[i][j][k][n] = 0.25*(W[i][j][k+1][n] - W[i][j][k-1][n]);
        }
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// LIMITER
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// Limiter (1D)
// **************************************************************************//
inline void Solver1D::Limiter(const int& i){

  // DESCRIPTION
  // ----------------------
  // Compute the slope limiter at point (i) by using the 2 closest neighboring 
  // cells.


  // INPUTS
  // ----------------------
  // i       - index


  // VARIABLE DECLARATION
  // ----------------------
  const static int nb = 2, nr = 2;
  static double Wb[nb], Wr[nr];
  double y, Wmax, Wmin, phi_temp, dWb;


  // 1D LIMITER
  // ----------------------
  for(int n=0; n<nvar; n++){
    Wb[0] = W[i-1][n]; 
    Wb[1] = W[i+1][n];
    Wr[0] = W[i][n] - dWdx[i][n]; 
    Wr[1] = W[i][n] + dWdx[i][n]; 
    Wmax = fmax(Wb[0], Wb[1]);
    Wmin = fmin(Wb[0], Wb[1]);
    phi[i][n] = 1.0; 
    for(int r=0; r<nr; r++){
      dWb = Wr[r] - W[i][n];
      switch (LimiterType){
        case 0: // Zero
          phi[i][n] = 0.0;
          break;
        case 1: // One
          phi[i][n] = 1.0;
          break;
        case 2: // Barth-Jespersen
	  if (dWb>0.0){
            y = fmax(0.0, (Wmax-W[i][n])/dWb);
            phi_temp = fmin(y, 1.0);
	  }else if (dWb<0.0){
            y = fmax(0.0, (Wmin-W[i][n])/dWb);
            phi_temp = fmin(y, 1.0);
          }else{
            phi_temp = 1.0;
          }
          phi[i][n] = fmin(phi[i][n], phi_temp);
          break;
        case 3: // Venkatakrishnan 
          if (dWb>0.0){
            y = (Wmax-W[i][n])/dWb;
            phi_temp = (y*y + 2.0*y)/(y*y + y + 2.0);
          }else if (dWb<0.0){
            y = (Wmin-W[i][n])/dWb;
            phi_temp = (y*y + 2.0*y)/(y*y + y + 2.0);
	  }else{
            phi_temp = 1.0;
          }
          phi[i][n] = fmin(phi[i][n], phi_temp);
          break;
        }
    }
  }
}

// ------------------------------------------------------------------------- //

inline void Solver1D::Limiter(void){
  if(SolverType == 2) return;
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    Limiter(i);
  }
}


// **************************************************************************//
// Limiter (2D)
// **************************************************************************//
inline void Solver2D::Limiter(const int& i, const int& j){

  // DESCRIPTION
  // ----------------------
  // Compute the slope limiter at point (i,j) by using the 8 closest neighboring 
  // cells.


  // INPUTS
  // ----------------------
  // i       - index
  // j       - index


  // VARIABLE DECLARATION
  // ----------------------
  const static int nb = 8, nr = 4;
  static double Wb[nb], Wr[nr];
  double y, Wmax, Wmin, phi_temp, dWb;


  // 2D LIMITER
  // ----------------------
  for(int n=0; n<nvar; n++){
    Wb[0] = W[i-1][j][n]; 
    Wb[1] = W[i+1][j][n];
    Wb[2] = W[i][j-1][n]; 
    Wb[3] = W[i][j+1][n];
    Wb[4] = W[i-1][j-1][n]; 
    Wb[5] = W[i+1][j+1][n];
    Wb[6] = W[i-1][j+1][n]; 
    Wb[7] = W[i+1][j-1][n];
    Wr[0] = W[i][j][n] - dWdx[i][j][n]; 
    Wr[1] = W[i][j][n] + dWdx[i][j][n]; 
    Wr[2] = W[i][j][n] - dWdy[i][j][n]; 
    Wr[3] = W[i][j][n] + dWdy[i][j][n]; 
    Wmax = W[i][j][n]; Wmin=W[i][j][n];
    for(int p=0; p<nb; p++){
      Wmax = fmax(Wmax, Wb[p]);
      Wmin = fmin(Wmin, Wb[p]);
    }
    phi[i][j][n] = 1.0; 
    for(int r=0; r<nr; r++){
      dWb = Wr[r] - W[i][j][n];
      switch (LimiterType){
        case 0: // Zero
          phi[i][j][n] = 0.0;
          break;
        case 1: // One
          phi[i][j][n] = 1.0;
          break;
        case 2: // Barth-Jespersen
	  if (dWb>0.0){
            y = fmax(0.0, (Wmax-W[i][j][n])/dWb);
            phi_temp = fmin(y, 1.0);
	  }else if (dWb<0.0){
            y = fmax(0.0, (Wmin-W[i][j][n])/dWb);
            phi_temp = fmin(y, 1.0);
          }else{
            phi_temp = 1.0;
          }
          phi[i][j][n] = fmin(phi[i][j][n], phi_temp);
          break;
        case 3: // Venkatakrishnan 
          if (dWb>0.0){
            y = (Wmax-W[i][j][n])/dWb;
            phi_temp = (y*y + 2.0*y)/(y*y + y + 2.0);
          }else if (dWb<0.0){
            y = (Wmin-W[i][j][n])/dWb;
            phi_temp = (y*y + 2.0*y)/(y*y + y + 2.0);
	  }else{
            phi_temp = 1.0;
          }
          phi[i][j][n] = fmin(phi[i][j][n], phi_temp);
          break;
        }
    }
  }
}

// ------------------------------------------------------------------------- //

inline void Solver2D::Limiter(void){
  if(SolverType == 2) return;
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      Limiter(i, j);
    }
  }
}


// **************************************************************************//
// Limiter (3D)
// **************************************************************************//
inline void Solver3D::Limiter(const int& i, const int& j, const int& k){

  // DESCRIPTION
  // ----------------------
  // Compute the slope limiter at point (i,j,k) by using the 26 closest
  // neighboring cells.


  // INPUTS
  // ----------------------
  // i       - index
  // j       - index
  // k       - index


  // VARIABLE DECLARATION
  // ----------------------
  const static int nb = 26, nr = 6;
  static double Wb[nb], Wr[nr];
  double y, Wmax, Wmin, phi_temp, dWb;


  // 3D LIMITER
  // ----------------------
  for(int n=0; n<nvar; n++){
    Wb[0] = W[i-1][j][k][n]; 
    Wb[1] = W[i+1][j][k][n];
    Wb[2] = W[i][j-1][k][n]; 
    Wb[3] = W[i][j+1][k][n];
    Wb[4] = W[i][j][k-1][n]; 
    Wb[5] = W[i][j][k+1][n];
    Wb[6] = W[i-1][j-1][k][n]; 
    Wb[7] = W[i+1][j+1][k][n];
    Wb[8] = W[i-1][j+1][k][n]; 
    Wb[9] = W[i+1][j-1][k][n];
    Wb[10] = W[i][j-1][k-1][n]; 
    Wb[11] = W[i][j+1][k+1][n];
    Wb[12] = W[i][j-1][k+1][n]; 
    Wb[13] = W[i][j+1][k-1][n];
    Wb[14] = W[i-1][j][k-1][n]; 
    Wb[15] = W[i+1][j][k+1][n];
    Wb[16] = W[i-1][j][k+1][n]; 
    Wb[17] = W[i+1][j][k-1][n];
    Wb[18] = W[i-1][j-1][k-1][n]; 
    Wb[19] = W[i+1][j+1][k+1][n];
    Wb[20] = W[i-1][j+1][k+1][n]; 
    Wb[21] = W[i+1][j-1][k-1][n];
    Wb[22] = W[i-1][j+1][k-1][n]; 
    Wb[23] = W[i+1][j-1][k+1][n];
    Wb[24] = W[i-1][j-1][k-1][n]; 
    Wb[25] = W[i+1][j+1][k+1][n];
    Wr[0] = W[i][j][k][n] - dWdx[i][j][k][n]; 
    Wr[1] = W[i][j][k][n] + dWdx[i][j][k][n]; 
    Wr[2] = W[i][j][k][n] - dWdy[i][j][k][n]; 
    Wr[3] = W[i][j][k][n] + dWdy[i][j][k][n]; 
    Wr[4] = W[i][j][k][n] - dWdz[i][j][k][n]; 
    Wr[5] = W[i][j][k][n] + dWdz[i][j][k][n]; 
    Wmax = W[i][j][k][n]; Wmin=W[i][j][k][n];
    for(int p=0; p<nb; p++){
      Wmax = fmax(Wmax, Wb[p]);
      Wmin = fmin(Wmin, Wb[p]);
    }
    phi[i][j][k][n] = 1.0; 
    for(int r=0; r<nr; r++){
      dWb = Wr[r] - W[i][j][k][n];
      switch (LimiterType){
        case 0: // Zero
          phi[i][j][k][n] = 0.0;
          break;
        case 1: // One
          phi[i][j][k][n] = 1.0;
          break;
        case 2: // Barth-Jespersen
	  if (dWb>0.0){
            y = fmax(0.0, (Wmax-W[i][j][k][n])/dWb);
            phi_temp = fmin(y, 1.0);
	  }else if (dWb<0.0){
            y = fmax(0.0, (Wmin-W[i][j][k][n])/dWb);
            phi_temp = fmin(y, 1.0);
          }else{
            phi_temp = 1.0;
          }
          phi[i][j][k][n] = fmin(phi[i][j][k][n], phi_temp);
          break;
        case 3: // Venkatakrishnan 
          if (dWb>0.0){
            y = (Wmax-W[i][j][k][n])/dWb;
            phi_temp = (y*y + 2.0*y)/(y*y + y + 2.0);
          }else if (dWb<0.0){
            y = (Wmin-W[i][j][k][n])/dWb;
            phi_temp = (y*y + 2.0*y)/(y*y + y + 2.0);
	  }else{
            phi_temp = 1.0;
          }
          phi[i][j][k][n] = fmin(phi[i][j][k][n], phi_temp);
          break;
        }
    }
  }
}

// ------------------------------------------------------------------------- //

inline void Solver3D::Limiter(void){
  if(SolverType == 2) return;
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
        Limiter(i, j, k);
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// RESET RESIDUAL
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// ResetResidual (1D)
// **************************************************************************//
inline void Solver1D::ResetResidual(const int& index){

  // Reset residual
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      for (int n=0; n<nvar; n++){
        dUdt[K*i+ip][index][n] = 0.0;
      }
    }
  }
}


// **************************************************************************//
// ResetResidual (2D)
// **************************************************************************//
inline void Solver2D::ResetResidual(const int& index){

  // Reset residual
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int ip=0; ip<K; ip++){
        for (int jp=0; jp<K; jp++){
          for (int n=0; n<nvar; n++){
            dUdt[K*i+ip][K*j+jp][index][n] = 0.0;
          }
        }
      }
    }
  }
}


// **************************************************************************//
// ResetResidual (3D)
// **************************************************************************//
inline void Solver3D::ResetResidual(const int& index){

  // Reset residual
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int ip=0; ip<K; ip++){
          for (int jp=0; jp<K; jp++){
            for (int kp=0; kp<K; kp++){
              for (int n=0; n<nvar; n++){
                dUdt[K*i+ip][K*j+jp][K*k+kp][index][n] = U[K*i+ip][K*j+jp][K*k+kp][n];
              }
            }
          }
        }
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// RESIDUAL
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// Residual (1D)
// **************************************************************************//
inline void Solver1D::Residual(const int& index){

  // DESCRIPTION
  // ----------------------
  // Compute the residual for the linear wave equation, the non-linear wave
  // equation, or for euler.


  // INPUTS
  // ----------------------
  // index   - Residual index (0, nResidual-1)


  // RESIDUAL
  // ----------------------
  switch (EquationType){
    case 0: // Test
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          dUdt[K*i+ip][index][0] = U[K*i+ip][0];
        }
      }
      return;
  }

  switch (SolverType){
    case 1: // FVM
      FVMResidual(index);
      break;
    case 2: // FR
      FRResidual(index);
      break;
    default:
      throw invalid_argument("Unknown solver type"); 
  }
}


// **************************************************************************//
// Residual (2D)
// **************************************************************************//
inline void Solver2D::Residual(const int& index){

  // DESCRIPTION
  // ----------------------
  // Compute the residual for the linear wave equation, the non-linear wave
  // equation, or for euler.


  // INPUTS
  // ----------------------
  // index   - Residual index (0, nResidual-1)


  // RESIDUAL
  // ----------------------
  switch (EquationType){
    case 0: // Test
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int j=0; j<mesh.ny+2*mesh.ng; j++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){
              dUdt[K*i+ip][K*j+jp][index][0] = U[K*i+ip][K*j+jp][0];
            }
          }
        }
      }
      return;
  }

  switch (SolverType){
    case 1: // FVM
      FVMResidual(index);
      break;
    case 2: // FR
      FRResidual(index);
      break;
    default:
      throw invalid_argument("Unknown solver type"); 
      break;
  }
}


// **************************************************************************//
// Residual (3D)
// **************************************************************************//
inline void Solver3D::Residual(const int& index){

  // DESCRIPTION
  // ----------------------
  // Compute the residual for the linear wave equation, the non-linear wave
  // equation, or for euler.


  // INPUTS
  // ----------------------
  // index   - Residual index (0, nResidual-1)


  // RESIDUAL
  // ----------------------
  switch (EquationType){
    case 0: // Test
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int j=0; j<mesh.ny+2*mesh.ng; j++){
          for (int k=0; k<mesh.nz+2*mesh.ng; k++){
            for (int ip=0; ip<K; ip++){
              for (int jp=0; jp<K; jp++){
                for (int kp=0; kp<K; kp++){
                  dUdt[K*i+ip][K*j+jp][K*k+kp][index][0] = U[K*i+ip][K*j+jp][K*k+kp][0];
                }
              }
            }
          }
        }
      }
      return;
  }

  switch (SolverType){
    case 1: // FVM
      FVMResidual(index);
      break;
    case 2: // FR
      FRResidual(index);
      break;
    default:
      throw invalid_argument("Unknown solver type"); 
      break;
  }
}


///////////////////////////////////////////////////////////////////////////////
// FVMRESIDUAL
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// FVMResidual (1D)
// **************************************************************************//
inline void Solver1D::FVMResidual(const int& index){

  // Compute boundary solution
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int n=0; n<nvar; n++){
      Wl[i][n] = W[i][n]-phi[i][n]*dWdx[i][n];
      Wr[i][n] = W[i][n]+phi[i][n]*dWdx[i][n];
    }
  }

  // Compute boundary upwind flux
  for (int i=0; i<mesh.nx+mesh.ng; i++){
    if (EquationType==1){ // Linear advection
      if(wavespeed>0.0) Fupw[i][0] = wavespeed*Wr[i][0];
      else              Fupw[i][0] = wavespeed*Wl[i+1][0];
    }else if (EquationType==2){ // Euler
      Euler1D.RoeSolver(&Wr[i][0], &Wl[i+1][0], &Fupw[i][0]);
    }
  }

  // Compute residual
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int n=0; n<nvar; n++){
      dUdt[i][index][n] = -(Fupw[i][n]-Fupw[i-1][n])/mesh.dx[i];
    }
  }
}


// **************************************************************************//
// FVMResidual (2D)
// **************************************************************************//
inline void Solver2D::FVMResidual(const int& index){

  // X-DIRECTION
  // Compute boundary solution
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int n=0; n<nvar; n++){
        Wl[i][j][n] = W[i][j][n]-phi[i][j][n]*dWdx[i][j][n];
        Wr[i][j][n] = W[i][j][n]+phi[i][j][n]*dWdx[i][j][n];
      }
    }
  }

  // Compute boundary upwind flux
  for (int i=0; i<mesh.nx+mesh.ng; i++){
    for (int j=0; j<mesh.ny+mesh.ng; j++){
      if (EquationType==1){
        if(wavespeed>0.0) Fupw[i][j][0] = wavespeed*Wr[i][j][0];
        else              Fupw[i][j][0] = wavespeed*Wl[i+1][j][0];
      }else if (EquationType==2){
        Euler2D.RoeSolver(&Wr[i][j][0], &Wl[i+1][j][0], &Fupw[i][j][0]);
      }
    }
  }

  // Compute residual
  double xl, xr;
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int n=0; n<nvar; n++){
        dUdt[i][j][index][n] = -(Fupw[i][j][n]-Fupw[i-1][j][n])/mesh.dx[i][j];
      }
    }
  }

  // Y-DIRECTION
  // Compute boundary solution
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int n=0; n<nvar; n++){
        Wl[i][j][n] = W[i][j][n]-phi[i][j][n]*dWdy[i][j][n];
        Wr[i][j][n] = W[i][j][n]+phi[i][j][n]*dWdy[i][j][n];
      }
    }
  }

  // Compute boundary upwind flux
  double ur, vr, ul, vl, Fu, Fv;
  for (int i=0; i<mesh.nx+mesh.ng; i++){
    for (int j=0; j<mesh.ny+mesh.ng; j++){
      if (EquationType==1){
        if(wavespeed>0.0) Fupw[i][j][0] = wavespeed*Wr[i][j][0];
        else              Fupw[i][j][0] = wavespeed*Wl[i][j+1][0];
      }else if (EquationType==2){
        ur = Wr[i][j][1]; vr = Wr[i][j][2]; ul = Wl[i][j+1][1]; vl = Wl[i][j+1][2];
        Wr[i][j][1] = vr; Wr[i][j][2] = -ur; Wl[i][j+1][1] = vl; Wl[i][j+1][2] = -ul;
        Euler2D.RoeSolver(&Wr[i][j][0], &Wl[i][j+1][0], &Fupw[i][j][0]);
        Fu = Fupw[i][j][1]; Fv = Fupw[i][j][2];
        Fupw[i][j][1] = -Fv; Fupw[i][j][2] = Fu;
      }
    }
  }

  // Compute residual
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int n=0; n<nvar; n++){
        dUdt[i][j][index][n] += -(Fupw[i][j][n]-Fupw[i][j-1][n])/mesh.dy[i][j];
      }
    }
  }
}


// **************************************************************************//
// FVMResidual (3D)
// **************************************************************************//
inline void Solver3D::FVMResidual(const int& index){

  // X-DIRECTION
  // Compute boundary solution
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int n=0; n<nvar; n++){
          Wl[i][j][k][n] = W[i][j][k][n]-phi[i][j][k][n]*dWdx[i][j][k][n];
          Wr[i][j][k][n] = W[i][j][k][n]+phi[i][j][k][n]*dWdx[i][j][k][n];
        }
      }
    }
  }

  // Compute boundary upwind flux
  for (int i=0; i<mesh.nx+mesh.ng; i++){
    for (int j=0; j<mesh.ny+mesh.ng; j++){
      for (int k=0; k<mesh.nz+mesh.ng; k++){
        if (EquationType==1){
          if(wavespeed>0.0) Fupw[i][j][k][0] = wavespeed*Wr[i][j][k][0];
          else              Fupw[i][j][k][0] = wavespeed*Wl[i+1][j][k][0];
        }else if (EquationType==2){
          Euler3D.RoeSolver(&Wr[i][j][k][0], &Wl[i+1][j][k][0], &Fupw[i][j][k][0]);
        }
      }
    }
  }

  // Compute residual
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
        for (int n=0; n<nvar; n++){
          dUdt[i][j][k][index][n] = -(Fupw[i][j][k][n]-Fupw[i-1][j][k][n])/mesh.dx[i][j][k];
        }
      }
    }
  }

  // Y-DIRECTION
  // Compute boundary solution
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int n=0; n<nvar; n++){
          Wl[i][j][k][n] = W[i][j][k][n]-phi[i][j][k][n]*dWdy[i][j][k][n];
          Wr[i][j][k][n] = W[i][j][k][n]+phi[i][j][k][n]*dWdy[i][j][k][n];
        }
      }
    }
  }

  // Compute boundary upwind flux
  double ur, vr, wr, ul, vl, wl, Fu, Fv, Fw;
  for (int i=0; i<mesh.nx+mesh.ng; i++){
    for (int j=0; j<mesh.ny+mesh.ng; j++){
      for (int k=0; k<mesh.nz+mesh.ng; k++){
        if (EquationType==1){
          if(wavespeed>0.0) Fupw[i][j][k][0] = wavespeed*Wr[i][j][k][0];
          else              Fupw[i][j][k][0] = wavespeed*Wl[i][j+1][k][0];
        }else if (EquationType==2){
          ur = Wr[i][j][k][1]; vr = Wr[i][j][k][2]; ul = Wl[i][j+1][k][1]; vl = Wl[i][j+1][k][2];
          Wr[i][j][k][1] = vr; Wr[i][j][k][2] =-ur; Wl[i][j+1][k][1] = vl; Wl[i][j+1][k][2] =-ul;
          Euler3D.RoeSolver(&Wr[i][j][k][0], &Wl[i][j+1][k][0], &Fupw[i][j][k][0]);
          Fu = Fupw[i][j][k][1]; Fv = Fupw[i][j][k][2];
          Fupw[i][j][k][1] =-Fv; Fupw[i][j][k][2] = Fu;
        }
      }
    }
  }

  // Compute residual
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
        for (int n=0; n<nvar; n++){
          dUdt[i][j][k][index][n] += -(Fupw[i][j][k][n]-Fupw[i][j-1][k][n])/mesh.dy[i][j][k];
        }
      }
    }
  }

  // Z-DIRECTION
  // Compute boundary solution
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int n=0; n<nvar; n++){
          Wl[i][j][k][n] = W[i][j][k][n]-phi[i][j][k][n]*dWdz[i][j][k][n];
          Wr[i][j][k][n] = W[i][j][k][n]+phi[i][j][k][n]*dWdz[i][j][k][n];
        }
      }
    }
  }

  // Compute boundary upwind flux
  for (int i=0; i<mesh.nx+mesh.ng; i++){
    for (int j=0; j<mesh.ny+mesh.ng; j++){
      for (int k=0; k<mesh.nz+mesh.ng; k++){
        if (EquationType==1){
          if(wavespeed>0.0) Fupw[i][j][k][0] = wavespeed*Wr[i][j][k][0];
          else              Fupw[i][j][k][0] = wavespeed*Wl[i][j][k+1][0];
        }else if (EquationType==2){
          ur = Wr[i][j][k][1]; wr = Wr[i][j][k][3]; ul = Wl[i][j][k+1][1]; wl = Wl[i][j][k+1][3];
          Wr[i][j][k][1] = wr; Wr[i][j][k][3] =-ur; Wl[i][j][k+1][1] = wl; Wl[i][j][k+1][3] =-ul;
          Euler3D.RoeSolver(&Wr[i][j][k][0], &Wl[i][j][k+1][0], &Fupw[i][j][k][0]);
          Fu = Fupw[i][j][k][1]; Fw = Fupw[i][j][k][3];
          Fupw[i][j][k][1] =-Fw; Fupw[i][j][k][3] = Fu;
        }
      }
    }
  }

  // Compute residual
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
        for (int n=0; n<nvar; n++){
          dUdt[i][j][k][index][n] += -(Fupw[i][j][k][n]-Fupw[i][j][k-1][n])/mesh.dz[i][j][k];
        }
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// FRRESIDUAL
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// FRResidual (1D)
// **************************************************************************//
inline void Solver1D::FRResidual(const int& index){

  // Variable declaration
  double Flux[K][nvar];
  static double* Fl = new double[nvar];
  static double* Fr = new double[nvar];
  static double* Flcom = new double[nvar];
  static double* Frcom = new double[nvar];

  // X-DIRECTION
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    // Compute solution flux
    for (int ip=0; ip<K; ip++){
      if (EquationType==1){
        Flux[ip][0] = wavespeed*U[K*i+ip][0];
      }else if (EquationType==2){
        Euler1D.InviscidFlux(&U[K*i+ip][0], &Flux[ip][0]);
      }
    }

    // Compute boundary flux
    for (int n=0; n<nvar; n++){
      Fl[n] = 0.0;
      Fr[n] = 0.0;
      for (int ip=0; ip<K; ip++){
        Fl[n] += FR.Lvec[ip]*Flux[ip][n];
        Fr[n] += FR.Rvec[ip]*Flux[ip][n];
      }
    }

    // Compute common flux
    if (i==mesh.ng) CommonFluxXdir(Frcom, i-1);
    for (int n=0; n<nvar; n++) Flcom[n] = Frcom[n];
    CommonFluxXdir(Frcom, i);

    // Compute residual
    for (int ip=0; ip<K; ip++){
      for (int n=0; n<nvar; n++){
        // Reset residual
        dUdt[K*i+ip][index][n] = 0.0;
        // Uncorrected flux
        for (int s=0; s<K; s++){
          dUdt[K*i+ip][index][n] += FR.Dmat[ip][s]*Flux[s][n];
        }
        // Flux correction
        dUdt[K*i+ip][index][n] += (Flcom[n]-Fl[n])*FR.gl_prime[ip] + 
                                  (Frcom[n]-Fr[n])*FR.gr_prime[ip];
        // Cell Scaling
        dUdt[K*i+ip][index][n] *= -2.0/mesh.dx[i];
      } 
    }
  }
}

// **************************************************************************//
// CommonFluxXdir
// **************************************************************************//
inline void Solver1D::CommonFluxXdir(double* commonflux, int i){

  // Variable declaration
  static double* Wl = new double[nvar];
  static double* Wr = new double[nvar];
  static double* Ul = new double[nvar];
  static double* Ur = new double[nvar];

  // Compute boundary solution
  for (int n=0; n<nvar; n++){
    Ul[n] = 0.0;
    Ur[n] = 0.0;
    for (int ip=0; ip<K; ip++){
      Ul[n] += FR.Rvec[ip]*U[K*i+ip][n];
      Ur[n] += FR.Lvec[ip]*U[K*(i+1)+ip][n];
    }
  }

  // Update primitive variable
  if (EquationType==2){
    Euler1D.ConstoPrim(&Ul[0],&Wl[0]);
    Euler1D.ConstoPrim(&Ur[0],&Wr[0]);
  }

  // Compute upwind flux
  if (EquationType==1){
    if(wavespeed>0.0) commonflux[0] = wavespeed*Ul[0];
    else              commonflux[0] = wavespeed*Ur[0];
  }else if (EquationType==2){
    Euler1D.RoeSolver(&Wl[0], &Wr[0], &commonflux[0]);
  }
}


// **************************************************************************//
// FRResidual (2D)
// **************************************************************************//
inline void Solver2D::FRResidual(const int& index){

  // Variable declaration
  double u, v, ul, vl, ur, vr, Fu, Fv;
  double Flux[K][nvar];
  static double* Fl = new double[nvar];
  static double* Fr = new double[nvar];
  static double* Flcom = new double[nvar];
  static double* Frcom = new double[nvar];

  // X-DIRECTION
  for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
    for (int jp=0; jp<K; jp++){
      // Compute solution flux
      for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (EquationType==1){
            Flux[ip][0] = wavespeed*U[K*i+ip][K*j+jp][0];
          }else if (EquationType==2){
            Euler2D.InviscidFlux(&U[K*i+ip][K*j+jp][0], &Flux[ip][0]);
          }
        }

        // Compute boundary flux
        for (int n=0; n<nvar; n++){
          Fl[n] = 0.0;
          Fr[n] = 0.0;
          for (int ip=0; ip<K; ip++){
            Fl[n] += FR.Lvec[ip]*Flux[ip][n];
            Fr[n] += FR.Rvec[ip]*Flux[ip][n];
          }
        }

        // Compute common flux
        if (i==mesh.ng) CommonFluxXdir(Frcom, i-1, j, jp);
        for (int n=0; n<nvar; n++) Flcom[n] = Frcom[n];
        CommonFluxXdir(Frcom, i, j, jp);

        // Compute residual
        for (int ip=0; ip<K; ip++){
          for (int n=0; n<nvar; n++){
            // Reset residual
            dUdt[K*i+ip][K*j+jp][index][n] = 0.0;
            // Uncorrected flux
            for (int s=0; s<K; s++){
              dUdt[K*i+ip][K*j+jp][index][n] -= 2.0/mesh.dx[i][j]*FR.Dmat[ip][s]*Flux[s][n];
            }
            // Flux correction
            dUdt[K*i+ip][K*j+jp][index][n] -= 2.0/mesh.dx[i][j]*((Flcom[n]-Fl[n])*FR.gl_prime[ip] + 
                                                                 (Frcom[n]-Fr[n])*FR.gr_prime[ip]);
          }
        }
      }
    }
  }

  // Y-DIRECTION
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      // Compute solution flux
      for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
        for (int jp=0; jp<K; jp++){
          if (EquationType==1){
            Flux[jp][0] = wavespeed*U[K*i+ip][K*j+jp][0];
          }else if (EquationType==2){
            u = U[K*i+ip][K*j+jp][1]; v = U[K*i+ip][K*j+jp][2];
            U[K*i+ip][K*j+jp][1] = v; U[K*i+ip][K*j+jp][2] = -u;
            Euler2D.InviscidFlux(&U[K*i+ip][K*j+jp][0], &Flux[jp][0]);
            Fu = Flux[jp][1];  Fv = Flux[jp][2];
            Flux[jp][1] = -Fv; Flux[jp][2] = Fu;
            U[K*i+ip][K*j+jp][1] = u;   U[K*i+ip][K*j+jp][2] = v;
          }
        }

        // Compute boundary flux
        for (int n=0; n<nvar; n++){
          Fl[n] = 0.0;
          Fr[n] = 0.0;
          for (int jp=0; jp<K; jp++){
            Fl[n] += FR.Lvec[jp]*Flux[jp][n];
            Fr[n] += FR.Rvec[jp]*Flux[jp][n];
          }
        }

        // Compute common flux
        if (j==mesh.ng) CommonFluxYdir(Frcom, i, j-1, ip);
        for (int n=0; n<nvar; n++) Flcom[n] = Frcom[n];
        CommonFluxYdir(Frcom, i, j, ip);

        // Compute residual
        for (int jp=0; jp<K; jp++){
          for (int n=0; n<nvar; n++){
            // Uncorrected flux
            for (int s=0; s<K; s++){
              dUdt[K*i+ip][K*j+jp][index][n] -= 2.0/mesh.dy[i][j]*FR.Dmat[jp][s]*Flux[s][n];
            }
            // Flux correction
            dUdt[K*i+ip][K*j+jp][index][n] -= 2.0/mesh.dy[i][j]*((Flcom[n]-Fl[n])*FR.gl_prime[jp] + 
                                                                 (Frcom[n]-Fr[n])*FR.gr_prime[jp]);
          }
        }
      }
    }
  }
}

// **************************************************************************//
// CommonFluxXdir
// **************************************************************************//
inline void Solver2D::CommonFluxXdir(double* commonflux, int i, int j, int jp){

  // Variable declaration
  static double* Wl = new double[nvar];
  static double* Wr = new double[nvar];
  static double* Ul = new double[nvar];
  static double* Ur = new double[nvar];

  // Compute boundary solution
  for (int n=0; n<nvar; n++){
    Ul[n] = 0.0;
    Ur[n] = 0.0;
    for (int ip=0; ip<K; ip++){
      Ul[n] += FR.Rvec[ip]*U[K*i+ip][K*j+jp][n];
      Ur[n] += FR.Lvec[ip]*U[K*(i+1)+ip][K*j+jp][n];
    }
  }

  // Update primitive variable
  if (EquationType==2){
    Euler2D.ConstoPrim(&Ul[0],&Wl[0]);
    Euler2D.ConstoPrim(&Ur[0],&Wr[0]);
  }

  // Compute upwind flux
  if (EquationType==1){
    if(wavespeed>0.0) commonflux[0] = wavespeed*Ul[0];
    else              commonflux[0] = wavespeed*Ur[0];
  }else if (EquationType==2){
    Euler2D.RoeSolver(&Wl[0], &Wr[0], &commonflux[0]);
  }
}

// **************************************************************************//
// CommonFluxYdir
// **************************************************************************//
inline void Solver2D::CommonFluxYdir(double* commonflux, int i, int j, int ip){

  // Variable declaration
  double u, v, ul, vl, ur, vr, Fu, Fv;
  static double* Wl = new double[nvar];
  static double* Wr = new double[nvar];
  static double* Ul = new double[nvar];
  static double* Ur = new double[nvar];

  // Compute boundary solution
  for (int n=0; n<nvar; n++){
    Ul[n] = 0.0;
    Ur[n] = 0.0;
    for (int jp=0; jp<K; jp++){
      Ul[n] += FR.Rvec[jp]*U[K*i+ip][K*j+jp][n];
      Ur[n] += FR.Lvec[jp]*U[K*i+ip][K*(j+1)+jp][n];
    }
  }

  // Update primitive variable
  if (EquationType==2){
    Euler2D.ConstoPrim(&Ul[0],&Wl[0]);
    Euler2D.ConstoPrim(&Ur[0],&Wr[0]);
  }

  // Compute upwind flux
  if (EquationType==1){
    if(wavespeed>0.0) commonflux[0] = wavespeed*Ul[0];
    else              commonflux[0] = wavespeed*Ur[0];
  }else if (EquationType==2){
    ur = Wr[1]; vr = Wr[2]; ul = Wl[1]; vl = Wl[2];
    Wr[1] = vr; Wr[2] = -ur; Wl[1] = vl; Wl[2] = -ul;
    Euler2D.RoeSolver(&Wl[0], &Wr[0], &commonflux[0]);
    Fu = commonflux[1];  Fv = commonflux[2];
    commonflux[1] = -Fv; commonflux[2] = Fu;
  }
}


// **************************************************************************//
// FRResidual (3D)
// **************************************************************************//
inline void Solver3D::FRResidual(const int& index){

  // Variable declaration
  double u, v, w, ul, vl, wl, ur, vr, wr, Fu, Fv, Fw;
  double Flux[K][nvar];
  static double* Fl = new double[nvar];
  static double* Fr = new double[nvar];
  static double* Flcom = new double[nvar];
  static double* Frcom = new double[nvar];

  // X-DIRECTION
  for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
    for (int kp=0; kp<K; kp++){
      for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
        for (int jp=0; jp<K; jp++){
          // Compute solution flux
          for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
            for (int ip=0; ip<K; ip++){
              if (EquationType==1){
                Flux[ip][0] = wavespeed*U[K*i+ip][K*j+jp][K*k+kp][0];
              }else if (EquationType==2){
                Euler3D.InviscidFlux(&U[K*i+ip][K*j+jp][K*k+kp][0], &Flux[ip][0]);
              }
            }

            // Compute boundary flux
            for (int n=0; n<nvar; n++){
              Fl[n] = 0.0;
              Fr[n] = 0.0;
              for (int ip=0; ip<K; ip++){
                Fl[n] += FR.Lvec[ip]*Flux[ip][n];
                Fr[n] += FR.Rvec[ip]*Flux[ip][n];
              }
            }

            // Compute common flux
            if (i==mesh.ng) CommonFluxXdir(Frcom, i-1, j, k, jp, kp);
            for (int n=0; n<nvar; n++) Flcom[n] = Frcom[n];
            CommonFluxXdir(Frcom, i, j, k, jp, kp);

            // Compute residual
            for (int ip=0; ip<K; ip++){
              for (int n=0; n<nvar; n++){
                // Reset residual
                dUdt[K*i+ip][K*j+jp][K*k+kp][index][n] = 0.0;
                // Uncorrected flux
                for (int s=0; s<K; s++){
                  dUdt[K*i+ip][K*j+jp][K*k+kp][index][n] -= 2.0/mesh.dx[i][j][k]*FR.Dmat[ip][s]*Flux[s][n];
                }
                // Flux correction
                dUdt[K*i+ip][K*j+jp][K*k+kp][index][n] -= 2.0/mesh.dx[i][j][k]*((Flcom[n]-Fl[n])*FR.gl_prime[ip] + 
                                                                                (Frcom[n]-Fr[n])*FR.gr_prime[ip]);
              }
            }
          }
        }
      }
    }
  }

  // Y-DIRECTION
  for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
    for (int kp=0; kp<K; kp++){
      for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          // Compute solution flux
          for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
            for (int jp=0; jp<K; jp++){
              if (EquationType==1){
                Flux[jp][0] = wavespeed*U[K*i+ip][K*j+jp][K*k+kp][0];
              }else if (EquationType==2){
                u = U[K*i+ip][K*j+jp][K*k+kp][1]; v = U[K*i+ip][K*j+jp][K*k+kp][2];
                U[K*i+ip][K*j+jp][K*k+kp][1] = v; U[K*i+ip][K*j+jp][K*k+kp][2] = -u;
                Euler3D.InviscidFlux(&U[K*i+ip][K*j+jp][K*k+kp][0], &Flux[jp][0]);
                Fu = Flux[jp][1]; Fv = Flux[jp][2];
                Flux[jp][1] = -Fv; Flux[jp][2] = Fu;
                U[K*i+ip][K*j+jp][K*k+kp][1] = u;   U[K*i+ip][K*j+jp][K*k+kp][2] = v;
              }
            }

            // Compute boundary flux
            for (int n=0; n<nvar; n++){
              Fl[n] = 0.0;
              Fr[n] = 0.0;
              for (int jp=0; jp<K; jp++){
                Fl[n] += FR.Lvec[jp]*Flux[jp][n];
                Fr[n] += FR.Rvec[jp]*Flux[jp][n];
              }
            }

            // Compute common flux
            if (j==mesh.ng) CommonFluxYdir(Frcom, i, j-1, k, ip, kp);
            for (int n=0; n<nvar; n++) Flcom[n] = Frcom[n];
            CommonFluxYdir(Frcom, i, j, k, ip, kp);

            // Compute residual
            for (int jp=0; jp<K; jp++){
              for (int n=0; n<nvar; n++){
                // Uncorrected flux
                for (int s=0; s<K; s++){
                  dUdt[K*i+ip][K*j+jp][K*k+kp][index][n] -= 2.0/mesh.dy[i][j][k]*FR.Dmat[jp][s]*Flux[s][n];
                }
                // Flux correction
                dUdt[K*i+ip][K*j+jp][K*k+kp][index][n] -= 2.0/mesh.dy[i][j][k]*((Flcom[n]-Fl[n])*FR.gl_prime[jp] + 
                                                                                (Frcom[n]-Fr[n])*FR.gr_prime[jp]);
              }
            }
          }
        }
      }
    }
  }

  // Z-DIRECTION
  for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
    for (int jp=0; jp<K; jp++){
      for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          // Compute solution flux
          for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
            for (int kp=0; kp<K; kp++){
              if (EquationType==1){
                Flux[kp][0] = wavespeed*U[K*i+ip][K*j+jp][K*k+kp][0];
              }else if (EquationType==2){
                u = U[K*i+ip][K*j+jp][K*k+kp][1]; w = U[K*i+ip][K*j+jp][K*k+kp][3];
                U[K*i+ip][K*j+jp][K*k+kp][1] = w; U[K*i+ip][K*j+jp][K*k+kp][3] = -u;
                Euler3D.InviscidFlux(&U[K*i+ip][K*j+jp][K*k+kp][0], &Flux[kp][0]);
                Fu = Flux[kp][1];  Fw = Flux[kp][3];
                Flux[kp][1] = -Fw; Flux[kp][3] = Fu;
                U[K*i+ip][K*j+jp][K*k+kp][1] = u;   U[K*i+ip][K*j+jp][K*k+kp][3] = w;
              }
            }

            // Compute boundary flux
            for (int n=0; n<nvar; n++){
              Fl[n] = 0.0;
              Fr[n] = 0.0;
              for (int kp=0; kp<K; kp++){
                Fl[n] += FR.Lvec[kp]*Flux[kp][n];
                Fr[n] += FR.Rvec[kp]*Flux[kp][n];
              }
            }

            // Compute common flux
            if (k==mesh.ng) CommonFluxZdir(Frcom, i, j, k-1, ip, jp);
            for (int n=0; n<nvar; n++) Flcom[n] = Frcom[n];
            CommonFluxZdir(Frcom, i, j, k, ip, jp);

            // Compute residual
            for (int kp=0; kp<K; kp++){
              for (int n=0; n<nvar; n++){
                // Uncorrected flux
                for (int s=0; s<K; s++){
                  dUdt[K*i+ip][K*j+jp][K*k+kp][index][n] -= 2.0/mesh.dz[i][j][k]*FR.Dmat[kp][s]*Flux[s][n];
                }
                // Flux correction
                dUdt[K*i+ip][K*j+jp][K*k+kp][index][n] -= 2.0/mesh.dz[i][j][k]*((Flcom[n]-Fl[n])*FR.gl_prime[kp] + 
                                                                                (Frcom[n]-Fr[n])*FR.gr_prime[kp]);
              }
            }
          }
        }
      }
    }
  }
}

// **************************************************************************//
// CommonFluxXdir
// **************************************************************************//
inline void Solver3D::CommonFluxXdir(double* commonflux, int i, int j, int k, int jp, int kp){

  // Variable declaration
  static double* Wl = new double[nvar];
  static double* Wr = new double[nvar];
  static double* Ul = new double[nvar];
  static double* Ur = new double[nvar];

  // Compute boundary solution
  for (int n=0; n<nvar; n++){
    Ul[n] = 0.0;
    Ur[n] = 0.0;
    for (int ip=0; ip<K; ip++){
      Ul[n] += FR.Rvec[ip]*U[K*i+ip][K*j+jp][K*k+kp][n];
      Ur[n] += FR.Lvec[ip]*U[K*(i+1)+ip][K*j+jp][K*k+kp][n];
    }
  }

  // Update primitive variable
  if (EquationType==2){
    Euler3D.ConstoPrim(&Ul[0],&Wl[0]);
    Euler3D.ConstoPrim(&Ur[0],&Wr[0]);
  }

  // Compute upwind flux
  if (EquationType==1){
    if(wavespeed>0.0) commonflux[0] = wavespeed*Ul[0];
    else              commonflux[0] = wavespeed*Ur[0];
  }else if (EquationType==2){
    Euler3D.RoeSolver(&Wl[0], &Wr[0], &commonflux[0]);
  }
}

// **************************************************************************//
// CommonFluxYdir
// **************************************************************************//
inline void Solver3D::CommonFluxYdir(double* commonflux, int i, int j, int k, int ip, int kp){

  // Variable declaration
  double u, v, w, ul, vl, wl, ur, vr, wr, Fu, Fv, Fw;
  static double* Wl = new double[nvar];
  static double* Wr = new double[nvar];
  static double* Ul = new double[nvar];
  static double* Ur = new double[nvar];

  // Compute boundary solution
  for (int n=0; n<nvar; n++){
    Ul[n] = 0.0;
    Ur[n] = 0.0;
    for (int jp=0; jp<K; jp++){
      Ul[n] += FR.Rvec[jp]*U[K*i+ip][K*j+jp][K*k+kp][n];
      Ur[n] += FR.Lvec[jp]*U[K*i+ip][K*(j+1)+jp][K*k+kp][n];
    }
  }

  // Update primitive variable
  if (EquationType==2){
    Euler3D.ConstoPrim(&Ul[0],&Wl[0]);
    Euler3D.ConstoPrim(&Ur[0],&Wr[0]);
  }

  // Compute upwind flux
  if (EquationType==1){
    if(wavespeed>0.0) commonflux[0] = wavespeed*Ul[0];
    else              commonflux[0] = wavespeed*Ur[0];
  }else if (EquationType==2){
    ur = Wr[1]; vr = Wr[2]; ul = Wl[1]; vl = Wl[2];
    Wr[1] = vr; Wr[2] = -ur; Wl[1] = vl; Wl[2] = -ul;
    Euler3D.RoeSolver(&Wl[0], &Wr[0], &commonflux[0]);
    Fu = commonflux[1];  Fv = commonflux[2];
    commonflux[1] = -Fv; commonflux[2] = Fu;
  }
}

// **************************************************************************//
// CommonFluxZdir
// **************************************************************************//
inline void Solver3D::CommonFluxZdir(double* commonflux, int i, int j, int k, int ip, int jp){

  // Variable declaration
  double u, v, w, ul, vl, wl, ur, vr, wr, Fu, Fv, Fw;
  static double* Wl = new double[nvar];
  static double* Wr = new double[nvar];
  static double* Ul = new double[nvar];
  static double* Ur = new double[nvar];

  // Compute boundary solution
  for (int n=0; n<nvar; n++){
    Ul[n] = 0.0;
    Ur[n] = 0.0;
    for (int kp=0; kp<K; kp++){
      Ul[n] += FR.Rvec[kp]*U[K*i+ip][K*j+jp][K*k+kp][n];
      Ur[n] += FR.Lvec[kp]*U[K*i+ip][K*j+jp][K*(k+1)+kp][n];
    }
  }

  // Update primitive variable
  if (EquationType==2){
    Euler3D.ConstoPrim(&Ul[0],&Wl[0]);
    Euler3D.ConstoPrim(&Ur[0],&Wr[0]);
  }

  // Compute upwind flux
  if (EquationType==1){
    if(wavespeed>0.0) commonflux[0] = wavespeed*Ul[0];
    else              commonflux[0] = wavespeed*Ur[0];
  }else if (EquationType==2){
    ur = Wr[1]; wr = Wr[3]; ul = Wl[1]; wl = Wl[3];
    Wr[1] = wr; Wr[3] = -ur; Wl[1] = wl; Wl[3] = -ul;
    Euler3D.RoeSolver(&Wl[0], &Wr[0], &commonflux[0]);
    Fu = commonflux[1];  Fw = commonflux[3];
    commonflux[1] = -Fw; commonflux[3] = Fu;
  }
}


///////////////////////////////////////////////////////////////////////////////
// TIMESTEPPING
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// AdamsBashforth (1D) - not working well
// **************************************************************************//
/*inline void Solver1D::AdamsBashforth(const int& nstep){

  // DESCRIPTION
  // ----------------------
  // Semi-discrete Adams-Bashforth multi-step time-marching approach for 1D 
  // problems.


  // INPUTS
  // ----------------------
  // nstep   - Number of steps (1-5) of the Adams-Bashforth method


  // ADAMS-BASHFORTH
  // ----------------------
  // Enforce BCs
  ApplyBCs();

  // Update primitive variables
  UpdatePrimitive();

  // Enforce positivity
  EnforcePositivity();

  // Update gradient
  Gradient();

  // Limiter
  Limiter();

  // Reset residual
  ResetResidual(0);

  // Update residual
  this->Residual(0);

  // Update solution
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      for (int n=0; n<nvar; n++){
        switch (nstep){
          case 1: // AB1
            U[K*i+ip][n] += dt[K*i+ip]*dUdt[K*i+ip][0][n];
            break;
          case 2: // AB2
            U[K*i+ip][n] += dt[K*i+ip]/2.0*(3.0*dUdt[K*i+ip][0][n] 
                                          - 1.0*dUdt[K*i+ip][1][n]);
            break;
          case 3: // AB3
            U[K*i+ip][n] += dt[K*i+ip]/12.0*(23.0*dUdt[K*i+ip][0][n] 
                                           - 16.0*dUdt[K*i+ip][1][n] 
                                           +  5.0*dUdt[K*i+ip][2][n]);
            break;
          case 4: // AB4
            U[K*i+ip][n] += dt[K*i+ip]/24.0*(55.0*dUdt[K*i+ip][0][n]
                                           - 59.0*dUdt[K*i+ip][1][n]
                                           + 37.0*dUdt[K*i+ip][2][n] 
                                           -  9.0*dUdt[K*i+ip][3][n]);
            break;
          case 5: // AB5
            U[K*i+ip][n] += dt[K*i+ip]/720.0*(1901.0*dUdt[K*i+ip][0][n] 
                                            - 2774.0*dUdt[K*i+ip][1][n] 
                                            + 2616.0*dUdt[K*i+ip][2][n] 
                                            - 1274.0*dUdt[K*i+ip][3][n]
                                            +  251.0*dUdt[K*i+ip][4][n]);
            break;
          default:
            throw invalid_argument("Step number out of range"); 
        }
        //for (int s=nstep-1; s>0; s--){
        for (int s=nstep; s>0; s--){
          dUdt[K*i+ip][s][n] = dUdt[K*i+ip][s-1][n];
        }
      }
    }
  }
}*/


// **************************************************************************//
// RungeKutta (1D)
// **************************************************************************//
inline void Solver1D::RungeKutta(const int& nstage){

  // DESCRIPTION
  // ----------------------
  // Semi-discrete Runge-Kutta multi-stage time-marching approach for 1D 
  // problems.


  // INPUTS
  // ----------------------
  // nstage  - Number of stages (1-6) of the Runge-Kutta method


  // RUNGKE KUTTA
  // ----------------------
  // Enforce BCs
  ApplyBCs();

  // Update primitive variables
  UpdatePrimitive();

  // Enforce positivity
  EnforcePositivity();

  // Stage residual
  for (int s=0; s<nstage; s++){
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      for (int ip=0; ip<K; ip++){
        for (int n=0; n<nvar; n++){
          if (s==0) Un[K*i+ip][n] = U[K*i+ip][n];
          // RK1
          if (nstage==1){
          // RK2
          }else if (nstage==2){
            if (s==1)      U[K*i+ip][n] = Un[K*i+ip][n] 
                                         + dt[K*i+ip]*dUdt[K*i+ip][s-1][n];
          // RK3
          }else if (nstage==3){
            if (s==1)      U[K*i+ip][n] = Un[K*i+ip][n] 
                                        + dt[K*i+ip]/3.0*dUdt[K*i+ip][s-1][n];
            else if (s==2) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 2.0*dt[K*i+ip]/3.0*dUdt[K*i+ip][s-1][n];
          // RK4
          }else if (nstage==4){
            if (s==1)      U[K*i+ip][n] = Un[K*i+ip][n] 
                                        + dt[K*i+ip]/2.0*dUdt[K*i+ip][s-1][n];
            else if (s==2) U[K*i+ip][n] = Un[K*i+ip][n] 
                                        + dt[K*i+ip]/2.0*dUdt[K*i+ip][s-1][n];
            else if (s==3) U[K*i+ip][n] = Un[K*i+ip][n]
                                        + dt[K*i+ip]*dUdt[K*i+ip][s-1][n];
          // RKM4
          }else if (nstage==5){
            if (s==1)      U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 1.0*dt[K*i+ip]/3.0*dUdt[K*i+ip][s-1][n];
            else if (s==2) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 1.0*dt[K*i+ip]/6.0*dUdt[K*i+ip][s-1][n] 
                                    + 1.0*dt[K*i+ip]/6.0*dUdt[K*i+ip][s-2][n];
            else if (s==3) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 3.0*dt[K*i+ip]/8.0*dUdt[K*i+ip][s-1][n] 
                                    + 1.0*dt[K*i+ip]/8.0*dUdt[K*i+ip][s-3][n];
            else if (s==4) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 2.0*dt[K*i+ip]/1.0*dUdt[K*i+ip][s-1][n] 
                                    - 3.0*dt[K*i+ip]/2.0*dUdt[K*i+ip][s-2][n]
                                    + 1.0*dt[K*i+ip]/2.0*dUdt[K*i+ip][s-4][n];
          // RKF5 - Formula I
          /*}else if (nstage==6){
            if (s==1)      U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 2.0*dt[K*i+ip]/9.0*dUdt[K*i+ip][s-1][n];
            else if (s==2) U[K*i+ip][n] = Un[K*i+ip][n]
                                    + 1.0*dt[K*i+ip]/4.0*dUdt[K*i+ip][s-1][n] 
                                    + 1.0*dt[K*i+ip]/12.0*dUdt[K*i+ip][s-2][n];
            else if (s==3) U[K*i+ip][n] = Un[K*i+ip][n] 
                                  + 135.0*dt[K*i+ip]/64.0*dUdt[K*i+ip][s-1][n] 
                                  - 243.0*dt[K*i+ip]/128.0*dUdt[K*i+ip][s-2][n]
                                  +  69.0*dt[K*i+ip]/128.0*dUdt[K*i+ip][s-3][n];
            else if (s==4) U[K*i+ip][n] = Un[K*i+ip][n] 
                                   + 16.0*dt[K*i+ip]/15.0*dUdt[K*i+ip][s-1][n] 
                                   - 27.0*dt[K*i+ip]/5.0*dUdt[K*i+ip][s-2][n]
                                   + 27.0*dt[K*i+ip]/4.0*dUdt[K*i+ip][s-3][n]
                                   - 17.0*dt[K*i+ip]/12.0*dUdt[K*i+ip][s-4][n];
            else if (s==5) U[K*i+ip][n] = Un[K*i+ip][n]
                                   +  5.0*dt[K*i+ip]/144.0*dUdt[K*i+ip][s-1][n] 
                                   +  4.0*dt[K*i+ip]/27.0*dUdt[K*i+ip][s-2][n]
                                   + 13.0*dt[K*i+ip]/16.0*dUdt[K*i+ip][s-3][n]
                                   -  5.0*dt[K*i+ip]/16.0*dUdt[K*i+ip][s-4][n]
                                   + 65.0*dt[K*i+ip]/432.0*dUdt[K*i+ip][s-5][n];*/
          /*// RKF5 - Formula II
          }else if (nstage==6){
            if (s==1)      U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 1.0*dt[K*i+ip]/4.0*dUdt[K*i+ip][s-1][n];
            else if (s==2) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 9.0*dt[K*i+ip]/32.0*dUdt[K*i+ip][s-1][n] 
                                    + 3.0*dt[K*i+ip]/32.0*dUdt[K*i+ip][s-2][n];
            else if (s==3) U[K*i+ip][n] = Un[K*i+ip][n]
                                 + 7296.0*dt[K*i+ip]/2197.0*dUdt[K*i+ip][s-1][n] 
                                 - 7200.0*dt[K*i+ip]/2197.0*dUdt[K*i+ip][s-2][n]
                                 + 1932.0*dt[K*i+ip]/2197.0*dUdt[K*i+ip][s-3][n];
            else if (s==4) U[K*i+ip][n] = Un[K*i+ip][n]
                                  - 845.0*dt[K*i+ip]/4104.0*dUdt[K*i+ip][s-1][n] 
                                  +3680.0*dt[K*i+ip]/513.0*dUdt[K*i+ip][s-2][n]
                                  -   8.0*dt[K*i+ip]*dUdt[K*i+ip][s-3][n]
                                  + 439.0*dt[K*i+ip]/216.0*dUdt[K*i+ip][s-4][n];
            else if (s==5) U[K*i+ip][n] = Un[K*i+ip][n] 
                                  -  11.0*dt[K*i+ip]/40.0*dUdt[K*i+ip][s-1][n] 
                                  +1859.0*dt[K*i+ip]/4104.0*dUdt[K*i+ip][s-2][n]
                                  -3544.0*dt[K*i+ip]/2565.0*dUdt[K*i+ip][s-3][n]
                                  +   2.0*dt[K*i+ip]*dUdt[K*i+ip][s-4][n]
                                  -   8.0*dt[K*i+ip]/27.0*dUdt[K*i+ip][s-5][n];*/
          /*// RKF5 - Formula IV
          }else if (nstage==6){
            if (s==1)      U[K*i+ip][n] = Un[K*i+ip][n]
                                    + 1.0*dt[K*i+ip]/2.0*dUdt[K*i+ip][s-1][n];
            else if (s==2) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 1.0*dt[K*i+ip]/4.0*dUdt[K*i+ip][s-1][n] 
                                    + 1.0*dt[K*i+ip]/4.0*dUdt[K*i+ip][s-2][n];
            else if (s==3) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 2.0*dt[K*i+ip]*dUdt[K*i+ip][s-1][n] 
                                    - 1.0*dt[K*i+ip]*dUdt[K*i+ip][s-2][n];
            else if (s==4) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 1.0*dt[K*i+ip]/27.0*dUdt[K*i+ip][s-1][n] 
                                   + 10.0*dt[K*i+ip]/27.0*dUdt[K*i+ip][s-3][n]
                                   +  7.0*dt[K*i+ip]/27.0*dUdt[K*i+ip][s-4][n];
            else if (s==5) U[K*i+ip][n] = Un[K*i+ip][n] 
                                  - 378.0*dt[K*i+ip]/625.0*dUdt[K*i+ip][s-1][n] 
                                  +  54.0*dt[K*i+ip]/625.0*dUdt[K*i+ip][s-2][n]
                                  + 546.0*dt[K*i+ip]/625.0*dUdt[K*i+ip][s-3][n]
                                  - 125.0*dt[K*i+ip]/625.0*dUdt[K*i+ip][s-4][n]
                                  +  28.0*dt[K*i+ip]/625.0*dUdt[K*i+ip][s-5][n];*/
          /*// RKF5 - Cash–Karp
          }else if (nstage==6){
            if (s==1)      U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 1.0*dt[K*i+ip]/5.0*dUdt[K*i+ip][s-1][n];
            else if (s==2) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 9.0*dt[K*i+ip]/40.0*dUdt[K*i+ip][s-1][n] 
                                    + 3.0*dt[K*i+ip]/40.0*dUdt[K*i+ip][s-2][n];
            else if (s==3) U[K*i+ip][n] = Un[K*i+ip][n] 
                                   + 12.0*dt[K*i+ip]/10.0*dUdt[K*i+ip][s-1][n] 
                                    - 9.0*dt[K*i+ip]/10.0*dUdt[K*i+ip][s-2][n]
                                    + 3.0*dt[K*i+ip]/10.0*dUdt[K*i+ip][s-3][n];
            else if (s==4) U[K*i+ip][n] = Un[K*i+ip][n] 
                                   + 70.0*dt[K*i+ip]/54.0*dUdt[K*i+ip][s-1][n] 
                                   -140.0*dt[K*i+ip]/54.0*dUdt[K*i+ip][s-2][n]
                                   +135.0*dt[K*i+ip]/54.0*dUdt[K*i+ip][s-3][n]
                                   - 11.0*dt[K*i+ip]/54.0*dUdt[K*i+ip][s-4][n];
            else if (s==5) U[K*i+ip][n] = Un[K*i+ip][n] 
                                  + 253.0*dt[K*i+ip]/4096.0*dUdt[K*i+ip][s-1][n] 
                                + 44275.0*dt[K*i+ip]/110592.0*dUdt[K*i+ip][s-2][n]
                                  + 575.0*dt[K*i+ip]/13824.0*dUdt[K*i+ip][s-3][n]
                                  + 175.0*dt[K*i+ip]/512.0*dUdt[K*i+ip][s-4][n]
                                 + 1631.0*dt[K*i+ip]/55296.0*dUdt[K*i+ip][s-5][n];*/
          // RK5
          }else if (nstage==6){
            if (s==1)      U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 1.0*dt[K*i+ip]/4.0*dUdt[K*i+ip][s-1][n];
            else if (s==2) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 1.0*dt[K*i+ip]/8.0*dUdt[K*i+ip][s-1][n] 
                                    + 1.0*dt[K*i+ip]/8.0*dUdt[K*i+ip][s-2][n];
            else if (s==3) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 2.0*dt[K*i+ip]/2.0*dUdt[K*i+ip][s-1][n] 
                                    - 1.0*dt[K*i+ip]/2.0*dUdt[K*i+ip][s-2][n];
            else if (s==4) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 9.0*dt[K*i+ip]/16.0*dUdt[K*i+ip][s-1][n] 
                                    + 3.0*dt[K*i+ip]/16.0*dUdt[K*i+ip][s-4][n];
            else if (s==5) U[K*i+ip][n] = Un[K*i+ip][n] 
                                    + 8.0*dt[K*i+ip]/7.0*dUdt[K*i+ip][s-1][n] 
                                    -12.0*dt[K*i+ip]/7.0*dUdt[K*i+ip][s-2][n]
                                    +12.0*dt[K*i+ip]/7.0*dUdt[K*i+ip][s-3][n]
                                    + 2.0*dt[K*i+ip]/7.0*dUdt[K*i+ip][s-4][n]
                                    - 3.0*dt[K*i+ip]/7.0*dUdt[K*i+ip][s-5][n];
          }else{
            throw invalid_argument("Stage number out of range"); 
          }
        }
      }
    }

    if (s>0){
      // Enforce BCs
      ApplyBCs();
 
      // Update primitive variables
      UpdatePrimitive();

      // Enforce positivity
      EnforcePositivity();
    }

    // Update gradient
    Gradient();

    // Limiter
    Limiter();

    // Reset residual
    ResetResidual(s);

    // Compute residual
    this->Residual(s);
  }

  // Update residual
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      for (int n=0; n<nvar; n++){
        switch (nstage){
          case 1: // RK1
            dUdt[K*i+ip][0][n] *= dt[K*i+ip];
            break;
          case 2: // RK2
            dUdt[K*i+ip][0][n] += dUdt[K*i+ip][1][n];
            dUdt[K*i+ip][0][n] *= dt[K*i+ip]/2.0;
            break;
          case 3: // RK3
            dUdt[K*i+ip][0][n] += 3.0*dUdt[K*i+ip][2][n];
            dUdt[K*i+ip][0][n] *= dt[K*i+ip]/4.0;
            break;
          case 4: // RK4
            dUdt[K*i+ip][0][n] += 2.0*dUdt[K*i+ip][1][n] 
                                + 2.0*dUdt[K*i+ip][2][n] 
                                    + dUdt[K*i+ip][3][n];
            dUdt[K*i+ip][0][n] *= dt[K*i+ip]/6.0;
            break;
          case 5: // RKM4
            dUdt[K*i+ip][0][n] += 4.0*dUdt[K*i+ip][3][n] 
                                    + dUdt[K*i+ip][4][n];
            dUdt[K*i+ip][0][n] *= dt[K*i+ip]/6.0;
            break;
          /*case 6: // RKF5 - Formula I
            dUdt[K*i+ip][0][n] = 47.0*dUdt[K*i+ip][0][n] 
                              + 216.0*dUdt[K*i+ip][2][n]
                               + 64.0*dUdt[K*i+ip][3][n] 
                               + 15.0*dUdt[K*i+ip][4][n] 
                              + 108.0*dUdt[K*i+ip][5][n];
            dUdt[K*i+ip][0][n] *= dt[K*i+ip]/450.0;
            break;*/
          /*case 6: // RKF5 - Formula II
            dUdt[K*i+ip][0][n] =    16.0/135.0*dUdt[K*i+ip][0][n]
                             +  6656.0/12825.0*dUdt[K*i+ip][2][n] 
                             + 28561.0/56430.0*dUdt[K*i+ip][3][n] 
                                    - 9.0/50.0*dUdt[K*i+ip][4][n] 
                                    + 2.0/55.0*dUdt[K*i+ip][5][n];
            dUdt[K*i+ip][0][n] *= dt[K*i+ip];
            break;*/
          /*case 6: // RKF5 - Formula IV
            dUdt[K*i+ip][0][n] =  1.0/24.0*dUdt[K*i+ip][0][n]
                                + 5.0/48.0*dUdt[K*i+ip][3][n] 
                               + 27.0/56.0*dUdt[K*i+ip][4][n] 
                             + 125.0/336.0*dUdt[K*i+ip][5][n];
            dUdt[K*i+ip][0][n] *= dt[K*i+ip];
            break;*/
          /*case 6: // RKF5 - Cash-Karp
            dUdt[K*i+ip][0][n] = 37.0/378.0*dUdt[K*i+ip][0][n]
                              + 250.0/621.0*dUdt[K*i+ip][2][n] 
                              + 125.0/594.0*dUdt[K*i+ip][3][n] 
                             + 512.0/1771.0*dUdt[K*i+ip][5][n];
            dUdt[K*i+ip][0][n] *= dt[K*i+ip];
            break;*/
          case 6: // RK5
            dUdt[K*i+ip][0][n] =   7.0/90.0*dUdt[K*i+ip][0][n]
                                + 32.0/90.0*dUdt[K*i+ip][2][n] 
                                + 12.0/90.0*dUdt[K*i+ip][3][n] 
                                + 32.0/90.0*dUdt[K*i+ip][4][n] 
                                 + 7.0/90.0*dUdt[K*i+ip][5][n];
            dUdt[K*i+ip][0][n] *= dt[K*i+ip];
            break;
          default:
            throw invalid_argument("Stage number out of range"); 
        }
      }
    }
  }

  // Update solution
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      for (int n=0; n<nvar; n++){
        U[K*i+ip][n] = Un[K*i+ip][n] + dUdt[K*i+ip][0][n];
      }
    }
  }
}


// **************************************************************************//
// RungeKutta (2D)
// **************************************************************************//
inline void Solver2D::RungeKutta(const int& nstage){

  // DESCRIPTION
  // ----------------------
  // Semi-discrete Runge-Kutta multi-stage time-marching approach for 2D 
  // problems.


  // INPUTS
  // ----------------------
  // nstage  - Number of stages (1-6) of the Runge-Kutta method


  // RUNGKE KUTTA
  // ----------------------
  // Enforce BCs
  ApplyBCs();

  // Update primitive variables
  UpdatePrimitive();

  // Enforce positivity
  EnforcePositivity();

  // Stage residual
  for (int s=0; s<nstage; s++){
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
        for (int ip=0; ip<K; ip++){
          for (int jp=0; jp<K; jp++){
            for (int n=0; n<nvar; n++){
              if (s==0) Un[K*i+ip][K*j+jp][n] = U[K*i+ip][K*j+jp][n];
              // RK1
              if (nstage==1){
              // RK2
              }else if (nstage==2){
                if (s==1)      U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                    + dt[K*i+ip][K*j+jp]*dUdt[K*i+ip][K*j+jp][s-1][n];
              // RK3
              }else if (nstage==3){
                if (s==1)      U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                    + dt[K*i+ip][K*j+jp]/3.0*dUdt[K*i+ip][K*j+jp][s-1][n];
                else if (s==2) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 2.0*dt[K*i+ip][K*j+jp]/3.0*dUdt[K*i+ip][K*j+jp][s-1][n];
              // RK4
              }else if (nstage==4){
                if (s==1)      U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] +
                                                      dt[K*i+ip][K*j+jp]/2.0*dUdt[K*i+ip][K*j+jp][s-1][n];
                else if (s==2) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                    + dt[K*i+ip][K*j+jp]/2.0*dUdt[K*i+ip][K*j+jp][s-1][n];
                else if (s==3) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                    + dt[K*i+ip][K*j+jp]*dUdt[K*i+ip][K*j+jp][s-1][n];
              // RKM4
              }else if (nstage==5){
                if (s==1)      U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/3.0*dUdt[K*i+ip][K*j+jp][s-1][n];
                else if (s==2) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/6.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/6.0*dUdt[K*i+ip][K*j+jp][s-2][n];
                else if (s==3) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 3.0*dt[K*i+ip][K*j+jp]/8.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/8.0*dUdt[K*i+ip][K*j+jp][s-3][n];
                else if (s==4) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 2.0*dt[K*i+ip][K*j+jp]/1.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                - 3.0*dt[K*i+ip][K*j+jp]/2.0*dUdt[K*i+ip][K*j+jp][s-2][n]
                                                + 1.0*dt[K*i+ip][K*j+jp]/2.0*dUdt[K*i+ip][K*j+jp][s-4][n];
              /*// RKF5 - Formula I
              }else if (nstage==6){
                if (s==1)      U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 2.0*dt[K*i+ip][K*j+jp]/9.0*dUdt[K*i+ip][K*j+jp][s-1][n];
                else if (s==2) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/4.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/12.0*dUdt[K*i+ip][K*j+jp][s-2][n];
                else if (s==3) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                              + 135.0*dt[K*i+ip][K*j+jp]/64.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                              - 243.0*dt[K*i+ip][K*j+jp]/128.0*dUdt[K*i+ip][K*j+jp][s-2][n]
                                               + 69.0*dt[K*i+ip][K*j+jp]128.0*dUdt[K*i+ip][K*j+jp][s-3][n];
                else if (s==4) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                               + 16.0*dt[K*i+ip][K*j+jp]/15.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                              -  27.0*dt[K*i+ip][K*j+jp]/5.0*dUdt[K*i+ip][K*j+jp][s-2][n]
                                              +  27.0*dt[K*i+ip][K*j+jp]/4.0*dUdt[K*i+ip][K*j+jp][s-3][n]
                                              -  17.0*dt[K*i+ip][K*j+jp]/12.0*dUdt[K*i+ip][K*j+jp][s-4][n];
                else if (s==5) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 5.0*dt[K*i+ip][K*j+jp]/144.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                              +   4.0*dt[K*i+ip][K*j+jp]/27.0*dUdt[K*i+ip][K*j+jp][s-2][n]
                                              +  13.0*dt[K*i+ip][K*j+jp]/16.0*dUdt[K*i+ip][K*j+jp][s-3][n]
                                              -   5.0*dt[K*i+ip][K*j+jp]/16.0*dUdt[K*i+ip][K*j+jp][s-4][n]
                                              +  65.0*dt[K*i+ip][K*j+jp]/432.0*dUdt[K*i+ip][K*j+jp][s-5][n];*/
              /*// RKF5 - Formula II
              }else if (nstage==6){
                if (s==1)      U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/4.0*dUdt[K*i+ip][K*j+jp][s-1][n];
                else if (s==2) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp]n]
                                                + 9.0*dt[K*i+ip][K*j+jp]/32.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                + 3.0*dt[K*i+ip][K*j+jp]/32.0*dUdt[K*i+ip][K*j+jp][s-2][n];
                else if (s==3) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                             + 7296.0*dt[K*i+ip][K*j+jp]/2197.0*dUdt[K*i+ip][K*j+jp]s-1][n] 
                                             - 7200.0*dt[K*i+ip][K*j+jp]/2197.0*dUdt[K*i+ip][K*j+jp][s-2][n]
                                             + 1932.0*dt[K*i+ip][K*j+jp]/2197.0*dUdt[K*i+ip][K*j+jp][s-3][n];
                else if (s==4) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                              - 845.0*dt[K*i+ip][K*j+jp]/4104.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                              +3680.0*dt[K*i+ip][K*j+jp]/513.0*dUdt[K*i+ip][K*j+jp][s-2][n]
                                              -   8.0*dt[K*i+ip][K*j+jp]*dUdt[K*i+ip][K*j+jp][s-3][n]
                                              + 439.0*dt[K*i+ip][K*j+jp]/216.0*dUdt[K*i+ip][K*j+jp][s-4][n];
                else if (s==5) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                              -  11.0*dt[K*i+ip][K*j+jp]/40.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                             + 1859.0*dt[K*i+ip][K*j+jp]/4104.0*dUdt[K*i+ip][K*j+jp][s-2][n]
                                              -3544.0*dt[K*i+ip][K*j+jp]/2565.0*dUdt[K*i+ip][K*j+jp][s-3][n]
                                              +   2.0*dt[K*i+ip][K*j+jp]*dUdt[K*i+ip][K*j+jp][s-4][n]
                                              -   8.0*dt[K*i+ip][K*j+jp]/27.0*dUdt[K*i+ip][K*j+jp][s-5][n];*/
              /*// RKF5 - Formula IV
              }else if (nstage==6){
                if (s==1)      U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/2.0*dUdt[K*i+ip][K*j+jp][s-1][n];
                else if (s==2) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/4.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/4.0*dUdt[K*i+ip][K*j+jp][s-2][n];
                else if (s==3) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 2.0*dt[K*i+ip][K*j+jp]*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                - 1.0*dt[K*i+ip][K*j+jp]*dUdt[K*i+ip][K*j+jp][s-2][n];
                else if (s==4) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/27.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                               + 10.0*dt[K*i+ip][K*j+jp]/27.0*dUdt[K*i+ip][K*j+jp][s-3][n]
                                               +  7.0*dt[K*i+ip][K*j+jp]/27.0*dUdt[K*i+ip][K*j+jp][s-4][n];
                else if (s==5) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                              - 378.0*dt[K*i+ip][K*j+jp]/625.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                              +  54.0*dt[K*i+ip][K*j+jp]/625.0*dUdt[K*i+ip][K*j+jp][s-2][n]
                                              + 546.0*dt[K*i+ip][K*j+jp]/625.0*dUdt[K*i+ip][K*j+jp][s-3][n]
                                              - 125.0*dt[K*i+ip][K*j+jp]/625.0*dUdt[K*i+ip][K*j+jp][s-4][n]
                                              +  28.0*dt[K*i+ip][K*j+jp]/625.0*dUdt[K*i+ip][K*j+jp][s-5][n];*/
              /*// RKF5 - Cash–Karp
              }else if (nstage==6){
                if (s==1)      U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/5.0*dUdt[K*i+ip][K*j+jp][s-1][n];
                else if (s==2) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 9.0*dt[K*i+ip][K*j+jp]/40.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                + 3.0*dt[K*i+ip][K*j+jp]/40.0*dUdt[K*i+ip][K*j+jp][s-2][n];
                else if (s==3) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                               + 12.0*dt[K*i+ip][K*j+jp]/10.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                - 9.0*dt[K*i+ip][K*j+jp]/10.0*dUdt[K*i+ip][K*j+jp][s-2][n]
                                                + 3.0*dt[K*i+ip][K*j+jp]/10.0*dUdt[K*i+ip][K*j+jp][s-3][n];
                else if (s==4) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                               + 70.0*dt[K*i+ip][K*j+jp]/54.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                               -140.0*dt[K*i+ip][K*j+jp]/54.0*dUdt[K*i+ip][K*j+jp][s-2][n]
                                               +135.0*dt[K*i+ip][K*j+jp]/54.0*dUdt[K*i+ip][K*j+jp][s-3][n]
                                               - 11.0*dt[K*i+ip][K*j+jp]/54.0*dUdt[K*i+ip][K*j+jp][s-4][n];
                else if (s==5) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n]
                                              + 253.0*dt[K*i+ip][K*j+jp]/4096.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                            + 44275.0*dt[K*i+ip][K*j+jp]/110592.0*dUdt[K*i+ip][K*j+jp][s-2][n]
                                              + 575.0*dt[K*i+ip][K*j+jp]/13824.0*dUdt[K*i+ip][K*j+jp][s-3][n]
                                              + 175.0*dt[K*i+ip][K*j+jp]/512.0*dUdt[K*i+ip][K*j+jp][s-4][n]
                                             + 1631.0*dt[K*i+ip][K*j+jp]/55296.0*dUdt[K*i+ip][K*j+jp][s-5][n];*/
              // RK5
              }else if (nstage==6){
                if (s==1)      U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/4.0*dUdt[K*i+ip][K*j+jp][s-1][n];
                else if (s==2) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/8.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                + 1.0*dt[K*i+ip][K*j+jp]/8.0*dUdt[K*i+ip][K*j+jp][s-2][n];
                else if (s==3) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 2.0*dt[K*i+ip][K*j+jp]/2.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                - 1.0*dt[K*i+ip][K*j+jp]/2.0*dUdt[K*i+ip][K*j+jp][s-2][n];
                else if (s==4) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 9.0*dt[K*i+ip][K*j+jp]/16.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                + 3.0*dt[K*i+ip][K*j+jp]/16.0*dUdt[K*i+ip][K*j+jp][s-4][n];
                else if (s==5) U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] 
                                                + 8.0*dt[K*i+ip][K*j+jp]/7.0*dUdt[K*i+ip][K*j+jp][s-1][n] 
                                                -12.0*dt[K*i+ip][K*j+jp]/7.0*dUdt[K*i+ip][K*j+jp][s-2][n]
                                                +12.0*dt[K*i+ip][K*j+jp]/7.0*dUdt[K*i+ip][K*j+jp][s-3][n]
                                                + 2.0*dt[K*i+ip][K*j+jp]/7.0*dUdt[K*i+ip][K*j+jp][s-4][n]
                                                - 3.0*dt[K*i+ip][K*j+jp]/7.0*dUdt[K*i+ip][K*j+jp][s-5][n];
              }else{
                throw invalid_argument("Stage number out of range"); 
              }
            }
          }
        }
      }
    }

    if (s>0){
      // Enforce BCs
      ApplyBCs();

      // Update primitive variables
      UpdatePrimitive();

      // Enforce positivity
      EnforcePositivity();
    }

    // Update gradient
    Gradient();

    // Limiter
    Limiter();

    // Reset residual
    ResetResidual(s);

    // Compute residual
    this->Residual(s);
  }

  // Update residual
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int ip=0; ip<K; ip++){
        for (int jp=0; jp<K; jp++){
          for (int n=0; n<nvar; n++){
            switch (nstage){
              case 1: //RK1
                dUdt[K*i+ip][K*j+jp][0][n] *= dt[K*i+ip][K*j+jp];
                break;
              case 2: //RK2
                dUdt[K*i+ip][K*j+jp][0][n] += dUdt[K*i+ip][K*j+jp][1][n];
                dUdt[K*i+ip][K*j+jp][0][n] *= dt[K*i+ip][K*j+jp]/2.0;
                break;
              case 3: //RK3
                dUdt[K*i+ip][K*j+jp][0][n] += 3.0*dUdt[K*i+ip][K*j+jp][2][n];
                dUdt[K*i+ip][K*j+jp][0][n] *= dt[K*i+ip][K*j+jp]/4.0;
                break;
              case 4: //RK4
                dUdt[K*i+ip][K*j+jp][0][n] += 2.0*dUdt[K*i+ip][K*j+jp][1][n] 
                                            + 2.0*dUdt[K*i+ip][K*j+jp][2][n] 
                                                + dUdt[K*i+ip][K*j+jp][3][n];
                dUdt[K*i+ip][K*j+jp][0][n] *= dt[K*i+ip][K*j+jp]/6.0;
                break;
              case 5: // RKM4
                dUdt[K*i+ip][K*j+jp][0][n] += 4.0*dUdt[K*i+ip][K*j+jp][3][n] 
                                                + dUdt[K*i+ip][K*j+jp][4][n];
                dUdt[K*i+ip][K*j+jp][0][n] *= dt[K*i+ip][K*j+jp]/6.0;
                break;
              /*case 6: // RKF5 - Formula I
                dUdt[K*i+ip][K*j+jp][0][n] = 47.0*dUdt[K*i+ip][K*j+jp][0][n] 
                                          + 216.0*dUdt[K*i+ip][K*j+jp][2][n] 
                                           + 64.0*dUdt[K*i+ip][K*j+jp][3][n] 
                                           + 15.0*dUdt[K*i+ip][K*j+jp][4][n] 
                                          + 108.0*dUdt[K*i+ip][K*j+jp][5][n];
                dUdt[K*i+ip][K*j+jp][0][n] *= dt[K*i+ip][K*j+jp]/450.0;
                break;*/
              /*case 6: // RKF5 - Formula II
                dUdt[K*i+ip][K*j+jp][0][n] =    16.0/135.0*dUdt[K*i+ip][K*j+jp][0][n]
                                         +  6656.0/12825.0*dUdt[K*i+ip][K*j+jp][2][n] 
                                         + 28561.0/56430.0*dUdt[K*i+ip][K*j+jp][3][n] 
                                                - 9.0/50.0*dUdt[K*i+ip][K*j+jp][4][n] 
                                                + 2.0/55.0*dUdt[K*i+ip][K*j+jp][5][n];
                dUdt[K*i+ip][K*j+jp][0][n] *= dt[K*i+ip][K*j+jp];
                break;*/
              /*case 6: // RKF5 - Formula IV
                dUdt[K*i+ip][K*j+jp][0][n] =  1.0/24.0*dUdt[K*i+ip][K*j+jp][0][n]
                                            + 5.0/48.0*dUdt[K*i+ip][K*j+jp][3][n] 
                                           + 27.0/56.0*dUdt[K*i+ip][K*j+jp][4][n] 
                                         + 125.0/336.0*dUdt[K*i+ip][K*j+jp][5][n];
                dUdt[K*i+ip][K*j+jp][0][n] *= dt[K*i+ip][K*j+jp];
                break;*/
              /*case 6: // RKF5 - Cash-Karp
                dUdt[K*i+ip][K*j+jp][0][n] = 37.0/378.0*dUdt[K*i+ip][K*j+jp][0][n]
                                          + 250.0/621.0*dUdt[K*i+ip][K*j+jp][2][n] 
                                          + 125.0/594.0*dUdt[K*i+ip][K*j+jp][3][n] 
                                         + 512.0/1771.0*dUdt[K*i+ip][K*j+jp][5][n];
                dUdt[K*i+ip][K*j+jp][0][n] *= dt[K*i+ip][K*j+jp];
                break;*/
              case 6: // RK5
                dUdt[K*i+ip][K*j+jp][0][n] =   7.0/90.0*dUdt[K*i+ip][K*j+jp][0][n]
                                            + 32.0/90.0*dUdt[K*i+ip][K*j+jp][2][n] 
                                            + 12.0/90.0*dUdt[K*i+ip][K*j+jp][3][n] 
                                            + 32.0/90.0*dUdt[K*i+ip][K*j+jp][4][n] 
                                             + 7.0/90.0*dUdt[K*i+ip][K*j+jp][5][n];
                dUdt[K*i+ip][K*j+jp][0][n] *= dt[K*i+ip][K*j+jp];
                break;
              default:
                throw invalid_argument("Stage number out of range"); 
            }
          }
        }
      }
    }
  }

  // Update solution
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int ip=0; ip<K; ip++){
        for (int jp=0; jp<K; jp++){
          for (int n=0; n<nvar; n++){
           U[K*i+ip][K*j+jp][n] = Un[K*i+ip][K*j+jp][n] + dUdt[K*i+ip][K*j+jp][0][n];
          }
        }
      }
    }
  }
}


// **************************************************************************//
// RungeKutta (3D)
// **************************************************************************//
inline void Solver3D::RungeKutta(const int& nstage){

  // DESCRIPTION
  // ----------------------
  // Semi-discrete Runge-Kutta multi-stage time-marching approach for 3D 
  // problems.


  // INPUTS
  // ----------------------
  // nstage  - Number of stages (1-6) of the Runge-Kutta method


  // RUNGKE KUTTA
  // ----------------------
  // Enforce BCs
  ApplyBCs();

  // Update primitive variables
  UpdatePrimitive();

  // Enforce positivity
  EnforcePositivity();

  // Stage residual
  for (int s=0; s<nstage; s++){
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
        for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){
              for (int kp=0; kp<K; kp++){
                for (int n=0; n<nvar; n++){
                  if (s==0) Un[K*i+ip][K*j+jp][K*k+kp][n] = U[K*i+ip][K*j+jp][K*k+kp][n];
                  // RK1
                  if (nstage==1){
                  // RK2
                  }else if (nstage==2){
                    if (s==1)      U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                                + dt[K*i+ip][K*j+jp][K*k+kp]*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n];
                  // RK3
                  }else if (nstage==3){
                    if (s==1)      U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n]
                                                                + dt[K*i+ip][K*j+jp][K*k+kp]/3.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n];
                    else if (s==2) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 2.0*dt[K*i+ip][K*j+jp][K*k+kp]/3.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n];
                  // RK4
                  }else if (nstage==4){
                    if (s==1)      U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                                + dt[K*i+ip][K*j+jp][K*k+kp]/2.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n];
                    else if (s==2) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                                + dt[K*i+ip][K*j+jp][K*k+kp]/2.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n];
                    else if (s==3) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                                + dt[K*i+ip][K*j+jp][K*k+kp]*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n];
                  // RKM4
                  }else if (nstage==5){
                    if (s==1)      U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/3.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n];
                    else if (s==2) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/6.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/6.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n];
                    else if (s==3) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 3.0*dt[K*i+ip][K*j+jp][K*k+kp]/8.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/8.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n];
                    else if (s==4) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 2.0*dt[K*i+ip][K*j+jp][K*k+kp]/1.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            - 3.0*dt[K*i+ip][K*j+jp][K*k+kp]/2.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n]
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/2.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-4][n];
                  /*// RKF5 - Formula I
                  }else if (nstage==6){
                    if (s==1)      U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 2.0*dt[K*i+ip][K*j+jp][K*k+kp]/9.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n];
                    else if (s==2) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/4.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/12.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n];
                    else if (s==3) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                          + 135.0*dt[K*i+ip][K*j+jp][K*k+kp]/64.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                          - 243.0*dt[K*i+ip][K*j+jp][K*k+kp]/128.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n]
                                                           + 69.0*dt[K*i+ip][K*j+jp][K*k+kp]/128.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n];
                    else if (s==4) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                           + 16.0*dt[K*i+ip][K*j+jp][K*k+kp]/15.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                          -  27.0*dt[K*i+ip][K*j+jp][K*k+kp]/5.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n]
                                                          +  27.0*dt[K*i+ip][K*j+jp][K*k+kp]/4.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n]
                                                          -  17.0*dt[K*i+ip][K*j+jp][K*k+kp]/12.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-4][n];
                    else if (s==5) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                          +   5.0*dt[K*i+ip][K*j+jp][K*k+kp]/144.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                          +   4.0*dt[K*i+ip][K*j+jp][K*k+kp]/27.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n]
                                                          +  13.0*dt[K*i+ip][K*j+jp][K*k+kp]/16.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n]
                                                          -   5.0*dt[K*i+ip][K*j+jp][K*k+kp]/16.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-4][n]
                                                          +  65.0*dt[K*i+ip][K*j+jp][K*k+kp]/432.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-5][n];*/
                  /*// RKF5 - Formula II
                  }else if (nstage==6){
                    if (s==1)      U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/4.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n];
                    else if (s==2) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 9.0*dt[K*i+ip][K*j+jp][K*k+kp]/32.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            + 3.0*dt[K*i+ip][K*j+jp][K*k+kp]/32.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n];
                    else if (s==3) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                         + 7296.0*dt[K*i+ip][K*j+jp][K*k+kp]/2197.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                         - 7200.0*dt[K*i+ip][K*j+jp][K*k+kp]/2197.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n]
                                                         + 1932.0*dt[K*i+ip][K*j+jp][K*k+kp]/2197.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n];
                    else if (s==4) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                          - 845.0*dt[K*i+ip][K*j+jp][K*k+kp]/4104.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                          +3680.0*dt[K*i+ip][K*j+jp][K*k+kp]/513.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n]
                                                          -   8.0*dt[K*i+ip][K*j+jp][K*k+kp]*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n]
                                                          + 439.0*dt[K*i+ip][K*j+jp][K*k+kp]/216.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-4][n];
                    else if (s==5) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                          -  11.0*dt[K*i+ip][K*j+jp][K*k+kp]/40.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                          +1859.0*dt[K*i+ip][K*j+jp][K*k+kp]/4104.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n]
                                                          -3544.0*dt[K*i+ip][K*j+jp][K*k+kp]/2565.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n]
                                                          +   2.0*dt[K*i+ip][K*j+jp][K*k+kp]*dUdt[K*i+ip][K*j+jp][K*k+kp][s-4][n]
                                                          -   8.0*dt[K*i+ip][K*j+jp][K*k+kp]/27.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-5][n];*/
                  /*// RKF5 - Formula IV
                  }else if (nstage==6){
                    if (s==1)      U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/2.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n];
                    else if (s==2) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/4.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/4.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n];
                    else if (s==3) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 2.0*dt[K*i+ip][K*j+jp][K*k+kp]*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            - 1.0*dt[K*i+ip][K*j+jp][K*k+kp]*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n];
                    else if (s==4) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] + 
                                                              1.0*dt[K*i+ip][K*j+jp][K*k+kp]/27.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                           + 10.0*dt[K*i+ip][K*j+jp][K*k+kp]/27.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n]
                                                           +  7.0*dt[K*i+ip][K*j+jp][K*k+kp]/27.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-4][n];
                    else if (s==5) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] - 
                                                            378.0*dt[K*i+ip][K*j+jp][K*k+kp]/625.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                          +  54.0*dt[K*i+ip][K*j+jp][K*k+kp]/625.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n]
                                                          + 546.0*dt[K*i+ip][K*j+jp][K*k+kp]/625.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n]
                                                          - 125.0*dt[K*i+ip][K*j+jp][K*k+kp]/625.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-4][n]
                                                          +  28.0*dt[K*i+ip][K*j+jp][K*k+kp]/625.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-5][n];*/
                  /*// RKF5 - Cash–Karp
                  }else if (nstage==6){
                    if (s==1)      U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/5.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n];
                    else if (s==2) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 9.0*dt[K*i+ip][K*j+jp][K*k+kp]/40.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            + 3.0*dt[K*i+ip][K*j+jp][K*k+kp]/40.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n];
                    else if (s==3) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                           + 12.0*dt[K*i+ip][K*j+jp][K*k+kp]/10.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            - 9.0*dt[K*i+ip][K*j+jp][K*k+kp]/10.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n]
                                                            + 3.0*dt[K*i+ip][K*j+jp][K*k+kp]/10.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n];
                    else if (s==4) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                           + 70.0*dt[K*i+ip][K*j+jp][K*k+kp]/54.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                           -140.0*dt[K*i+ip][K*j+jp][K*k+kp]/54.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n]
                                                           +135.0*dt[K*i+ip][K*j+jp][K*k+kp]/54.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n]
                                                           - 11.0*dt[K*i+ip][K*j+jp][K*k+kp]/54.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-4][n];
                    else if (s==5) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                          + 253.0*dt[K*i+ip][K*j+jp][K*k+kp]/4096.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                        + 44275.0*dt[K*i+ip][K*j+jp][K*k+kp]/110592.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n]
                                                          + 575.0*dt[K*i+ip][K*j+jp][K*k+kp]/13824.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n]
                                                          + 175.0*dt[K*i+ip][K*j+jp][K*k+kp]/512.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-4][n]
                                                         + 1631.0*dt[K*i+ip][K*j+jp][K*k+kp]/55296.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-5][n];*/
                  // RK5
                  }else if (nstage==6){
                    if (s==1)      U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/4.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n];
                    else if (s==2) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/8.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            + 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/8.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n];
                    else if (s==3) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 2.0*dt[K*i+ip][K*j+jp][K*k+kp]/2.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            - 1.0*dt[K*i+ip][K*j+jp][K*k+kp]/2.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n];
                    else if (s==4) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 9.0*dt[K*i+ip][K*j+jp][K*k+kp]/16.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            + 3.0*dt[K*i+ip][K*j+jp][K*k+kp]/16.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-4][n];
                    else if (s==5) U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                                            + 8.0*dt[K*i+ip][K*j+jp][K*k+kp]/7.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-1][n] 
                                                            -12.0*dt[K*i+ip][K*j+jp][K*k+kp]/7.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-2][n]
                                                            +12.0*dt[K*i+ip][K*j+jp][K*k+kp]/7.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-3][n]
                                                            + 2.0*dt[K*i+ip][K*j+jp][K*k+kp]/7.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-4][n]
                                                            - 3.0*dt[K*i+ip][K*j+jp][K*k+kp]/7.0*dUdt[K*i+ip][K*j+jp][K*k+kp][s-5][n];
                  }else{
                    throw invalid_argument("Stage number out of range"); 
                  }
                }
              }
            }
          }
        }
      }
    }

    if (s>0){
      // Enforce BCs
      ApplyBCs();

      // Update primitive variables
      UpdatePrimitive();

      // Enforce positivity
      EnforcePositivity();
    }

    // Update gradient
    Gradient();

    // Limiter
    Limiter();

    // Reset residual
    ResetResidual(s);

    // Compute residual
    this->Residual(s);
  }

  // Update residual
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
        for (int ip=0; ip<K; ip++){
          for (int jp=0; jp<K; jp++){
            for (int kp=0; kp<K; kp++){
              for (int n=0; n<nvar; n++){
                switch (nstage){
                  case 1: //RK1
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] *= dt[K*i+ip][K*j+jp][K*k+kp]; 
                    break;
                  case 2: //RK2
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] += dUdt[K*i+ip][K*j+jp][K*k+kp][1][n];
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] *= dt[K*i+ip][K*j+jp][K*k+kp]/2.0;
                    break;
                  case 3: //RK3
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] += 3.0*dUdt[K*i+ip][K*j+jp][K*k+kp][2][n];
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] *= dt[K*i+ip][K*j+jp][K*k+kp]/4.0;
                    break;
                  case 4: //RK4
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] += 2.0*dUdt[K*i+ip][K*j+jp][K*k+kp][1][n] 
                                                        + 2.0*dUdt[K*i+ip][K*j+jp][K*k+kp][2][n] 
                                                            + dUdt[K*i+ip][K*j+jp][K*k+kp][3][n];
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] *= dt[K*i+ip][K*j+jp][K*k+kp]/6.0;
                    break;
                  case 5: // RKM4
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] += 4.0*dUdt[K*i+ip][K*j+jp][K*k+kp][3][n] 
                                                            + dUdt[K*i+ip][K*j+jp][K*k+kp][4][n];
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] *= dt[K*i+ip][K*j+jp][K*k+kp]/6.0;
                    break;
                  /*case 6: // RKF5 - Formula I
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] = 47.0*dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] 
                                                      + 216.0*dUdt[K*i+ip][K*j+jp][K*k+kp][2][n] 
                                                       + 64.0*dUdt[K*i+ip][K*j+jp][K*k+kp][3][n] 
                                                       + 15.0*dUdt[K*i+ip][K*j+jp][K*k+kp][4][n] 
                                                      + 108.0*dUdt[K*i+ip][K*j+jp][K*k+kp][5][n];
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] *= dt[K*i+ip][K*j+jp][K*k+kp]/450.0;
                    break;*/
                  /*case 6: // RKF5 - Formula II
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] =    16.0/135.0*dUdt[K*i+ip][K*j+jp][K*k+kp][0][n]
                                                     +  6656.0/12825.0*dUdt[K*i+ip][K*j+jp][K*k+kp][2][n] 
                                                     + 28561.0/56430.0*dUdt[K*i+ip][K*j+jp][K*k+kp][3][n] 
                                                            - 9.0/50.0*dUdt[K*i+ip][K*j+jp][K*k+kp][4][n] 
                                                            + 2.0/55.0*dUdt[K*i+ip][K*j+jp][K*k+kp][5][n];
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] *= dt[K*i+ip][K*j+jp][K*k+kp];
                    break;*/
                  /*case 6: // RKF5 - Formula IV
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] =  1.0/24.0*dUdt[K*i+ip][K*j+jp][K*k+kp][0][n]
                                                        + 5.0/48.0*dUdt[K*i+ip][K*j+jp][K*k+kp][3][n] 
                                                       + 27.0/56.0*dUdt[K*i+ip][K*j+jp][K*k+kp][4][n] 
                                                     + 125.0/336.0*dUdt[K*i+ip][K*j+jp][K*k+kp][5][n];
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] *= dt[K*i+ip][K*j+jp][K*k+kp];
                    break;*/
                  /*case 6: // RKF5 - Cash-Karp
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] = 37.0/378.0*dUdt[K*i+ip][K*j+jp][K*k+kp][0][n]
                                                      + 250.0/621.0*dUdt[K*i+ip][K*j+jp][K*k+kp][2][n] 
                                                      + 125.0/594.0*dUdt[K*i+ip][K*j+jp][K*k+kp][3][n] 
                                                     + 512.0/1771.0*dUdt[K*i+ip][K*j+jp][K*k+kp][5][n];
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] *= dt[K*i+ip][K*j+jp][K*k+kp];
                    break;*/
                  case 6: // RK5
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] =   7.0/90.0*dUdt[K*i+ip][K*j+jp][K*k+kp][0][n]
                                                        + 32.0/90.0*dUdt[K*i+ip][K*j+jp][K*k+kp][2][n] 
                                                        + 12.0/90.0*dUdt[K*i+ip][K*j+jp][K*k+kp][3][n] 
                                                        + 32.0/90.0*dUdt[K*i+ip][K*j+jp][K*k+kp][4][n] 
                                                         + 7.0/90.0*dUdt[K*i+ip][K*j+jp][K*k+kp][5][n];
                    dUdt[K*i+ip][K*j+jp][K*k+kp][0][n] *= dt[K*i+ip][K*j+jp][K*k+kp];
                    break;
                  default:
                    throw invalid_argument("Stage number out of range"); 
                }
              }
            }
          }
        }
      }
    }
  }

  // Update solution
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
      for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
        for (int ip=0; ip<K; ip++){
          for (int jp=0; jp<K; jp++){
            for (int kp=0; kp<K; kp++){
              for (int n=0; n<nvar; n++){
                U[K*i+ip][K*j+jp][K*k+kp][n] = Un[K*i+ip][K*j+jp][K*k+kp][n] 
                                           + dUdt[K*i+ip][K*j+jp][K*k+kp][0][n];
              }
            }
          }
        }
      }
    }
  }
}
#endif
