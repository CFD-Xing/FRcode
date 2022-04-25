///////////////////////////////////////////////////////////////////////////////
/// \file LinearAdvectionSolver2D.cc
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 26 2021
///
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <iostream>
#include <fstream>
using namespace std;

// Include
#include "FlowSolver.h"


///////////////////////////////////////////////////////////////////////////////
// LINEAR ADVECTION 2D
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// LinearAdvectionInitialCondition
// **************************************************************************//
void Solver2D::LinearAdvectionInitialCondition(const int& option, const double& ncycle, double& t){

  // INPUTS
  // ----------------------
  // option  - Case option
  // nCycle  - Number of solution cycles


  // OUTPUTS
  // ----------------------
  // t       - Run time


  // VARIABLE DECLARATION
  // ----------------------
  static double L = 10.0;


  // INITIAL CONDITION
  // ----------------------
  // Set initial condition
  wavespeed = 10.0;
  t = ncycle*L/fabs(wavespeed);
  BCtype = 2;
  switch (option){
    case 1: // Smooth case
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int j=0; j<mesh.ny+2*mesh.ng; j++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){
              U[K*i+ip][K*j+jp][0] = exp(-0.4*pow(mesh.xc[i][j]+0.5*mesh.dx[i][j]*FR.eta[ip] - L/2.0, 2))*
                                     exp(-0.4*pow(mesh.yc[i][j]+0.5*mesh.dy[i][j]*FR.eta[jp] - L/2.0, 2));
            }
          }
        }
      }
      break;
    case 2: // Non-smooth case
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int j=0; j<mesh.ny+2*mesh.ng; j++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){
              if ((mesh.xc[i][j]+0.5*mesh.dx[i][j]*FR.eta[ip] > 10.0/3.0 && 
                   mesh.xc[i][j]+0.5*mesh.dx[i][j]*FR.eta[ip] < 20.0/3.0) && 
                  (mesh.yc[i][j]+0.5*mesh.dy[i][j]*FR.eta[jp] > 10.0/3.0 && 
                   mesh.yc[i][j]+0.5*mesh.dy[i][j]*FR.eta[jp] < 20.0/3.0)){
                U[K*i+ip][K*j+jp][0] = 1.0;
              } else{
                U[K*i+ip][K*j+jp][0] = 0.0;
              }
            }
          }
        }
      }
      break;
    default:
      throw invalid_argument("Unknown initial condition"); 
  }
  // Save initial solution
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int ip=0; ip<K; ip++){
        for (int jp=0; jp<K; jp++){ 
          W[K*i+ip][K*j+jp][0]  = U[K*i+ip][K*j+jp][0];
          U0[K*i+ip][K*j+jp][0] = U[K*i+ip][K*j+jp][0];
          W0[K*i+ip][K*j+jp][0] = W[K*i+ip][K*j+jp][0];
        }
      }
    }
  }
}


// **************************************************************************//
// LinearAdvectionExactSolution
// **************************************************************************//
void Solver2D::LinearAdvectionExactSolution(const int& option, const double& t){

  // INPUTS
  // ----------------------
  // option  - Case option
  // t       - Run time


  // VARIABLE DECLARATION
  // ----------------------
  static double L = 10.0;


  // EXACT SOLUTION
  // ----------------------
  switch (option){
    case 1: // Smooth case
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int j=0; j<mesh.ny+2*mesh.ng; j++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){ 
              double x = mesh.xc[i][j]+0.5*mesh.dx[i][j]*FR.eta[ip] - wavespeed*t;
              double y = mesh.yc[i][j]+0.5*mesh.dy[i][j]*FR.eta[jp] - wavespeed*t;
              Un[K*i+ip][K*j+jp][0] = exp(-0.4*pow(x-floor(x/L)*L-L/2.0, 2))*
                                      exp(-0.4*pow(y-floor(y/L)*L-L/2.0, 2));
            }
          }
        }
      }
      break;
    case 2: // Non-smooth case
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int j=0; j<mesh.ny+2*mesh.ng; j++){
          for (int ip=0; ip<K; ip++){
            for (int jp=0; jp<K; jp++){ 
              double x = mesh.xc[i][j]+0.5*mesh.dx[i][j]*FR.eta[ip] - wavespeed*t;
              double y = mesh.yc[i][j]+0.5*mesh.dy[i][j]*FR.eta[jp] - wavespeed*t;
              if ((x-floor(x/L)*L > 10.0/3.0 && x-floor(x/L)*L < 20.0/3.0) &&
                  (y-floor(y/L)*L > 10.0/3.0 && y-floor(y/L)*L < 20.0/3.0)){
                Un[K*i+ip][K*j+jp][0] = 1.0;
              } else{
                Un[K*i+ip][K*j+jp][0] = 0.0;
              }
            }
          }
        }
      }
      break;
    default:
      throw invalid_argument("Unknown exact solution"); 
  }
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int ip=0; ip<K; ip++){
        for (int jp=0; jp<K; jp++){ 
          Wn[K*i+ip][K*j+jp][0] = Un[K*i+ip][K*j+jp][0];
        }
      }
    }
  }
}


// **************************************************************************//
// LinearAdvectionOutput
// **************************************************************************//
void Solver2D::LinearAdvectionOutput(const int& option, const double& CFL, const int& nStage, const double& t){

  // INPUTS
  // ----------------------
  // option  - Case option
  // CFL     - CFL number
  // nStage  - Number of stage
  // t       - Run time


  // VARIABLE DECLARATION
  // ----------------------
  ofstream outfile;
  char filename[256];


  // OUTPUT SOLUTION
  // ----------------------
  // Exact solution
  LinearAdvectionExactSolution(option, t);

  // Numerical solution
  if (SolverType==1){ // FVM
    sprintf(filename, "LinearAdvection2D-Case%02d-FVM-%s-Grid=%dx%d-CFL=%0.2f-Stage=%d.dat",option,LimiterName,nx,ny,CFL,nStage);
    outfile.open (filename);
  }else{ // FR
    sprintf(filename, "LinearAdvection2D-Case%02d-FR%d-%s-Grid=%dx%d-CFL=%0.2f-Stage=%d.dat",option,K,FR.FRname,nx,ny,CFL,nStage);
    outfile.open (filename);
  }
  outfile << "VARIABLES = \"x\"\n";
  outfile << "\"y\"\n";
  outfile << "\"u\"\n";
  outfile << "ZONE T = \"Solution Points\"\n";
  outfile << "STRANDID=0, SOLUTIONTIME=0\n";
  outfile << "I=" << mesh.nx*K << ", J=" << mesh.ny*K << ", K=1, ZONETYPE=Ordered\n";
  outfile << "DATAPACKING=POINT\n";
  outfile << "DT=(DOUBLE DOUBLE DOUBLE)\n";
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
        for (int jp=0; jp<K; jp++){
          outfile.precision(16);
          outfile << fixed;
          outfile << mesh.xc[i][j] + 0.5*mesh.dx[i][j]*FR.eta[ip] << " " 
                  << mesh.yc[i][j] + 0.5*mesh.dy[i][j]*FR.eta[jp] << " " 
                  << W[K*i+ip][K*j+jp][0] << endl;
        }
      }
    }
    outfile << endl;
  }
  outfile << "VARIABLES = \"x\"\n";
  outfile << "\"y\"\n";
  outfile << "\"u\"\n";
  outfile << "ZONE T = \"Exact Solution\"\n";
  outfile << "STRANDID=0, SOLUTIONTIME=0\n";
  outfile << "I=" << mesh.nx*K << ", J=" << mesh.ny*K << ", K=1, ZONETYPE=Ordered\n";
  outfile << "DATAPACKING=POINT\n";
  outfile << "DT=(DOUBLE DOUBLE DOUBLE)\n";
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
        for (int jp=0; jp<K; jp++){
          outfile.precision(16);
          outfile << fixed;
          outfile << mesh.xc[i][j] + 0.5*mesh.dx[i][j]*FR.eta[ip] << " " 
                  << mesh.yc[i][j] + 0.5*mesh.dy[i][j]*FR.eta[jp] << " " 
                  << Wn[K*i+ip][K*j+jp][0] << endl;
        }
      }
    }
    outfile << endl;
  }
  outfile.close();
}
