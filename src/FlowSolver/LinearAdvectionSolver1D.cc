///////////////////////////////////////////////////////////////////////////////
/// \file LinearAdvectionSolver1D.cc
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 23 2021
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
// LINEAR ADVECTION 1D
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// LinearAdvectionInitialCondition
// **************************************************************************//
void Solver1D::LinearAdvectionInitialCondition(const int& option, const double& ncycle, double &t){

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
        for (int ip=0; ip<K; ip++){
          U[K*i+ip][0] = exp(-0.4*pow(mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] - L/2.0, 2));
        }
      }
      break;
    case 2: // Non-smooth case
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] > 10.0/3.0 && 
              mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] < 20.0/3.0){
            U[K*i+ip][0] = 1.0;
          } else{
            U[K*i+ip][0] = 0.0;
          }
        }
      }
      break;
    default:
      throw invalid_argument("Unknown initial condition"); 
  }
  // Save initial solution
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      W[K*i+ip][0]  = U[K*i+ip][0];
      U0[K*i+ip][0] = U[K*i+ip][0];
      W0[K*i+ip][0] = W[K*i+ip][0];
    }
  }
}


// **************************************************************************//
// LinearAdvectionExactSolution
// **************************************************************************//
void Solver1D::LinearAdvectionExactSolution(const int& option, const double& t){

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
        for (int ip=0; ip<K; ip++){
          double x = mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] - wavespeed*t;
          Un[K*i+ip][0] = exp(-0.4*pow(x-floor(x/L)*L-L/2.0, 2));
        }
      }
      break;
    case 2: // Non-smooth case
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          double x = mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] - wavespeed*t;
          if (x-floor(x/L)*L > 10.0/3.0 && x-floor(x/L)*L < 20.0/3.0){
            Un[K*i+ip][0] = 1.0;
          } else{
            Un[K*i+ip][0] = 0.0;
          }
        }
      }
      break;
    default:
      throw invalid_argument("Unknown exact solution"); 
  }
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      Wn[K*i+ip][0] = Un[K*i+ip][0];
    }
  }
}


// **************************************************************************//
// LinearAdvectionOutput
// **************************************************************************//
void Solver1D::LinearAdvectionOutput(const int& option, const double& CFL, const int& nStage, const double& t){

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
  double* f = new double[K];


  // OUTPUT SOLUTION
  // ----------------------
  // Exact solution
  LinearAdvectionExactSolution(option, t);

  // Numerical solution
  if (SolverType==1){ // FVM
    sprintf(filename, "LinearAdvection1D-Case%02d-FVM-%s-Grid=%d-CFL=%0.2f-Stage=%d.dat",option,LimiterName,nx,CFL,nStage);
    outfile.open (filename);
    outfile << "VARIABLES = \"x\"\n";
    outfile << "\"u\"\n";
    outfile << "ZONE T = \"FVM Solution\"\n";
    outfile << "STRANDID=0, SOLUTIONTIME=0\n";
    outfile << "I=" << mesh.nx << ", J=1, K=1, ZONETYPE=Ordered\n";
    outfile << "DATAPACKING=POINT\n";
    outfile << "DT=(DOUBLE DOUBLE)\n";
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      outfile.precision(16);
      outfile << fixed;
      outfile << mesh.xc[i] << " " 
              << W[i][0] << endl;
    }
    outfile << "VARIABLES = \"x\"\n";
    outfile << "\"u\"\n";
    outfile << "ZONE T = \"Exact Solution\"\n";
    outfile << "STRANDID=0, SOLUTIONTIME=0\n";
    outfile << "I=" << mesh.nx << ", J=1, K=1, ZONETYPE=Ordered\n";
    outfile << "DATAPACKING=POINT\n";
    outfile << "DT=(DOUBLE DOUBLE)\n";
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      outfile.precision(16);
      outfile << fixed;
      outfile << mesh.xc[i] << " " 
              << Wn[i][0] << endl;
    }
  }else{ // FR
    sprintf(filename, "LinearAdvection1D-Case%02d-FR%d-%s-Grid=%d-CFL=%0.2f-Stage=%d.dat",option,K,FR.FRname,nx,CFL,nStage);
    outfile.open (filename);
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      outfile << "VARIABLES = \"x\"\n";
      outfile << "\"u\"\n";
      outfile << "ZONE T = \"" << i <<"\"\n";
      outfile << "STRANDID=0, SOLUTIONTIME=0\n";
      outfile << "I=" << 11 << ", J=1, K=1, ZONETYPE=Ordered\n";
      outfile << "DATAPACKING=POINT\n";
      outfile << "DT=(DOUBLE DOUBLE)\n";
      outfile.precision(16);
      outfile << fixed;
      for (int ip=0; ip<K; ip++){
        f[ip] = W[K*i+ip][0];
      }
      for (int ip=-5; ip<=5; ip++){
        outfile << mesh.xc[i] + 0.1*ip*mesh.dx[i] << " ";
        outfile << LagrangeInterpol1D(K, FR.eta, &f[0], 0.2*ip) << " ";
        outfile << endl;
      }
    }
    outfile << endl;
    outfile << "VARIABLES = \"x\"\n";
    outfile << "\"u\"\n";
    outfile << "ZONE T = \"Solution Points\"\n";
    outfile << "STRANDID=0, SOLUTIONTIME=0\n";
    outfile << "I=" << mesh.nx*K << ", J=1, K=1, ZONETYPE=Ordered\n";
    outfile << "DATAPACKING=POINT\n";
    outfile << "DT=(DOUBLE DOUBLE)\n";
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      for (int ip=0; ip<K; ip++){
        outfile.precision(16);
        outfile << fixed;
        outfile << mesh.xc[i] + 0.5*mesh.dx[i]*FR.eta[ip] << " " 
                << W[K*i+ip][0] << " " << endl;
      }
    }
    outfile << endl;
    outfile << "VARIABLES = \"x\"\n";
    outfile << "\"u\"\n";
    outfile << "ZONE T = \"Exact Solution\"\n";
    outfile << "STRANDID=0, SOLUTIONTIME=0\n";
    outfile << "I=" << mesh.nx*K << ", J=1, K=1, ZONETYPE=Ordered\n";
    outfile << "DATAPACKING=POINT\n";
    outfile << "DT=(DOUBLE DOUBLE)\n";
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      for (int ip=0; ip<K; ip++){
        outfile.precision(16);
        outfile << fixed;
        outfile << mesh.xc[i] + 0.5*mesh.dx[i]*FR.eta[ip] << " " 
                << Wn[K*i+ip][0] << " " << endl;
      }
    }
  }
  outfile.close();

  // Deallocate
  if (f!=NULL) delete[] f;
}
