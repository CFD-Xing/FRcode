///////////////////////////////////////////////////////////////////////////////
/// \file EulerSolver2D.cc
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 27 2021
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
// EULER 2D
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// EulerInitialCondition
// **************************************************************************//
void Solver2D::EulerInitialCondition(const int& option, double& t){

  // INPUTS
  // ----------------------
  // option  - Case option


  // OUTPUTS
  // ----------------------
  // t       - Run time


  // VARIABLE DECLARATION
  // ----------------------
  double rhol, ul, vl, pl;
  double rhor, ur, vr, pr;
  double xm, ym;


  // INITIAL CONDITION
  // ----------------------
  // Set initial condition
  switch (option){
    case 1: // Shockbox
      rhol = 1*1.225; ul = 0.0; vl = 0.0; pl = 1*101325.0;
      rhor = 4*1.225; ur = 0.0; vr = 0.0; pr = 4*101325.0;
      xm = 5.0; ym = 5.0; t = 0.007;
      BCtype = 3;
      break;
    case 2: // Oblique shock
      rhol = 1.225; ul = 700.0*cos(M_PI/18.0); vl = -700.0*sin(M_PI/18.0); pl = 101325.0;
      rhor = 1.225; ur = 700.0*cos(M_PI/18.0); vr = -700.0*sin(M_PI/18.0); pr = 101325.0;
      xm = 2.0; ym = 2.0; t = 0.050;
      BCtype = 4;
      break;
    default:
      throw invalid_argument("Unknown initial condition"); 
  }
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int ip=0; ip<K; ip++){
        for (int jp=0; jp<K; jp++){
          if (mesh.xc[i][j]+0.5*mesh.dx[i][j]*FR.eta[ip] < xm && 
              mesh.yc[i][j]+0.5*mesh.dy[i][j]*FR.eta[jp] < ym){
            W[K*i+ip][K*j+jp][0]=rhol; 
            W[K*i+ip][K*j+jp][1]=ul; 
            W[K*i+ip][K*j+jp][2]=vl; 
            W[K*i+ip][K*j+jp][3]=pl;
            Euler2D.PrimtoCons(&W[K*i+ip][K*j+jp][0], &U[K*i+ip][K*j+jp][0]);
          }else{
            W[K*i+ip][K*j+jp][0]=rhor; 
            W[K*i+ip][K*j+jp][1]=ur; 
            W[K*i+ip][K*j+jp][2]=vr; 
            W[K*i+ip][K*j+jp][3]=pr;
            Euler2D.PrimtoCons(&W[K*i+ip][K*j+jp][0], &U[K*i+ip][K*j+jp][0]);
          }
        }
      }
    }
  }
  // Save initial solution
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int ip=0; ip<K; ip++){
        for (int jp=0; jp<K; jp++){
          for (int n=0; n<nvar; n++){
            U0[K*i+ip][K*j+jp][n] = U[K*i+ip][K*j+jp][n];
            W0[K*i+ip][K*j+jp][n] = W[K*i+ip][K*j+jp][n];
          }
        }
      }
    }
  }
}


// **************************************************************************//
// EulerExactSolution
// **************************************************************************//
void Solver2D::EulerExactSolution(const int& option){

  // INPUTS
  // ----------------------
  // option  - Case option


  // VARIABLE DECLARATION
  // ----------------------
  double rhol, ul, vl, pl;
  double rhor, ur, vr, pr;
  double xm, ym;


  // EXACT SOLUTION
  // ----------------------

}


// **************************************************************************//
// EulerOutput
// **************************************************************************//
void Solver2D::EulerOutput(const int& option, const double& CFL, const int& nStage){

  // INPUTS
  // ----------------------
  // option  - Case option


  // VARIABLE DECLARATION
  // ----------------------
  ofstream outfile;
  char filename[256];


  // OUTPUT SOLUTION
  // ----------------------
  // Exact solution
  EulerExactSolution(option);

  // Numerical solution
  if (SolverType==1){ // FVM
    sprintf(filename, "Euler2D-Case%02d-FVM-%s-Grid=%dx%d-CFL=%0.2f-Stage=%d.dat",option,LimiterName,nx,ny,CFL,nStage);
    outfile.open (filename);
  }else{ // FR
    sprintf(filename, "Euler2D-Case%02d-FR%d-%s-Grid=%dx%d-CFL=%0.2f-Stage=%d.dat",option,K,FR.FRname,nx,ny,CFL,nStage);
    outfile.open (filename);
  }
  outfile << "VARIABLES = \"x\"\n";
  outfile << "\"y\"\n";
  outfile << "\"rho\"\n";
  outfile << "\"u\"\n";
  outfile << "\"v\"\n";
  outfile << "\"p\"\n";
  outfile << "ZONE T = \"Solution Points\"\n";
  outfile << "STRANDID=0, SOLUTIONTIME=0\n";
  outfile << "I=" << mesh.nx*K << ", J=" << mesh.ny*K << ", K=1, ZONETYPE=Ordered\n";
  outfile << "DATAPACKING=POINT\n";
  outfile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n";
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
        for (int jp=0; jp<K; jp++){
          outfile.precision(16);
          outfile << fixed;
          outfile << mesh.xc[i][j] + 0.5*mesh.dx[i][j]*FR.eta[ip] << " " 
                  << mesh.yc[i][j] + 0.5*mesh.dy[i][j]*FR.eta[jp] << " " 
                  << W[K*i+ip][K*j+jp][0] << " "
                  << W[K*i+ip][K*j+jp][1] << " "
                  << W[K*i+ip][K*j+jp][2] << " "
                  << W[K*i+ip][K*j+jp][3] << endl;
        }
      }
    }
    outfile << endl;
  }
  outfile << "VARIABLES = \"x\"\n";
  outfile << "\"y\"\n";
  outfile << "\"rho\"\n";
  outfile << "\"u\"\n";
  outfile << "\"v\"\n";
  outfile << "\"p\"\n";
  outfile << "ZONE T = \"Exact Solution\"\n";
  outfile << "STRANDID=0, SOLUTIONTIME=0\n";
  outfile << "I=" << mesh.nx*K << ", J=" << mesh.ny*K << ", K=1, ZONETYPE=Ordered\n";
  outfile << "DATAPACKING=POINT\n";
  outfile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n";
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
        for (int jp=0; jp<K; jp++){
          outfile.precision(16);
          outfile << fixed;
          outfile << mesh.xc[i][j] + 0.5*mesh.dx[i][j]*FR.eta[ip] << " " 
                  << mesh.yc[i][j] + 0.5*mesh.dy[i][j]*FR.eta[jp] << " " 
                  << W0[K*i+ip][K*j+jp][0] << " "
                  << W0[K*i+ip][K*j+jp][1] << " "
                  << W0[K*i+ip][K*j+jp][2] << " "
                  << W0[K*i+ip][K*j+jp][3] << endl;
        }
      }
    }
    outfile << endl;
  }
  outfile.close();
}
