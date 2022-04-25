///////////////////////////////////////////////////////////////////////////////
/// \file EulerSolver3D.cc
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
// EULER 3D
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// EulerInitialCondition
// **************************************************************************//
void Solver3D::EulerInitialCondition(const int& option, double& t){

  // INPUTS
  // ----------------------
  // option  - Case option


  // OUTPUTS
  // ----------------------
  // t       - Run time


  // VARIABLE DECLARATION
  // ----------------------
  double rhol, ul, vl, wl, pl;
  double rhor, ur, vr, wr, pr;
  double xm, ym, zm;


  // INITIAL CONDITION
  // ----------------------
  // Set initial condition
  switch (option){
    case 1: // Shockbox
      rhol = 1*1.225; ul = 0.0; vl = 0.0; wl = 0.0; pl = 1*101325.0;
      rhor = 4*1.225; ur = 0.0; vr = 0.0; wr = 0.0; pr = 4*101325.0;
      xm = 5.0; ym = 5.0; zm = 5.0; t = 0.007;
      BCtype = 3;
      break;
    case 2: // Oblique shock
      rhol = 1.225; ul = 700.0*cos(M_PI/18.0); vl = 0.0; wl = -700.0*sin(M_PI/18.0); pl = 101325.0;
      rhor = 1.225; ur = 700.0*cos(M_PI/18.0); vr = 0.0; wr = -700.0*sin(M_PI/18.0); pr = 101325.0;
      xm = 2.0; ym = 2.0; zm = 2.0; t = 0.050;
      BCtype = 4;
      break;
    default:
      throw invalid_argument("Unknown initial condition"); 
  }
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int ip=0; ip<K; ip++){
          for (int jp=0; jp<K; jp++){
            for (int kp=0; kp<K; kp++){
              if (mesh.xc[i][j][k]+0.5*mesh.dx[i][j][k]*FR.eta[ip] < xm && 
                  mesh.yc[i][j][k]+0.5*mesh.dy[i][j][k]*FR.eta[jp] < ym && 
                  mesh.zc[i][j][k]+0.5*mesh.dz[i][j][k]*FR.eta[kp] < zm){
                W[K*i+ip][K*j+jp][K*k+kp][0]=rhol; 
                W[K*i+ip][K*j+jp][K*k+kp][1]=ul; 
                W[K*i+ip][K*j+jp][K*k+kp][2]=vl; 
                W[K*i+ip][K*j+jp][K*k+kp][3]=wl; 
                W[K*i+ip][K*j+jp][K*k+kp][4]=pl;
                Euler3D.PrimtoCons(&W[K*i+ip][K*j+jp][K*k+kp][0], &U[K*i+ip][K*j+jp][K*k+kp][0]);
              }else{
                W[K*i+ip][K*j+jp][K*k+kp][0]=rhor; 
                W[K*i+ip][K*j+jp][K*k+kp][1]=ur; 
                W[K*i+ip][K*j+jp][K*k+kp][2]=vr; 
                W[K*i+ip][K*j+jp][K*k+kp][3]=wr; 
                W[K*i+ip][K*j+jp][K*k+kp][4]=pr;
                Euler3D.PrimtoCons(&W[K*i+ip][K*j+jp][K*k+kp][0], &U[K*i+ip][K*j+jp][K*k+kp][0]);
              }
            }
          }
        }
      }
    }
  }
  // Save initial solution
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int j=0; j<mesh.ny+2*mesh.ng; j++){
      for (int k=0; k<mesh.nz+2*mesh.ng; k++){
        for (int ip=0; ip<K; ip++){
          for (int jp=0; jp<K; jp++){
            for (int kp=0; kp<K; kp++){
              for (int n=0; n<nvar; n++){
                U0[K*i+ip][K*j+jp][K*k+kp][n] = U[K*i+ip][K*j+jp][K*k+kp][n];
                W0[K*i+ip][K*j+jp][K*k+kp][n] = W[K*i+ip][K*j+jp][K*k+kp][n];
              }
            }
          }
        }
      }
    }
  }
}


// **************************************************************************//
// EulerExactSolution
// **************************************************************************//
void Solver3D::EulerExactSolution(const int& option){

  // INPUTS
  // ----------------------
  // option  - Case option


  // VARIABLE DECLARATION
  // ----------------------
  double rhol, ul, vl, wl, pl;
  double rhor, ur, vr, wr, pr;
  double xm, ym, zm;


  // EXACT SOLUTION
  // ----------------------

}


// **************************************************************************//
// EulerOutput
// **************************************************************************//
void Solver3D::EulerOutput(const int& option, const double& CFL, const int& nStage){

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
    sprintf(filename, "Euler3D-Case%02d-FVM-%s-Grid=%dx%dx%d-CFL=%0.2f-Stage=%d.dat",option,LimiterName,nx,ny,nz,CFL,nStage);
    outfile.open (filename);
  }else{ // FR
    sprintf(filename, "Euler3D-Case%02d-FR%d-%s-Grid=%dx%dx%d-CFL=%0.2f-Stage=%d.dat",option,K,FR.FRname,nx,ny,nz,CFL,nStage);
    outfile.open (filename);
  }
  outfile << "VARIABLES = \"x\"\n";
  outfile << "\"y\"\n";
  outfile << "\"z\"\n";
  outfile << "\"rho\"\n";
  outfile << "\"u\"\n";
  outfile << "\"v\"\n";
  outfile << "\"w\"\n";
  outfile << "\"p\"\n";
  outfile << "ZONE T = \"Solution Points\"\n";
  outfile << "STRANDID=0, SOLUTIONTIME=0\n";
  outfile << "I=" << mesh.nx*K << ", J=" << mesh.ny*K << ", K=" << mesh.nz*K << ", ZONETYPE=Ordered\n";
  outfile << "DATAPACKING=POINT\n";
  outfile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n";
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
        for (int jp=0; jp<K; jp++){
          for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
            for (int kp=0; kp<K; kp++){
              outfile.precision(16);
              outfile << fixed;
              outfile << mesh.xc[i][j][k] + 0.5*mesh.dx[i][j][k]*FR.eta[ip]  << " " 
                      << mesh.yc[i][j][k] + 0.5*mesh.dy[i][j][k]*FR.eta[jp]  << " " 
                      << mesh.zc[i][j][k] + 0.5*mesh.dz[i][j][k]*FR.eta[kp]  << " " 
                      << W[K*i+ip][K*j+jp][K*k+kp][0] << " "
                      << W[K*i+ip][K*j+jp][K*k+kp][1] << " "
                      << W[K*i+ip][K*j+jp][K*k+kp][2] << " "
                      << W[K*i+ip][K*j+jp][K*k+kp][3] << " "
                      << W[K*i+ip][K*j+jp][K*k+kp][4] << endl;
            }
          }
        }
      }
      outfile << endl;
    }
    outfile << endl;
  }
  outfile << "VARIABLES = \"x\"\n";
  outfile << "\"y\"\n";
  outfile << "\"z\"\n";
  outfile << "\"rho\"\n";
  outfile << "\"u\"\n";
  outfile << "\"v\"\n";
  outfile << "\"w\"\n";
  outfile << "\"p\"\n";
  outfile << "ZONE T = \"Exact Solution\"\n";
  outfile << "STRANDID=0, SOLUTIONTIME=0\n";
  outfile << "I=" << mesh.nx*K << ", J=" << mesh.ny*K << ", K=" << mesh.nz*K << ", ZONETYPE=Ordered\n";
  outfile << "DATAPACKING=POINT\n";
  outfile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n";
  for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      for (int j=mesh.ng; j<mesh.ny+mesh.ng; j++){
        for (int jp=0; jp<K; jp++){
          for (int k=mesh.ng; k<mesh.nz+mesh.ng; k++){
            for (int kp=0; kp<K; kp++){
              outfile.precision(16);
              outfile << fixed;
              outfile << mesh.xc[i][j][k] + 0.5*mesh.dx[i][j][k]*FR.eta[ip]  << " " 
                      << mesh.yc[i][j][k] + 0.5*mesh.dy[i][j][k]*FR.eta[jp]  << " " 
                      << mesh.zc[i][j][k] + 0.5*mesh.dz[i][j][k]*FR.eta[kp]  << " " 
                      << W0[K*i+ip][K*j+jp][K*k+kp][0] << " "
                      << W0[K*i+ip][K*j+jp][K*k+kp][1] << " "
                      << W0[K*i+ip][K*j+jp][K*k+kp][2] << " "
                      << W0[K*i+ip][K*j+jp][K*k+kp][3] << " "
                      << W0[K*i+ip][K*j+jp][K*k+kp][4] << endl;
            }
          }
        }
      }
      outfile << endl;
    }
    outfile << endl;
  }
  outfile.close();
}
