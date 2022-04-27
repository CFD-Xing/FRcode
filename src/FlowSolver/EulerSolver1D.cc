///////////////////////////////////////////////////////////////////////////////
/// \file EulerSolver1D.cc
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
// EULER 1D
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// EulerInitialCondition
// **************************************************************************//
void Solver1D::EulerInitialCondition(const int& option, double& t){

  // INPUTS
  // ----------------------
  // option  - Case option


  // OUTPUTS
  // ----------------------
  // t       - Run time


  // VARIABLE DECLARATION
  // ----------------------
  double rhol, ul, pl;
  double rhor, ur, pr;
  double rho, u, p;
  double xm;
  double g=Euler1D.getHeatRatio();


  // INITIAL CONDITION
  // ----------------------
  // Set initial condition
  BCtype = 3; // No BCs
  if (option<=8){
    switch (option){
    case 1:
      rhol = 2.281; ul = 164.83; pl = 201170.0;
      rhor = 1.408; ur = 0.0;    pr = 101100.0;
      xm = 2.0; t = 0.012;
      break;
    case 2:
      rhol = 1.408; ul = 0.0;    pl = 101100.0;
      rhor = 2.281; ur =-164.83; pr = 201170.0;
      xm = 8.0; t = 0.012;
      break;
    case 3:
      rhol = 1.045; ul = 200.00; pl = 300000.0;
      rhor = 3.483; ur = 200.00; pr = 300000.0;
      xm = 2.0; t = 0.025;
      break;
    case 4:
      rhol = 3.483; ul =-200.00; pl = 300000.0;
      rhor = 1.045; ur =-200.00; pr = 300000.0;
      xm = 8.0; t = 0.025;
      break;
    case 5:
      rhol = 1.598; ul =-383.64; pl =  91880.0;
      rhor = 2.787; ur =-216.97; pr = 200000.0;
      xm = 5.0; t = 0.035;
      break;
    case 6:
      rhol = 2.787; ul = 216.97; pl = 200000.0;
      rhor = 1.598; ur = 383.64; pr =  91880.0;
      xm = 5.0; t = 0.035;
      break;
    case 7:
      rhol = 4.696; ul = 0.0;    pl = 404400.0;
      rhor = 1.408; ur = 0.0;    pr = 101100.0;
      xm = 5.0; t = 0.007;
      break;
    case 8:
      rhol = 1.408; ul = 0.0;    pl = 101100.0;
      rhor = 4.696; ur = 0.0;    pr = 404400.0;
      xm = 5.0; t = 0.007;
      break;
    default:
      throw invalid_argument("Unknown initial condition"); 
    }
    for (int i=0; i<mesh.nx+2*mesh.ng; i++){
      for (int ip=0; ip<K; ip++){
        if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] < xm){
          W[K*i+ip][0]=rhol; 
          W[K*i+ip][1]=ul; 
          W[K*i+ip][2]=pl;
          Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
        }else{
          W[K*i+ip][0]=rhor; 
          W[K*i+ip][1]=ur; 
          W[K*i+ip][2]=pr;
          Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
        }
      }
    }
  }
  /*switch (option){
    case 9:
      rhol = 4.696; ul = 0.0;    pl = 404400.0;
      rhor = 1.744; ur = 311.689;pr = 101100.0;
      xm = 5.0; t = 0.014;
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] < 6.0){
            W[K*i+ip][0]=rhol; 
            W[K*i+ip][1]=ul; 
            W[K*i+ip][2]=pl;
            Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
          }else if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] > 8.5){
            W[K*i+ip][0]=rhor; 
            W[K*i+ip][1]=ur; 
            W[K*i+ip][2]=pr;
            Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
          }else{
            p = (pr+pl)/2.0 - (pr-pl)/2.0*cos(2.0*M_PI*(mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip]-6.0)/5.0)
            rho = rhor*pow((p/pr),(1.0/g));
            u =  2.0/(g-1.0)*sqrt(g*pl/rhol)*(1.0 - pow((p/pl),(g-1.0)/(2.0*g)));
            W[K*i+ip][0]=rho; 
            W[K*i+ip][1]=u; 
            W[K*i+ip][2]=p;
            Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
          }
        }
      }
      break;
    case 10:
      rhol = 1.744; ul =-311.689;pl = 101100.0;
      rhor = 4.696; ur = 0.0;    pr = 404400.0;
      xm = 5.0; t = 0.014;
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] > 4.0){
            W[K*i+ip][0]=rhor; 
            W[K*i+ip][1]=ur; 
            W[K*i+ip][2]=pr;
            Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
          }else if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] < 1.5){
            W[K*i+ip][0]=rhol; 
            W[K*i+ip][1]=ul; 
            W[K*i+ip][2]=pl;
            Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
          }else{
            p = (pr+pl)/2.0 - (pl-pr)/2.0*cos(2.0*M_PI*(mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip]-4.0)/5.0);
            rho = rhol*pow((p/pl),(1.0/g));
            u =  2.0/(g-1.0)*sqrt(g*pr/rhor)*(1.0 - pow((p/pr),(g-1.0)/(2.0*g)));
            W[K*i+ip][0]=rho; 
            W[K*i+ip][1]=u; 
            W[K*i+ip][2]=p;
            Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
          }
        }
      }
      break;
    default:
      throw invalid_argument("Unknown initial condition"); 
  }*/ 
  // Save initial solution
  for (int i=0; i<mesh.nx+2*mesh.ng; i++){
    for (int ip=0; ip<K; ip++){
      for (int n=0; n<nvar; n++){
        U0[K*i+ip][n] = U[K*i+ip][n];
        W0[K*i+ip][n] = W[K*i+ip][n];
      }
    }
  }
}


// **************************************************************************//
// EulerExactSolution
// **************************************************************************//
void Solver1D::EulerExactSolution(const int& option){

  // INPUTS
  // ----------------------
  // option  - Case option


  // VARIABLE DECLARATION
  // ----------------------
  double rhol, ul, pl;
  double rhor, ur, pr;
  double c_pos, c_neg, rho, u, p, c;
  double xm, x0, t;
  double df, dp_pl, err, cl, cr;
  double g=Euler1D.getHeatRatio();
  int iter;


  // EXACT SOLUTION
  // ----------------------
  if (option<=8){
    switch (option){
    case 1:
      rhol = 2.281; ul = 164.83; pl = 201170.0;
      rhor = 1.408; ur = 0.0;    pr = 101100.0;
      xm = 7.1723; t = 0.012;
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] < xm){
            Wn[K*i+ip][0]=rhol; 
            Wn[K*i+ip][1]=ul; 
            Wn[K*i+ip][2]=pl;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else{
            Wn[K*i+ip][0]=rhor; 
            Wn[K*i+ip][1]=ur; 
            Wn[K*i+ip][2]=pr;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }
        }
      }
      break;
    case 2:
      rhol = 1.408; ul = 0.0;    pl = 101100.0;
      rhor = 2.281; ur =-164.83; pr = 201170.0;
      xm = 2.8277; t = 0.012;
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] < xm){
            Wn[K*i+ip][0]=rhol; 
            Wn[K*i+ip][1]=ul; 
            Wn[K*i+ip][2]=pl;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else{
            Wn[K*i+ip][0]=rhor; 
            Wn[K*i+ip][1]=ur; 
            Wn[K*i+ip][2]=pr;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }
        }
      }
      break;
    case 3:
      rhol = 1.045; ul = 200.00; pl = 300000.0;
      rhor = 3.483; ur = 200.00; pr = 300000.0;
      xm = 7.0; t = 0.025;
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] < xm){
            Wn[K*i+ip][0]=rhol; 
            Wn[K*i+ip][1]=ul; 
            Wn[K*i+ip][2]=pl;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else{
            Wn[K*i+ip][0]=rhor; 
            Wn[K*i+ip][1]=ur; 
            Wn[K*i+ip][2]=pr;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }
        }
      }
      break;
    case 4:
      rhol = 3.483; ul =-200.00; pl = 300000.0;
      rhor = 1.045; ur =-200.00; pr = 300000.0;
      xm = 3.0; t = 0.025;
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] < xm){
            Wn[K*i+ip][0]=rhol; 
            Wn[K*i+ip][1]=ul; 
            Wn[K*i+ip][2]=pl;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else{
            Wn[K*i+ip][0]=rhor; 
            Wn[K*i+ip][1]=ur; 
            Wn[K*i+ip][2]=pr;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }
        }
      }
      break;
    case 5:
      rhol = 1.598; ul =-383.64; pl =  91880.0;
      rhor = 2.787; ur =-216.97; pr = 200000.0;
      xm = 5.0; t = 0.035;
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] <= 1.5){
            Wn[K*i+ip][0]=rhol; 
            Wn[K*i+ip][1]=ul; 
            Wn[K*i+ip][2]=pl;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] >= 8.5){
            Wn[K*i+ip][0]=Wn[K*i+ip-1][0]; 
            Wn[K*i+ip][1]=Wn[K*i+ip-1][1];
            Wn[K*i+ip][2]=Wn[K*i+ip-1][2];
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else{
            c_pos = (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip]-xm)/0.035;
            u = (g-1.0)/(g+1.0)*(2.0/(g-1.0)*(c_pos-sqrt(g*pr/rhor)) + ur);
            c = (g-1.0)/(g+1.0)*(2.0/(g-1.0)*sqrt(g*pr/rhor) - ur + c_pos);
            p = pr*pow(c/sqrt(g*pr/rhor),2.0*g/(g-1.0));
            Wn[K*i+ip][0]=g*p/(c*c); 
            Wn[K*i+ip][1]=u; 
            Wn[K*i+ip][2]=p;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }
        }
      }
      break;
    case 6:
      rhol = 2.787; ul = 216.97; pl = 200000.0;
      rhor = 1.598; ur = 383.64; pr =  91880.0;
      xm = 5.0; t = 0.035;
      for (int i=mesh.nx+2*mesh.ng-1; i>=0; i--){
        for (int ip=K-1; ip>=0; ip--){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] >= 8.5){
            Wn[K*i+ip][0]=rhor; 
            Wn[K*i+ip][1]=ur; 
            Wn[K*i+ip][2]=pr;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] <= 1.5){
            Wn[K*i+ip][0]=Wn[K*i+ip+1][0]; 
            Wn[K*i+ip][1]=Wn[K*i+ip+1][1];
            Wn[K*i+ip][2]=Wn[K*i+ip+1][2];
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else{
            c_neg = (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip]-xm)/0.035;
            u = (g-1.0)/(g+1.0)*(2.0/(g-1.0)*(c_neg+sqrt(g*pl/rhol)) + ul);
            c = (g-1.0)/(g+1.0)*(2.0/(g-1.0)*sqrt(g*pl/rhol) + ul - c_neg);
            p = pl*pow(c/sqrt(g*pl/rhol),2.0*g/(g-1.0));
            Wn[K*i+ip][0]=g*p/(c*c); 
            Wn[K*i+ip][1]=u; 
            Wn[K*i+ip][2]=p;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }
        }
      }
      break;
    case 7:
      rhol = 4.696; ul = 0.0;    pl = 404400.0;
      rhor = 1.408; ur = 0.0;    pr = 101100.0;
      xm = 5.0; t = 0.007;
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] <= 2.569455){
            Wn[K*i+ip][0]=rhol; 
            Wn[K*i+ip][1]=ul; 
            Wn[K*i+ip][2]=pl;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] >= 8.01725){
            Wn[K*i+ip][0]=rhor; 
            Wn[K*i+ip][1]=ur;
            Wn[K*i+ip][2]=pr;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] >= 6.15393){
            Wn[K*i+ip][0]=2.28000; 
            Wn[K*i+ip][1]=164.8445;
            Wn[K*i+ip][2]=2.0115e5;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] >= 3.95417){
            Wn[K*i+ip][0]=2.8516;
            Wn[K*i+ip][1]=164.8445;
            Wn[K*i+ip][2]=2.0115e5;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else{
            c_neg = (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip]-xm)/0.007;
            u = (g-1.0)/(g+1.0)*(2.0/(g-1.0)*(c_neg+sqrt(g*pl/rhol)) + ul);
            c = (g-1.0)/(g+1.0)*(2.0/(g-1.0)*sqrt(g*pl/rhol) + ul - c_neg);
            p = pl*pow(c/sqrt(g*pl/rhol),2.0*g/(g-1.0));
            Wn[K*i+ip][0]=g*p/(c*c); 
            Wn[K*i+ip][1]=u; 
            Wn[K*i+ip][2]=p;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }
        }
      }
      break;
    case 8:
      rhol = 1.408; ul = 0.0;    pl = 101100.0;
      rhor = 4.696; ur = 0.0;    pr = 404400.0;
      xm = 5.0; t = 0.007;
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] >= 7.430545){
            Wn[K*i+ip][0]=rhor; 
            Wn[K*i+ip][1]=ur; 
            Wn[K*i+ip][2]=pr;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] <= 1.98275){
            Wn[K*i+ip][0]=rhol; 
            Wn[K*i+ip][1]=ul;
            Wn[K*i+ip][2]=pl;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] <= 3.84607){
            Wn[K*i+ip][0]=2.28000; 
            Wn[K*i+ip][1]=-164.8445;
            Wn[K*i+ip][2]=2.0115e5;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] <= 6.0458259){
            Wn[K*i+ip][0]=2.8516;
            Wn[K*i+ip][1]=-164.8445;
            Wn[K*i+ip][2]=2.0115e5;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }else{
            c_pos = (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip]-xm)/0.007;
            u = (g-1.0)/(g+1.0)*(2.0/(g-1.0)*(c_pos-sqrt(g*pr/rhor)) + ur);
            c = (g-1.0)/(g+1.0)*(2.0/(g-1.0)*sqrt(g*pr/rhor) - ur + c_pos);
            p = pr*pow(c/sqrt(g*pr/rhor),2.0*g/(g-1.0));
            Wn[K*i+ip][0]=g*p/(c*c); 
            Wn[K*i+ip][1]=u; 
            Wn[K*i+ip][2]=p;
            Euler1D.PrimtoCons(&Wn[K*i+ip][0], &Un[K*i+ip][0]);
          }
        }
      }
      break;
    /*case 9:
      rhol = 4.696; ul = 0.0;    pl = 404400.0;
      rhor = 1.744; ur = 311.689;pr = 101100.0;
      xm = 5.0; t = 0.014;
      cl = sqrt(g*pl/rhol);
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] < 1.1389){
            W[K*i+ip][0]=rhol; 
            W[K*i+ip][1]=ul; 
            W[K*i+ip][2]=pl;
            Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
          }else if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] > 8.8753){
            W[K*i+ip][0]=rhor; 
            W[K*i+ip][1]=ur; 
            W[K*i+ip][2]=pr;
            Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
          }else{
            x0 = (8.5-6.0)/(8.8753-1.1389)*(mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip]-1.1389) + 6.0;
            err = (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip]-x0)/t - cl/(g-1.0)*(2.0-(g+1.0)*pow( (pr/pl+1.0)/2.0 - (pr/pl-1.0)/2.0*cos(2.0*M_PI*(x0-6.0)/5.0) ,(g-1.0)/(2.0*g)));
            iter = 0;
            while (abs(err)>1e-5 && iter<10){
              dp_pl = (pr/pl-1.0)/2.0*sin(2.0*M_PI*(x0-6.0)/5.0)*2.0*M_PI/5.0;
              df = -1.0/t + cl*(g+1.0)/(2.0*g)*pow( (pr/pl+1.0)/2.0 - (pr/pl-1)/2.0*cos(2.0*M_PI*(x0-6.0)/5) , -(g+1.0)/(2.0*g))*dp_pl;
              x0 = x0 - err/df;
              err = (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip]-x0)/t - cl/(g-1.0)*(2.0-(g+1.0)*pow( (pr/pl+1.0)/2.0 - (pr/pl-1)/2.0*cos(2.0*M_PI*(x0-6.0)/5) ,(g-1.0)/(2.0*g)));
              iter++;
            }
            p = (pr+pl)/2.0 - (pr-pl)/2.0*cos(2.0*M_PI*(x0-6.0)/5.0);
            rho = rhor*pow((p/pr),(1.0/g));
            u =  2.0/(g-1.0)*sqrt(g*pl/rhol)*(1.0 - pow((p/pl),(g-1.0)/(2.0*g)));
            W[K*i+ip][0]=rho; 
            W[K*i+ip][1]=u; 
            W[K*i+ip][2]=p;
            Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
          }
        }
      }
      break;
    case 10:
      rhol = 1.744; ul =-311.689;pl = 101100.0;
      rhor = 4.696; ur = 0.0;    pr = 404400.0;
      xm = 5.0; t = 0.014;
      for (int i=0; i<mesh.nx+2*mesh.ng; i++){
        for (int ip=0; ip<K; ip++){
          if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] > 8.8611){
            W[K*i+ip][0]=rhor; 
            W[K*i+ip][1]=ur; 
            W[K*i+ip][2]=pr;
            Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
          }else if (mesh.xc[i]+0.5*mesh.dx[i]*FR.eta[ip] < 1.1247){
            W[K*i+ip][0]=rhol; 
            W[K*i+ip][1]=ul; 
            W[K*i+ip][2]=pl;
            Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
          }else{
            p = (pr+pl)/2.0 - (pl-pr)/2.0*cos(2.0*M_PI*(x0-4.0)/5.0);
            rho = rhol*pow((p/pl),(1.0/g));
            u =  2.0/(g-1.0)*sqrt(g*pr/rhor)*(1.0 - pow((p/pr),(g-1.0)/(2.0*g)));
            W[K*i+ip][0]=rho; 
            W[K*i+ip][1]=u; 
            W[K*i+ip][2]=p;
            Euler1D.PrimtoCons(&W[K*i+ip][0], &U[K*i+ip][0]);
          }
        }
      }
      break;*/
    default:
      throw invalid_argument("Unknown initial condition"); 
    }
  }
}


// **************************************************************************//
// EulerOutput
// **************************************************************************//
void Solver1D::EulerOutput(const int& option, const double& CFL, const int& nStage){

  // INPUTS
  // ----------------------
  // option  - Case option


  // VARIABLE DECLARATION
  // ----------------------
  ofstream outfile;
  char filename[256];
  double** f = new double*[nvar];
  for (int n=0; n<nvar; n++){
    f[n] = new double[K];
  }


  // OUTPUT SOLUTION
  // ----------------------
  // Exact solution
  EulerExactSolution(option);

  // Numerical solution
  if (SolverType==1){ // FVM
    sprintf(filename, "Euler1D-Case%02d-FVM-%s-Grid=%d-CFL=%0.2f-Stage=%d.dat",option,LimiterName,nx,CFL,nStage);
    outfile.open (filename);
    outfile << "VARIABLES = \"x\"\n";
    outfile << "\"rho\"\n";
    outfile << "\"u\"\n";
    outfile << "\"p\"\n";
    outfile << "ZONE T = \"FVM Solution\"\n";
    outfile << "STRANDID=0, SOLUTIONTIME=0\n";
    outfile << "I=" << mesh.nx << ", J=1, K=1, ZONETYPE=Ordered\n";
    outfile << "DATAPACKING=POINT\n";
    outfile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE)\n";
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      outfile.precision(16);
      outfile << fixed;
      outfile << mesh.xc[i] << " " 
              << W[i][0] << " "
              << W[i][1] << " "
              << W[i][2] << endl;
    }
    outfile << "VARIABLES = \"x\"\n";
    outfile << "\"rho\"\n";
    outfile << "\"u\"\n";
    outfile << "\"p\"\n";
    outfile << "ZONE T = \"Exact Solution\"\n";
    outfile << "STRANDID=0, SOLUTIONTIME=0\n";
    outfile << "I=" << mesh.nx << ", J=1, K=1, ZONETYPE=Ordered\n";
    outfile << "DATAPACKING=POINT\n";
    outfile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE)\n";
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      outfile.precision(16);
      outfile << fixed;
      outfile << mesh.xc[i] << " " 
              << Wn[i][0] << " "
              << Wn[i][1] << " "
              << Wn[i][2] << endl;
    }
  }else{ // FR
    sprintf(filename, "Euler1D-Case%02d-FR%d-%s-Grid=%d-CFL=%0.2f-Stage=%d.dat",option,K,FR.FRname,nx,CFL,nStage);
    outfile.open (filename);
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      outfile << "VARIABLES = \"x\"\n";
      outfile << "\"rho\"\n";
      outfile << "\"u\"\n";
      outfile << "\"p\"\n";
      outfile << "ZONE T = \"" << i <<"\"\n";
      outfile << "STRANDID=0, SOLUTIONTIME=0\n";
      outfile << "I=" << 11 << ", J=1, K=1, ZONETYPE=Ordered\n";
      outfile << "DATAPACKING=POINT\n";
      outfile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE)\n";
      outfile.precision(16);
      outfile << fixed;
      for (int ip=0; ip<K; ip++){
        for (int n=0; n<nvar; n++){
          f[n][ip] = W[K*i+ip][n];
        }
      }
      for (int ip=-5; ip<=5; ip++){
        outfile << mesh.xc[i] + 0.1*ip*mesh.dx[i] << " ";
        for (int n=0; n<nvar; n++){
          outfile << LagrangeInterpol1D(K, FR.eta, &f[n][0], 0.2*ip) << " ";
        }
        outfile << endl;
      }
    }
    outfile << endl;
    outfile << "VARIABLES = \"x\"\n";
    outfile << "\"rho\"\n";
    outfile << "\"u\"\n";
    outfile << "\"p\"\n";
    outfile << "ZONE T = \"Solution Points\"\n";
    outfile << "STRANDID=0, SOLUTIONTIME=0\n";
    outfile << "I=" << mesh.nx*K << ", J=1, K=1, ZONETYPE=Ordered\n";
    outfile << "DATAPACKING=POINT\n";
    outfile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE)\n";
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      for (int ip=0; ip<K; ip++){
        outfile.precision(16);
        outfile << fixed;
        outfile << mesh.xc[i] + 0.5*mesh.dx[i]*FR.eta[ip] << " " 
                << W[K*i+ip][0] << " "
                << W[K*i+ip][1] << " "
                << W[K*i+ip][2] << " " << endl;
      }
    }
    outfile << endl;
    outfile << "VARIABLES = \"x\"\n";
    outfile << "\"rho\"\n";
    outfile << "\"u\"\n";
    outfile << "\"p\"\n";
    outfile << "ZONE T = \"Exact Solution\"\n";
    outfile << "STRANDID=0, SOLUTIONTIME=0\n";
    outfile << "I=" << mesh.nx*K << ", J=1, K=1, ZONETYPE=Ordered\n";
    outfile << "DATAPACKING=POINT\n";
    outfile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE)\n";
    for (int i=mesh.ng; i<mesh.nx+mesh.ng; i++){
      for (int ip=0; ip<K; ip++){
        outfile.precision(16);
        outfile << fixed;
        outfile << mesh.xc[i] + 0.5*mesh.dx[i]*FR.eta[ip] << " " 
                << Wn[K*i+ip][0] << " "
                << Wn[K*i+ip][1] << " "
                << Wn[K*i+ip][2] << " " << endl;
      }
    }
  }
  outfile.close();

  // Deallocate
  for (int n=0; n<nvar; n++){
    if (f!=NULL) delete[] f[n];
  }
  if (f!=NULL) delete[] f;
}
