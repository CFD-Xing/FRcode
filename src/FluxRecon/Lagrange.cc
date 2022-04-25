///////////////////////////////////////////////////////////////////////////////
/// \file Lagrange.cc
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 20 2021
///
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <iostream>
using namespace std;

// Include
#include "Lagrange.h"


///////////////////////////////////////////////////////////////////////////////
// LAGRANGE POLYNOMIAL INTERPOLATION
///////////////////////////////////////////////////////////////////////////////
// **************************************************************************//
// LagrangeBasis1D
// **************************************************************************//
void LagrangeBasis1D(const int& nx,
                     const double* xi,
                     const double& xp,
                     double* li){

  // DESCRIPTION
  // ----------------------
  // Compute the 1D Lagrange interpolation basis functions, li, at interpolation 
  // point xp.


  // INPUTS
  // ----------------------
  // nx      - Number of interpolation data points
  // xi      - Interpolation data points
  // xp      - Interpolation point


  // OUTPUTS
  // ----------------------
  // li      - Lagrange basis functions at interpolation point xp


  // VARIABLE DECLARATION
  // ----------------------
  double tmp;


  // 1D LAGRANGIAN INTERP.
  // ----------------------
  for (int i=0; i<nx; i++){
    tmp = 1.0;
    for (int m=0; m<nx; m++){
      if (m!=i) tmp*=(xp-xi[m])/(xi[i]-xi[m]);
    }
    li[i] = tmp;
  }
}


// **************************************************************************//
// LagrangeBasisDerivative1D
// **************************************************************************//
void LagrangeBasisDerivative1D(const int& nx,
                               const double* xi,
                               const double& xp,
                               double* dlidx){

  // DESCRIPTION
  // ----------------------
  // Compute the 1D Lagrange interpolation basis function derivatives, dlidx, at 
  // interpolation point xp.


  // INPUTS
  // ----------------------
  // nx      - Number of interpolation data points
  // xi      - Interpolation data points
  // xp      - Interpolation point


  // OUTPUTS
  // ----------------------
  // dlidx   - Lagrange basis derivatives at interpolation point xp


  // VARIABLE DECLARATION
  // ----------------------
  double tmp;


  // 1D LAGRANGIAN INTERP.
  // ----------------------
  for (int i=0; i<nx; i++){
    dlidx[i] = 0.0;
    for (int m=0; m<nx; m++){
      if (m==i) continue;
      tmp = 1.0;
      for (int n=0; n<nx; n++){
        if (n!=m && n!=i) tmp*=(xp-xi[n])/(xi[i]-xi[n]);
      }
      dlidx[i] += tmp/(xi[i]-xi[m]);
    }
  }
}


// **************************************************************************//
// LagrangeInterpol1D
// **************************************************************************//
double LagrangeInterpol1D(const int& nx,
                          const double* xi,
                          const double* f, 
                          const double& xp){

  // DESCRIPTION
  // ----------------------
  // Compute the 1D Lagrange interpolation, fp, at interpolation point xp.


  // INPUTS
  // ----------------------
  // nx      - Number of interpolation data points
  // xi      - Interpolation data points
  // f       - Interpolation data points (f=f(x))
  // xp      - Interpolation point


  // OUTPUTS
  // ----------------------
  // fp      - Interpolated value at point xp


  // VARIABLE DECLARATION
  // ----------------------
  double fp, tmp;


  // 1D LAGRANGIAN INTERP.
  // ----------------------
  fp = 0.0;
  for (int i=0; i<nx; i++){
    tmp = f[i];
    for (int m=0; m<nx; m++){
      if (m!=i) tmp*=(xp-xi[m])/(xi[i]-xi[m]);
    }
    fp += tmp;
  }


  // RETURN VARIABLE
  // ----------------------
  return fp;
}


// **************************************************************************//
// LagrangeBasis2D
// **************************************************************************//
void LagrangeBasis2D(const int& nx,
                     const int& ny,
                     const double* xi, 
                     const double* yj,
                     const double& xp,
                     const double& yp,
                     double** lij){

  // DESCRIPTION
  // ----------------------
  // Compute the 2D Lagrange interpolation basis functions, lij, at interpolation 
  // point xp, yp.


  // INPUTS
  // ----------------------
  // nx      - Number of data points
  // ny      - Number of data points
  // xi      - Interpolation data points
  // yj      - Interpolation data points
  // xp      - Interpolation point
  // yp      - Interpolation point


  // OUTPUTS
  // ----------------------
  // lij      - Lagrange basis functions at interpolation point xp, yp


  // VARIABLE DECLARATION
  // ----------------------
  double tmp;


  // 2D LAGRANGIAN INTERP.
  // ----------------------
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      tmp = 1.0;
      for (int m=0; m<nx; m++){
        if (m!=i) tmp*=(xp-xi[m])/(xi[i]-xi[m]);
      }
      for (int m=0; m<ny; m++){
        if (m!=j) tmp*=(yp-yj[m])/(yj[j]-yj[m]);
      }
      lij[i][j] = tmp;
    }
  }
}


// **************************************************************************//
// LagrangeBasisDerivative2D
// **************************************************************************//
void LagrangeBasisDerivative2D(const int& nx,
                               const int& ny,
                               const double* xi, 
                               const double* yj,
                               const double& xp,
                               const double& yp,
                               double** dlijdx,
                               double** dlijdy){

  // DESCRIPTION
  // ----------------------
  // Compute the 2D Lagrange interpolation basis function derivatives, dlijdx 
  // and dlijdy, at interpolation point xp, yp.


  // INPUTS
  // ----------------------
  // nx      - Number of data points
  // ny      - Number of data points
  // xi      - Interpolation data points
  // yj      - Interpolation data points
  // xp      - Interpolation point
  // yp      - Interpolation point


  // OUTPUTS
  // ----------------------
  // dlijdx  - Lagrange basis derivatives at interpolation point xp, yp
  // dlijdy  - Lagrange basis derivatives at interpolation point xp, yp


  // VARIABLE DECLARATION
  // ----------------------
  double tmp;


  // 2D LAGRANGIAN INTERP.
  // ----------------------
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      // dlijdx
      dlijdx[i][j] = 0.0;
      for (int m=0; m<nx; m++){
        if (m==i) continue;
        tmp = 1.0;
        for (int n=0; n<nx; n++){
          if (n!=m && n!=i) tmp*=(xp-xi[n])/(xi[i]-xi[n]);
        }
        dlijdx[i][j] += tmp/(xi[i]-xi[m]);
      }
      for (int m=0; m<ny; m++){
        if (m!=j) dlijdx[i][j]*=(yp-yj[m])/(yj[j]-yj[m]);
      }
      // dlijdy
      dlijdy[i][j] = 0.0;
      for (int m=0; m<ny; m++){
        if (m==j) continue;
        tmp = 1.0;
        for (int n=0; n<ny; n++){
          if (n!=m && n!=j) tmp*=(yp-yj[n])/(yj[j]-yj[n]);
        }
        dlijdy[i][j] += tmp/(yj[j]-yj[m]);
      }
      for (int m=0; m<nx; m++){
        if (m!=i) dlijdy[i][j]*=(xp-xi[m])/(xi[i]-xi[m]);
      }
    }
  }
}


// **************************************************************************//
// LagrangeInterpol2D
// **************************************************************************//
double LagrangeInterpol2D(const int& nx,
                          const int& ny,
                          const double* xi, 
                          const double* yj,
                          const double** f,
                          const double& xp,
                          const double& yp){

  // DESCRIPTION
  // ----------------------
  // Compute the 2D Lagrange interpolation, fp, at interpolation point xp, yp.


  // INPUTS
  // ----------------------
  // nx      - Number of data points
  // ny      - Number of data points
  // xi      - Interpolation data points
  // yj      - Interpolation data points
  // f       - Interpolation data points (f=f(x,y))
  // xp      - Interpolation point
  // yp      - Interpolation point


  // OUTPUTS
  // ----------------------
  // fp      - Interpolated value at point xp, yp


  // VARIABLE DECLARATION
  // ----------------------
  double fp, tmp;


  // 2D LAGRANGIAN INTERP.
  // ----------------------
  fp = 0.0;
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      tmp = f[i][j];
      for (int m=0; m<nx; m++){
        if (m!=i) tmp*=(xp-xi[m])/(xi[i]-xi[m]);
      }
      for (int m=0; m<ny; m++){
        if (m!=j) tmp*=(yp-yj[m])/(yj[j]-yj[m]);
      }
      fp += tmp;
    }
  }


  // RETURN VARIABLE
  // ----------------------
  return fp;
}


// **************************************************************************//
// LagrangeBasis3D
// **************************************************************************//
void LagrangeBasis3D(const int& nx,
                     const int& ny,
                     const int& nz,
                     const double* xi, 
                     const double* yj,
                     const double* zk,
                     const double& xp,
                     const double& yp,
                     const double& zp,
                     double*** lijk){

  // DESCRIPTION
  // ----------------------
  // Compute the 3D Lagrange interpolation basis functions, lijk, at interpolation 
  // point xp, yp, zp.


  // INPUTS
  // ----------------------
  // nx      - Number of data points
  // ny      - Number of data points
  // nz      - Number of data points
  // xi      - Interpolation data points
  // yj      - Interpolation data points
  // zk      - Interpolation data points
  // xp      - Interpolation point
  // yp      - Interpolation point
  // zp      - Interpolation point


  // OUTPUTS
  // ----------------------
  // lijk    - Lagrange basis functions at interpolation point xp, yp, zp


  // VARIABLE DECLARATION
  // ----------------------
  double tmp;


  // 3D LAGRANGIAN INTERP.
  // ----------------------
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
        tmp=1.0;
        for (int m=0; m<nx; m++){
          if (m!=i) tmp*=(xp-xi[m])/(xi[i]-xi[m]);
        }
        for (int m=0; m<ny; m++){
          if (m!=j) tmp*=(yp-yj[m])/(yj[j]-yj[m]);
        }
        for (int m=0; m<nz; m++){
          if (m!=k) tmp*=(zp-zk[m])/(zk[k]-zk[m]);
        }
        lijk[i][j][k] = tmp;
      }
    }
  }
}


// **************************************************************************//
// LagrangeBasisDerivative3D
// **************************************************************************//
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
                               double*** dlijkdz){

  // DESCRIPTION
  // ----------------------
  // Compute the 3D Lagrange interpolation basis function derivatives, dlijkdx,
  // dlijkdy, and dlijkdz at interpolation point xp, yp, zp.


  // INPUTS
  // ----------------------
  // nx      - Number of data points
  // ny      - Number of data points
  // nz      - Number of data points
  // xi      - Interpolation data points
  // yj      - Interpolation data points
  // zk      - Interpolation data points
  // xp      - Interpolation point
  // yp      - Interpolation point
  // zp      - Interpolation point


  // OUTPUTS
  // ----------------------
  // dlijkdx - Lagrange basis derivatives at interpolation point xp, yp, zp
  // dlijkdy - Lagrange basis derivatives at interpolation point xp, yp, zp
  // dlijkdz - Lagrange basis derivatives at interpolation point xp, yp, zp


  // VARIABLE DECLARATION
  // ----------------------
  double tmp;


  // 3D LAGRANGIAN INTERP.
  // ----------------------
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
        //dlijkdx
        dlijkdx[i][j][k]=0.0;
        for (int m=0; m<nx; m++){
          if (m==i) continue;
          tmp = 1.0;
          for (int n=0; n<nx; n++){
            if (n!=m && n!=i) tmp*=(xp-xi[n])/(xi[i]-xi[n]);
          }
          dlijkdx[i][j][k] += tmp/(xi[i]-xi[m]);
        }
        for (int m=0; m<ny; m++){
          if (m!=j) dlijkdx[i][j][k]*=(yp-yj[m])/(yj[j]-yj[m]);
        }
        for (int m=0; m<nz; m++){
          if (m!=k) dlijkdx[i][j][k]*=(zp-zk[m])/(zk[k]-zk[m]);
        }
        //dlijkdy
        dlijkdy[i][j][k]=0.0;
        for (int m=0; m<ny; m++){
          if (m==j) continue;
          tmp = 1.0;
          for (int n=0; n<ny; n++){
            if (n!=m && n!=j) tmp*=(yp-yj[n])/(yj[j]-yj[n]);
          }
          dlijkdy[i][j][k] += tmp/(yj[j]-yj[m]);
        }
        for (int m=0; m<nx; m++){
          if (m!=i) dlijkdy[i][j][k]*=(xp-xi[m])/(xi[i]-xi[m]);
        }
        for (int m=0; m<nz; m++){
          if (m!=k) dlijkdy[i][j][k]*=(zp-zk[m])/(zk[k]-zk[m]);
        }
        //dlijkdz
        dlijkdz[i][j][k]=0.0;
        for (int m=0; m<nz; m++){
          if (m==k) continue;
          tmp = 1.0;
          for (int n=0; n<nz; n++){
            if (n!=m && n!=k) tmp*=(zp-zk[n])/(zk[k]-zk[n]);
          }
          dlijkdz[i][j][k] += tmp/(zk[k]-zk[m]);
        }
        for (int m=0; m<nx; m++){
          if (m!=i) dlijkdz[i][j][k]*=(xp-xi[m])/(xi[i]-xi[m]);
        }
        for (int m=0; m<ny; m++){
          if (m!=j) dlijkdz[i][j][k]*=(yp-yj[m])/(yj[j]-yj[m]);
        }
      }
    }
  }
}


// **************************************************************************//
// LagrangeInterpol3D
// **************************************************************************//
double LagrangeInterpol3D(const int& nx,
                          const int& ny,
                          const int& nz,
                          const double* xi, 
                          const double* yj,
                          const double* zk,
                          const double*** f,
                          const double& xp,
                          const double& yp,
                          const double& zp){

  // DESCRIPTION
  // ----------------------
  // Compute the 3D Lagrange interpolation, fp, at interpolation point xp, yp, zp.


  // INPUTS
  // ----------------------
  // nx      - Number of data points
  // ny      - Number of data points
  // nz      - Number of data points
  // xi      - Interpolation data points
  // yj      - Interpolation data points
  // zk      - Interpolation data points
  // f       - Interpolation data points (f=f(x,y,z))
  // xp      - Interpolation point
  // yp      - Interpolation point
  // zp      - Interpolation point


  // OUTPUTS
  // ----------------------
  // fp      - Interpolated value at point xp, yp, zp


  // VARIABLE DECLARATION
  // ----------------------
  double fp, tmp;


  // 3D LAGRANGIAN INTERP.
  // ----------------------
  fp = 0.0;
  for (int i=0; i<nx; i++){
    for (int j=0; j<ny; j++){
      for (int k=0; k<nz; k++){
        tmp = f[i][j][k];
        for (int m=0; m<nx; m++){
          if (m!=i) tmp*=(xp-xi[m])/(xi[i]-xi[m]);
        }
        for (int m=0; m<ny; m++){
          if (m!=j) tmp*=(yp-yj[m])/(yj[j]-yj[m]);
        }
        for (int m=0; m<nz; m++){
          if (m!=k) tmp*=(zp-zk[m])/(zk[k]-zk[m]);
        }
        fp += tmp;
      }
    }
  }


  // RETURN VARIABLE
  // ----------------------
  return fp;
}
