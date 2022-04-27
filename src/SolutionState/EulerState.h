///////////////////////////////////////////////////////////////////////////////
/// \file EulerState.h
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 24 2021
///
///////////////////////////////////////////////////////////////////////////////
#ifndef EULERSTATE_H
#define EULERSTATE_H


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <iostream>
using namespace std;


///////////////////////////////////////////////////////////////////////////////
// EULER
///////////////////////////////////////////////////////////////////////////////
class EulerState{
  protected:
  // Air properties (perfect gas assumption)
  static constexpr double g=1.4;
  static constexpr double R=287.0;
  public:
  // Accessors
  static double getHeatRatio(void){return g;}
  static double getGasConstant(void){return R;}
};


///////////////////////////////////////////////////////////////////////////////
// EULER (1D)
///////////////////////////////////////////////////////////////////////////////
class EulerState1D: public EulerState{
  public:
  // Default constructor
  EulerState1D(){}
  // Default destructor
  ~EulerState1D(){}
  // Functions
  void PrimtoCons(const double& rho, const double& u, const double& p, double* U);
  void PrimtoCons(const double* W, double* U);
  void ConstoPrim(const double* U, double& rho, double& u, double& p);
  void ConstoPrim(const double* U, double* W);
  void InviscidFlux(const double* U, double* F);
  void InviscidFlux(const double& rho, const double& u, const double& p, double* F);
  void RoeSolver(const double* Wl, const double* Wr, double* Fupw);
  void HLLESolver(const double* Wl, const double* Wr, double* Fupw);
  void EigenValue(const double& u, const double& c, double* Lambda);
  void EntropyCorrection(const double& ul, const double& cl, 
                         const double& ur, const double& cr,
                         double* Lambda);
  void LeftPrimEigenVector(const double& rho, 
                           const double& c, 
                           const int& i,
                           double* Lp);
  void LeftConsEigenVector(const double& u, 
                           const double& c, 
                           const int& i,
                           double* Rc);
  void RightPrimEigenVector(const double& rho, 
                            const double& c, 
                            const int& i,
                            double* Rp);
  void RightConsEigenVector(const double& u, 
                            const double& c, 
                            const int& i,
                            double* Rp);
  void FluxJacobian(const double* W, double** dFdU);
};


///////////////////////////////////////////////////////////////////////////////
// EULER (2D)
///////////////////////////////////////////////////////////////////////////////
class EulerState2D: public EulerState{
  public:
  // Default constructor
  EulerState2D(){}
  // Default destructor
  ~EulerState2D(){}
  // Functions
  void PrimtoCons(const double& rho, const double& u, const double&v, const double& p, double* U);
  void PrimtoCons(const double* W, double* U);
  void ConstoPrim(const double* U, double& rho, double& u, double& v, double& p);
  void ConstoPrim(const double* U, double* W);
  void InviscidFlux(const double* U, double* F);
  void InviscidFlux(const double& rho, const double& u, const double& v, const double& p, double* F);
  void RoeSolver(const double* Wl, const double* Wr, double* Fupw);
  void HLLESolver(const double* Wl, const double* Wr, double* Fupw);
  void EigenValue(const double& u, const double& c, double* Lambda);
  void EntropyCorrection(const double& ul, const double& cl, 
                         const double& ur, const double& cr,
                         double* Lambda);
  void LeftPrimEigenVector(const double& rho, 
                           const double& c, 
                           const int& i,
                           double* Lp);
  void LeftConsEigenVector(const double& u, const double& v, 
                           const double& c, 
                           const int& i,
                           double* Rc);
  void RightPrimEigenVector(const double& rho, 
                            const double& c, 
                            const int& i,
                            double* Rp);
  void RightConsEigenVector(const double& u, const double& v,
                            const double& c, 
                            const int& i,
                            double* Rp);
  void FluxJacobian(const double* W, double** dFdU);
};


///////////////////////////////////////////////////////////////////////////////
// EULER (3D)
///////////////////////////////////////////////////////////////////////////////
class EulerState3D: public EulerState{
  public:
  // Default constructor
  EulerState3D(){}
  // Default destructor
  ~EulerState3D(){}
  // Functions
  void PrimtoCons(const double& rho, const double& u, const double& v, const double& w, const double& p, double* U);
  void PrimtoCons(const double* W, double* U);
  void ConstoPrim(const double* U, double& rho, double& u, double& v, double& w, double& p);
  void ConstoPrim(const double* U, double* W);
  void InviscidFlux(const double* U, double* F);
  void InviscidFlux(const double& rho, const double& u, const double& v, const double& w, const double& p, double* F);
  void RoeSolver(const double* Wl, const double* Wr, double* Fupw);
  void HLLESolver(const double* Wl, const double* Wr, double* Fupw);
  void EigenValue(const double& u, const double& c, double* Lambda);
  void EntropyCorrection(const double& ul, const double& cl, 
                         const double& ur, const double& cr,
                         double* Lambda);
  void LeftPrimEigenVector(const double& rho, 
                           const double& c, 
                           const int& i,
                           double* Lp);
  void LeftConsEigenVector(const double& u, const double& v, const double& w, 
                           const double& c, 
                           const int& i,
                           double* Rc);
  void RightPrimEigenVector(const double& rho, 
                            const double& c, 
                            const int& i,
                            double* Rp);
  void RightConsEigenVector(const double& u, const double& v, const double& w,
                            const double& c, 
                            const int& i,
                            double* Rp);
  void FluxJacobian(const double* W, double** dFdU);
};
#endif
