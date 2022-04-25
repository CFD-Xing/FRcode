///////////////////////////////////////////////////////////////////////////////
/// \file FRelement.h
///
/// \author Jacques Y. Xing
///
///   The main changes are the addition of:
///
/// \date November 21 2021
///
///////////////////////////////////////////////////////////////////////////////
#ifndef FR_ELEMENT_H
#define FR_ELEMENT_H


///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////
// Include required C++ libraries
#include <iostream>
using namespace std;

// Includes
#include "FluxCorrection.h"
#include "Lagrange.h"
#include <cstring>


///////////////////////////////////////////////////////////////////////////////
// CLASS DEFINITIONS
///////////////////////////////////////////////////////////////////////////////
class FRelement1D{
  public:
  int K;
  int Ktype;
  int FRtype;
  char FRname[256];
  double* eta;
  double* Lvec;
  double* Rvec;
  double* gl_prime;
  double* gr_prime;
  double** Dmat;
  // Default constructor
  FRelement1D(){ Nullify(); }
  // Constructor
  FRelement1D(int K, int Ktype, int FRtype){
    this->K = K;
    this->Ktype = Ktype;
    this->FRtype = FRtype;
    Allocate();
    setElement();
  }
  // Default destructor
  ~FRelement1D(){ Deallocate(); Nullify(); }
  // Some functions
  void setFRelement1D(int K, int Ktype, int FRtype){
    Deallocate();
    this->K = K;
    this->Ktype = Ktype;
    this->FRtype = FRtype;
    Allocate();
    setElement();
  }
  void setSolutionPoints();
  void setFluxCorrectionDerivative();
  void setElement();
  // Allocate
  void Allocate(){
    eta = new double[K];
    Lvec = new double[K];
    Rvec = new double[K];
    gl_prime = new double[K];
    gr_prime = new double[K];
    Dmat = new double*[K];
    for (int i=0; i<K; i++){
      Dmat[i] = new double[K];
    }
  }
  // Deallocate
  void Deallocate(){
    if (eta!=NULL) delete[] eta;
    if (Lvec!=NULL) delete[] Lvec;
    if (Rvec!=NULL) delete[] Rvec;
    if (gl_prime!=NULL) delete[] gl_prime;
    if (gr_prime!=NULL) delete[] gr_prime;
    for (int i=0; i<K; i++){
      if (Dmat!=NULL) delete[] Dmat[i];
    }
    if (Dmat!=NULL) delete[] Dmat;
  }
  // Nullify
  void Nullify(){
    this->K = -1;
    this->Ktype = -1;
    this->FRtype = -1;
    eta = NULL;
    Lvec = NULL;
    Rvec = NULL;
    gl_prime = NULL;
    gr_prime = NULL;
    Dmat = NULL;
  }
};
#endif
