// Include
#include <chrono>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;
#include "../src/FlowSolver/FlowSolver.h"
typedef std::chrono::high_resolution_clock Clock;

int main(void){

  auto start_time = Clock::now();
  // Case definition
  int EquationType=1; // 1 - Linear wave, 2 - Euler
  int option=2; // 1 - Case 1, 2 - Case 2 

  // Solver 
  int SolverType=2; // 1 - FVM, 2 - FR

  // Time stepping
  int nStage=4;
  double CFL=0.10;

  // Mesh
  int MeshType=1; // 1 - Uniform
  int nx=128;

  // FVM
  int LimiterType=1; // 0 - Zero, 1 - Unlimieted, 2 - Barth-Jespersen, 3 - Venkatakrishnan

  // Flux Reconstruction
  int K=4;
  int Ktype=1;
  int FRtype=1;

  // Other
  double t=0.0, tfinal, dt;

  // Declare 1D linear advection solver
  Solver1D LinearAdvection1D(K, nx, SolverType, EquationType, LimiterType, Ktype, FRtype);

  // Set mesh
  LinearAdvection1D.setMesh(MeshType, 0.0, 10.0);

  // Set initial condition
  LinearAdvection1D.LinearAdvectionInitialCondition(option, 10.0, tfinal);

  // Loop over time
  while(t<tfinal){
    cout << "Time = " << t << endl;
    dt = LinearAdvection1D.TimeStep(CFL);
    LinearAdvection1D.RungeKutta(nStage);
    t+=dt;
  }

  // Output solution
  LinearAdvection1D.LinearAdvectionOutput(option, CFL, nStage, t);
  
  auto end_time = Clock::now();
  cout << "Total time is: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()*1e-9 << " sec" << std::endl;
  return 0;
}
