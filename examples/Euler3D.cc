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
  int EquationType=2; // 1 - Linear wave, 2 - Euler
  int option=1; // 1 - Case 1, 2 - Case 2 

  // Solver 
  int SolverType=2; // 1 - FVM, 2 - FR

  // Time stepping
  int nStage=4;
  double CFL=0.10;

  // Mesh
  int MeshType=1; // 1 - Uniform
  int nx=32;
  int ny=32;
  int nz=32;

  // FVM
  int LimiterType=3; // 0 - Zero, 1 - Unlimieted, 2 - Barth-Jespersen, 3 - Venkatakrishnan

  // Flux Reconstruction
  int K=4;
  int Ktype=1;
  int FRtype=1;

  // Other
  double t=0.0, tfinal, dt;

  // Declare 3D Euler solver
  Solver3D Euler3D(K, nx, ny, nz, SolverType, EquationType, LimiterType, Ktype, FRtype);

  // Set mesh
  Euler3D.setMesh(MeshType, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0);

  // Set initial condition
  Euler3D.EulerInitialCondition(option, tfinal);

  // Loop over time
  while(t<tfinal){
    cout << "Time = " << t << endl;
    dt = Euler3D.TimeStep(CFL);
    Euler3D.RungeKutta(nStage);
    t+=dt;
  }

  // Output solution
  Euler3D.EulerOutput(option, CFL, nStage);
  
  auto end_time = Clock::now();
  cout << "Total time is: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()*1e-9 << " sec" << std::endl;
  return 0;
}
