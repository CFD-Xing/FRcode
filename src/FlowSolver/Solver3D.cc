// Include
#include <chrono>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;
#include "FlowSolver.h"
typedef std::chrono::high_resolution_clock Clock;

int main(int argc, char * argv[]){

  // Declares
  auto start_time = Clock::now();
  int MeshType=1, option=-1, EquationType=-1, SolverType=-1, nStage=-1, LimiterType=-1, K=4, Ktype=1, FRtype=1;
  int nx=-1, ny=-1, nz=-1;
  double CFL=-1;
  int i=1;

  // Arguments
  while (i<argc){
    if (strcmp(argv[i], "-EquationType")==0){
      i++;
      EquationType = atoi(argv[i]);
      if (EquationType!=1 && EquationType!=2){
        cout << "Error! Unknown EquationType, 1 - Linear advection, 2 - Euler" << endl;
        return 1;
      }
    }else if (strcmp(argv[i], "-Case")==0){
      i++;
      option = atoi(argv[i]);
      //if (option!=1 && option!=2){
      //  cout << "Error! Unknown Case, 1 - Smooth case, 2 - Non-smooth case" << endl;
      //  return 1;
      //}
    }else if (strcmp(argv[i], "-SolverType")==0){
      i++;
      SolverType = atoi(argv[i]);
      if (SolverType!=1 && SolverType!=2){
        cout << "Error! Unknown SolverType, 1 - FVM, 2 - FR" << endl;
        return 1;
      }
    }else if (strcmp(argv[i], "-LimiterType")==0){
      i++;
      LimiterType = atoi(argv[i]);
      if (LimiterType!=0 && LimiterType!=1 && LimiterType!=2 && LimiterType!=3){
        cout << "Error! Unknown LimiterType, 0 - Zero, 1 - Unlimited, 2 - Barth-Jespersen, 3 - Venkatakrishnan" << endl;
        return 1;
      }
    }else if (strcmp(argv[i], "-nStage")==0){
      i++;
      nStage = atoi(argv[i]);
      if (nStage<=0 || nStage>6){
        cout << "Error! nStage, 0 - 6" << endl;
        return 1;
      }
    }else if (strcmp(argv[i], "-CFL")==0){
      i++;
      CFL = atof(argv[i]);
      if (CFL<=0.0){
        cout << "Error! CFL must be larger than 0" << endl;
        return 1;
      }
    }else if (strcmp(argv[i], "-nx")==0){
      i++;
      nx = atoi(argv[i]);
      if (nx<=0){
        cout << "Error! nx must be larger than 0" << endl;
        return 1;
      }
    }else if (strcmp(argv[i], "-ny")==0){
      i++;
      ny = atoi(argv[i]);
      if (ny<=0){
        cout << "Error! ny must be larger than 0" << endl;
        return 1;
      }
    }else if (strcmp(argv[i], "-nz")==0){
      i++;
      nz = atoi(argv[i]);
      if (nz<=0){
        cout << "Error! nz must be larger than 0" << endl;
        return 1;
      }
    }else{
      cout << "Error! Unknow parameter " << argv[i] << endl;
      return 1;
    }
    i++;
  }

  // Check argument
  if(EquationType==-1){
    cout << "Error! EquationType not specified, 1 - Linear advection, 2 - Euler" << endl;
    return 1;
  }else if(option==-1){
    cout << "Error! Case not specified, 1 - Smooth case, 2 - Non-smooth case" << endl;
    return 1;
  }else if (SolverType==-1){
    cout << "Error! SolverType not specified, 1 - FVM, 2 - FR" << endl;
    return 1;
  }else if (SolverType==1 && LimiterType==-1){
    cout << "Error! LimiterType not specified, 0 - Zero, 1 - Unlimited, 2 - Barth-Jespersen, 3 - Venkatakrishnan" << endl;
    return 1;
  }else if (nStage==-1){
    cout << "Error! nStage not specified, 0 - 6" << endl;
    return 1;
  }else if (CFL==-1){
    cout << "Error! CFL not specified" << endl;
    return 1;
  }else if (nx==-1){
    cout << "Error! nx not specified" << endl;
    return 1;
  }else if (ny==-1){
    cout << "Error! ny not specified" << endl;
    return 1;
  }else if (nz==-1){
    cout << "Error! nz not specified" << endl;
    return 1;
  }

  // Other
  double t=0.0, tfinal, dt;

  // Declare 3D solver
  Solver3D solver3D(K, nx, ny, nz, SolverType, EquationType, LimiterType, Ktype, FRtype);

  // Set mesh
  solver3D.setMesh(MeshType, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0);

  // Set initial condition
  if (EquationType==1){
    solver3D.LinearAdvectionInitialCondition(option, 10.0, tfinal);
  }else if (EquationType==2){
    solver3D.EulerInitialCondition(option, tfinal);
  }

  // Loop over time
  while(t<tfinal){
    cout << "Time = " << t << endl;
    dt = solver3D.TimeStep(CFL);
    solver3D.RungeKutta(nStage);
    t+=dt;
  }

  // Output solution
  if (EquationType==1){
    solver3D.LinearAdvectionOutput(option, CFL, nStage, t);
  }else if (EquationType==2){
    solver3D.EulerOutput(option, CFL, nStage);
  }
  
  auto end_time = Clock::now();
  cout << "Total time is: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count()*1e-9 << " sec" << std::endl;
  return 0;
}
