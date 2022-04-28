#!/bin/bash
../src/solver3d -EquationType 1 -Case 1 -SolverType 1 -LimiterType 1 -nStage 4 -CFL 0.9 -nx 32 -ny 32 -nz 32
../src/solver3d -EquationType 1 -Case 2 -SolverType 1 -LimiterType 1 -nStage 4 -CFL 0.9 -nx 32 -ny 32 -nz 32
../src/solver3d -EquationType 1 -Case 1 -SolverType 2 -nStage 4 -CFL 0.1 -nx 32 -ny 32 -nz 32
../src/solver3d -EquationType 1 -Case 2 -SolverType 2 -nStage 4 -CFL 0.1 -nx 32 -ny 32 -nz 32
../src/solver3d -EquationType 2 -Case 1 -SolverType 1 -LimiterType 2 -nStage 4 -CFL 0.9 -nx 32 -ny 32 -nz 32
../src/solver3d -EquationType 2 -Case 2 -SolverType 1 -LimiterType 2 -nStage 4 -CFL 0.9 -nx 32 -ny 32 -nz 32
../src/solver3d -EquationType 2 -Case 1 -SolverType 2 -nStage 4 -CFL 0.1 -nx 32 -ny 32 -nz 32
../src/solver3d -EquationType 2 -Case 2 -SolverType 2 -nStage 4 -CFL 0.1 -nx 32 -ny 32 -nz 32
