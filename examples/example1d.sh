#!/bin/bash
../src/solver1d -EquationType 1 -Case 1 -SolverType 1 -LimiterType 1 -nStage 4 -CFL 0.9 -nx 128
../src/solver1d -EquationType 1 -Case 2 -SolverType 1 -LimiterType 1 -nStage 4 -CFL 0.9 -nx 128
../src/solver1d -EquationType 1 -Case 1 -SolverType 2 -nStage 4 -CFL 0.1 -nx 128
../src/solver1d -EquationType 1 -Case 2 -SolverType 2 -nStage 4 -CFL 0.1 -nx 128
../src/solver1d -EquationType 2 -Case 1 -SolverType 1 -LimiterType 2 -nStage 4 -CFL 0.9 -nx 128
../src/solver1d -EquationType 2 -Case 2 -SolverType 1 -LimiterType 2 -nStage 4 -CFL 0.9 -nx 128
../src/solver1d -EquationType 2 -Case 3 -SolverType 1 -LimiterType 2 -nStage 4 -CFL 0.9 -nx 128
../src/solver1d -EquationType 2 -Case 4 -SolverType 1 -LimiterType 2 -nStage 4 -CFL 0.9 -nx 128
../src/solver1d -EquationType 2 -Case 5 -SolverType 1 -LimiterType 2 -nStage 4 -CFL 0.9 -nx 128
../src/solver1d -EquationType 2 -Case 6 -SolverType 1 -LimiterType 2 -nStage 4 -CFL 0.9 -nx 128
../src/solver1d -EquationType 2 -Case 7 -SolverType 1 -LimiterType 2 -nStage 4 -CFL 0.9 -nx 128
../src/solver1d -EquationType 2 -Case 8 -SolverType 1 -LimiterType 2 -nStage 4 -CFL 0.9 -nx 128
../src/solver1d -EquationType 2 -Case 1 -SolverType 2 -nStage 4 -CFL 0.1 -nx 128
../src/solver1d -EquationType 2 -Case 2 -SolverType 2 -nStage 4 -CFL 0.1 -nx 128
../src/solver1d -EquationType 2 -Case 3 -SolverType 2 -nStage 4 -CFL 0.1 -nx 128
../src/solver1d -EquationType 2 -Case 4 -SolverType 2 -nStage 4 -CFL 0.1 -nx 128
../src/solver1d -EquationType 2 -Case 5 -SolverType 2 -nStage 4 -CFL 0.1 -nx 128
../src/solver1d -EquationType 2 -Case 6 -SolverType 2 -nStage 4 -CFL 0.1 -nx 128
../src/solver1d -EquationType 2 -Case 7 -SolverType 2 -nStage 4 -CFL 0.1 -nx 128
../src/solver1d -EquationType 2 -Case 8 -SolverType 2 -nStage 4 -CFL 0.1 -nx 128
