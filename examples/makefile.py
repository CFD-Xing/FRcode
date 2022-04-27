#!/usr/bin/env python

import os
# solver1d
cmd = 'rm -r solver1d'
failure = os.system(cmd)

cmd = 'g++ ' + \
      '-fpermissive ' + \
      '-o solver1d ' + \
      '../src/Mesh/Mesh.cc ' + \
      'Solver1D.cc ' + \
      '../src/FlowSolver/LinearAdvectionSolver1D.cc ' + \
      '../src/FlowSolver/EulerSolver1D.cc ' + \
      '../src/FluxRecon/FRelement.cc ' + \
      '../src/FluxRecon/Lagrange.cc ' + \
      '../src/FluxRecon/FluxCorrection.cc ' + \
      '../src/SolutionState/EulerState1D.cc ' 
failure = os.system(cmd)

# solver2d
cmd = 'rm -r solver2d'
failure = os.system(cmd)

cmd = 'g++ ' + \
      '-fpermissive ' + \
      '-o solver2d ' + \
      '../src/Mesh/Mesh.cc ' + \
      'Solver2D.cc ' + \
      '../src/FlowSolver/LinearAdvectionSolver2D.cc ' + \
      '../src/FlowSolver/EulerSolver2D.cc ' + \
      '../src/FluxRecon/FRelement.cc ' + \
      '../src/FluxRecon/Lagrange.cc ' + \
      '../src/FluxRecon/FluxCorrection.cc ' + \
      '../src/SolutionState/EulerState2D.cc ' 
failure = os.system(cmd)

# solver3d
cmd = 'rm -r solver3d'
failure = os.system(cmd)

cmd = 'g++ ' + \
      '-fpermissive ' + \
      '-o solver3d ' + \
      '../src/Mesh/Mesh.cc ' + \
      'Solver3D.cc ' + \
      '../src/FlowSolver/LinearAdvectionSolver3D.cc ' + \
      '../src/FlowSolver/EulerSolver3D.cc ' + \
      '../src/FluxRecon/FRelement.cc ' + \
      '../src/FluxRecon/Lagrange.cc ' + \
      '../src/FluxRecon/FluxCorrection.cc ' + \
      '../src/SolutionState/EulerState3D.cc ' 
failure = os.system(cmd)
