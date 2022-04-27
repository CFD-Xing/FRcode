#!/usr/bin/env python

import os
# solver1d
cmd = 'rm -r solver1d'
failure = os.system(cmd)

cmd = 'g++ ' + \
      '-fpermissive ' + \
      '-o solver1d ' + \
      'Mesh/Mesh.cc ' + \
      'FlowSolver/Solver1D.cc ' + \
      'FlowSolver/LinearAdvectionSolver1D.cc ' + \
      'FlowSolver/EulerSolver1D.cc ' + \
      'FluxRecon/FRelement.cc ' + \
      'FluxRecon/Lagrange.cc ' + \
      'FluxRecon/FluxCorrection.cc ' + \
      'SolutionState/EulerState1D.cc ' 
failure = os.system(cmd)

# solver2d
cmd = 'rm -r solver2d'
failure = os.system(cmd)

cmd = 'g++ ' + \
      '-fpermissive ' + \
      '-o solver2d ' + \
      'Mesh/Mesh.cc ' + \
      'FlowSolver/Solver2D.cc ' + \
      'FlowSolver/LinearAdvectionSolver2D.cc ' + \
      'FlowSolver/EulerSolver2D.cc ' + \
      'FluxRecon/FRelement.cc ' + \
      'FluxRecon/Lagrange.cc ' + \
      'FluxRecon/FluxCorrection.cc ' + \
      'SolutionState/EulerState2D.cc ' 
failure = os.system(cmd)

# solver3d
cmd = 'rm -r solver3d'
failure = os.system(cmd)

cmd = 'g++ ' + \
      '-fpermissive ' + \
      '-o solver3d ' + \
      'Mesh/Mesh.cc ' + \
      'FlowSolver/Solver3D.cc ' + \
      'FlowSolver/LinearAdvectionSolver3D.cc ' + \
      'FlowSolver/EulerSolver3D.cc ' + \
      'FluxRecon/FRelement.cc ' + \
      'FluxRecon/Lagrange.cc ' + \
      'FluxRecon/FluxCorrection.cc ' + \
      'SolutionState/EulerState3D.cc ' 
failure = os.system(cmd)
