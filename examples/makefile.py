#!/usr/bin/env python

import os
# linearAdvection1d
cmd = 'rm -r linearAdvection1d'
failure = os.system(cmd)

cmd = 'g++ ' + \
      '-fpermissive ' + \
      '-o linearAdvection1d ' + \
      '../src/Mesh/Mesh.cc ' + \
      'LinearAdvection1D.cc ' + \
      '../src/FlowSolver/LinearAdvectionSolver1D.cc ' + \
      '../src/FlowSolver/EulerSolver1D.cc ' + \
      '../src/FluxRecon/FRelement.cc ' + \
      '../src/FluxRecon/Lagrange.cc ' + \
      '../src/FluxRecon/FluxCorrection.cc ' + \
      '../src/SolutionState/EulerState1D.cc ' 
failure = os.system(cmd)

# linearAdvection2d
cmd = 'rm -r linearAdvection2d'
failure = os.system(cmd)

cmd = 'g++ ' + \
      '-fpermissive ' + \
      '-o linearAdvection2d ' + \
      '../src/Mesh/Mesh.cc ' + \
      'LinearAdvection2D.cc ' + \
      '../src/FlowSolver/LinearAdvectionSolver2D.cc ' + \
      '../src/FlowSolver/EulerSolver2D.cc ' + \
      '../src/FluxRecon/FRelement.cc ' + \
      '../src/FluxRecon/Lagrange.cc ' + \
      '../src/FluxRecon/FluxCorrection.cc ' + \
      '../src/SolutionState/EulerState2D.cc ' 
failure = os.system(cmd)

# linearAdvection3d
cmd = 'rm -r linearAdvection3d'
failure = os.system(cmd)

cmd = 'g++ ' + \
      '-fpermissive ' + \
      '-o linearAdvection3d ' + \
      '../src/Mesh/Mesh.cc ' + \
      'LinearAdvection3D.cc ' + \
      '../src/FlowSolver/LinearAdvectionSolver3D.cc ' + \
      '../src/FlowSolver/EulerSolver3D.cc ' + \
      '../src/FluxRecon/FRelement.cc ' + \
      '../src/FluxRecon/Lagrange.cc ' + \
      '../src/FluxRecon/FluxCorrection.cc ' + \
      '../src/SolutionState/EulerState3D.cc ' 
failure = os.system(cmd)

# euler1d
cmd = 'rm -r euler1d'
failure = os.system(cmd)

cmd = 'g++ ' + \
      '-fpermissive ' + \
      '-o euler1d ' + \
      '../src/Mesh/Mesh.cc ' + \
      'Euler1D.cc ' + \
      '../src/FlowSolver/LinearAdvectionSolver1D.cc ' + \
      '../src/FlowSolver/EulerSolver1D.cc ' + \
      '../src/FluxRecon/FRelement.cc ' + \
      '../src/FluxRecon/Lagrange.cc ' + \
      '../src/FluxRecon/FluxCorrection.cc ' + \
      '../src/SolutionState/EulerState1D.cc ' 
failure = os.system(cmd)

# euler2d
cmd = 'rm -r euler2d'
failure = os.system(cmd)

cmd = 'g++ ' + \
      '-fpermissive ' + \
      '-o euler2d ' + \
      '../src/Mesh/Mesh.cc ' + \
      'Euler2D.cc ' + \
      '../src/FlowSolver/LinearAdvectionSolver2D.cc ' + \
      '../src/FlowSolver/EulerSolver2D.cc ' + \
      '../src/FluxRecon/FRelement.cc ' + \
      '../src/FluxRecon/Lagrange.cc ' + \
      '../src/FluxRecon/FluxCorrection.cc ' + \
      '../src/SolutionState/EulerState2D.cc ' 
failure = os.system(cmd)

# euler3d
cmd = 'rm -r euler3d'
failure = os.system(cmd)

cmd = 'g++ ' + \
      '-fpermissive ' + \
      '-o euler3d ' + \
      '../src/Mesh/Mesh.cc ' + \
      'Euler3D.cc ' + \
      '../src/FlowSolver/LinearAdvectionSolver3D.cc ' + \
      '../src/FlowSolver/EulerSolver3D.cc ' + \
      '../src/FluxRecon/FRelement.cc ' + \
      '../src/FluxRecon/Lagrange.cc ' + \
      '../src/FluxRecon/FluxCorrection.cc ' + \
      '../src/SolutionState/EulerState3D.cc ' 
failure = os.system(cmd)
