all:
	g++ -fpermissive -o solver1d  Mesh/Mesh.cc FlowSolver/Solver1D.cc  FlowSolver/LinearAdvectionSolver1D.cc FlowSolver/EulerSolver1D.cc FluxRecon/FRelement.cc FluxRecon/Lagrange.cc FluxRecon/FluxCorrection.cc SolutionState/EulerState1D.cc -I. 
	g++ -fpermissive -o solver2d  Mesh/Mesh.cc FlowSolver/Solver2D.cc  FlowSolver/LinearAdvectionSolver2D.cc FlowSolver/EulerSolver2D.cc FluxRecon/FRelement.cc FluxRecon/Lagrange.cc FluxRecon/FluxCorrection.cc SolutionState/EulerState2D.cc -I. 
	g++ -fpermissive -o solver3d  Mesh/Mesh.cc FlowSolver/Solver3D.cc  FlowSolver/LinearAdvectionSolver3D.cc FlowSolver/EulerSolver3D.cc FluxRecon/FRelement.cc FluxRecon/Lagrange.cc FluxRecon/FluxCorrection.cc SolutionState/EulerState3D.cc -I. 

solver1d: 
	g++ -fpermissive -o solver1d  Mesh/Mesh.cc FlowSolver/Solver1D.cc  FlowSolver/LinearAdvectionSolver1D.cc FlowSolver/EulerSolver1D.cc FluxRecon/FRelement.cc FluxRecon/Lagrange.cc FluxRecon/FluxCorrection.cc SolutionState/EulerState1D.cc -I. 

solver2d: 
	g++ -fpermissive -o solver2d  Mesh/Mesh.cc FlowSolver/Solver2D.cc  FlowSolver/LinearAdvectionSolver2D.cc FlowSolver/EulerSolver2D.cc FluxRecon/FRelement.cc FluxRecon/Lagrange.cc FluxRecon/FluxCorrection.cc SolutionState/EulerState2D.cc -I. 

solver3d: 
	g++ -fpermissive -o solver3d  Mesh/Mesh.cc FlowSolver/Solver3D.cc  FlowSolver/LinearAdvectionSolver3D.cc FlowSolver/EulerSolver3D.cc FluxRecon/FRelement.cc FluxRecon/Lagrange.cc FluxRecon/FluxCorrection.cc SolutionState/EulerState3D.cc -I. 

clean:
	rm -rf solver1d
	rm -rf solver2d
	rm -rf solver3d
