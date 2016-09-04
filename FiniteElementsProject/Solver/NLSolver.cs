using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class NLSolver
    {
        private double[] residual;
        double residualNorm;
        private double[] deltaU;
        public double[] internalForcesTotalVector;
        private int maxIterations;
        public double tolerance = 1e-5;
        private int[] boundaryDof;
        public double[] solutionVector;
        private double[] forceVector;
        private Discretization2DFrame discretization;

        public NLSolver(Discretization2DFrame discretization, int maxIterations, InputData inputData)
        {
            this.discretization = discretization;
            this.forceVector = inputData.externalForcesVector;
            solutionVector = new double[inputData.nodesX.Length*3];
            this.maxIterations = maxIterations;
            boundaryDof = inputData.boundaryDof;
            this.residual = new double[inputData.nodesX.Length * 3 - boundaryDof.Length];
        }
        public void NewtonRaphson()
        {
            internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();
            
            residual = VectorOperations.VectorVectorSubtraction(BoundaryConditionsImposition.ReducedVector(internalForcesTotalVector, boundaryDof), forceVector);

            residualNorm = VectorOperations.VectorNorm2(residual);
            int iteration = 0;
            while (residualNorm > tolerance && iteration < maxIterations)
            {
                Array.Clear(discretization.TotalStiffnessMatrix, 0, discretization.TotalStiffnessMatrix.Length);
                discretization.CreateTotalStiffnessMatrix();
                double[,] reducedStiffnessMatrix = BoundaryConditionsImposition.ReducedTotalStiff(discretization.TotalStiffnessMatrix, boundaryDof);
                
                IterativeSolver linearSolution = new IterativeSolver(reducedStiffnessMatrix, residual, 1000);
                linearSolution.SolveWithMethod("PCG");
                deltaU = linearSolution.GetSolutionVector;

                double[] reducedSolutionVector = BoundaryConditionsImposition.ReducedVector(solutionVector, boundaryDof);
                reducedSolutionVector = VectorOperations.VectorVectorSubtraction(reducedSolutionVector, deltaU);
                solutionVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(reducedSolutionVector, boundaryDof);

                discretization.UpdateValues(solutionVector);
                Array.Clear(discretization.internalForcesTotalVector, 0, discretization.internalForcesTotalVector.Length);
                internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();

                residual = VectorOperations.VectorVectorSubtraction(BoundaryConditionsImposition.ReducedVector(internalForcesTotalVector, boundaryDof), forceVector);
                
                residualNorm = VectorOperations.VectorNorm2(residual);

                iteration = iteration + 1;
            }
        }

        public void SolveWithMethod(string method)
        {
            switch (method)
            {
                case "Newton-Raphson":
                    NewtonRaphson();
                    break;
            }
        }
    }
}
