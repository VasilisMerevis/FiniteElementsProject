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
        private int numberOfLoadSteps = 10;
        private double[] dU;
        double[] checking;
        private readonly double[] incrementDf;
        double lambda;

        public NLSolver(Discretization2DFrame discretization, int maxIterations, InputData inputData)
        {
            this.discretization = discretization;
            this.forceVector = inputData.externalForcesVector;
            solutionVector = new double[inputData.nodesX.Length*3];
            this.maxIterations = maxIterations;
            boundaryDof = inputData.boundaryDof;
            this.residual = new double[inputData.nodesX.Length * 3 - boundaryDof.Length];
            //this.df = new double[forceVector.Length];
            this.lambda = 1 / Convert.ToDouble(numberOfLoadSteps);
            this.incrementDf = VectorOperations.VectorScalarProduct(forceVector, lambda);
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
                
                DirectSolver linearSolution = new DirectSolver(reducedStiffnessMatrix, residual);
                linearSolution.SolveWithMethod("Gauss");
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

        public void LoadControlledNewtonRaphson()
        {
            double[] incrementalExternalForcesVector = new double[forceVector.Length];
            double[] tempSolutionVector = new double[solutionVector.Length];
            //double lambda;
            //double[] df = new double[forceVector.Length];
            deltaU = new double[solutionVector.Length - boundaryDof.Length];
            //numberOfLoadSteps = 10;
            //double[] reducedTempSolutionVector = BoundaryConditionsImposition.ReducedVector(tempSolutionVector, boundaryDof);
            //lambda = 1 / Convert.ToDouble(numberOfLoadSteps);
            //df = VectorOperations.VectorScalarProduct(forceVector, lambda);
            for (int i = 0; i < numberOfLoadSteps; i++)
            {
                //lambda = 1 / Convert.ToDouble(numberOfLoadSteps);
                //df = VectorOperations.VectorScalarProduct(forceVector, lambda);
                incrementalExternalForcesVector = VectorOperations.VectorVectorAddition(incrementalExternalForcesVector, incrementDf);


                discretization.UpdateValues(solutionVector);
                Array.Clear(discretization.internalForcesTotalVector, 0, discretization.internalForcesTotalVector.Length);
                internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();

                Array.Clear(discretization.TotalStiffnessMatrix, 0, discretization.TotalStiffnessMatrix.Length);
                discretization.CreateTotalStiffnessMatrix();
                double[,] reducedStiffnessMatrix = BoundaryConditionsImposition.ReducedTotalStiff(discretization.TotalStiffnessMatrix, boundaryDof);

                DirectSolver linearSolution = new DirectSolver(reducedStiffnessMatrix, incrementDf);
                linearSolution.SolveWithMethod("Cholesky");
                dU = linearSolution.GetSolutionVector;

                double[] reducedSolutionVector = BoundaryConditionsImposition.ReducedVector(solutionVector, boundaryDof);
                reducedSolutionVector = VectorOperations.VectorVectorAddition(reducedSolutionVector, dU);
                solutionVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(reducedSolutionVector, boundaryDof);


               

                residual = VectorOperations.VectorVectorSubtraction(BoundaryConditionsImposition.ReducedVector(internalForcesTotalVector, boundaryDof), incrementalExternalForcesVector);

                residualNorm = VectorOperations.VectorNorm2(residual);

                //tempSolutionVector = solutionVector;
                int iteration = 0;
                Array.Clear(deltaU, 0, deltaU.Length);
                while (residualNorm > tolerance && iteration < maxIterations)
                {
                    Array.Clear(discretization.TotalStiffnessMatrix, 0, discretization.TotalStiffnessMatrix.Length);
                    discretization.CreateTotalStiffnessMatrix();
                    reducedStiffnessMatrix = BoundaryConditionsImposition.ReducedTotalStiff(discretization.TotalStiffnessMatrix, boundaryDof);

                    double[] tempResidual = residual;
                    DirectSolver tempLinearSolution = new DirectSolver(reducedStiffnessMatrix, residual);
                    tempLinearSolution.SolveWithMethod("Gauss");
                    double[] addFactorVector = tempLinearSolution.GetSolutionVector;
                    deltaU = VectorOperations.VectorVectorSubtraction(deltaU, addFactorVector);

                    double[] reducedTempSolutionVector = BoundaryConditionsImposition.ReducedVector(tempSolutionVector, boundaryDof);
                    reducedTempSolutionVector = VectorOperations.VectorVectorAddition(reducedSolutionVector, deltaU);
                    tempSolutionVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(reducedTempSolutionVector, boundaryDof);

                    discretization.UpdateValues(tempSolutionVector);
                    Array.Clear(discretization.internalForcesTotalVector, 0, discretization.internalForcesTotalVector.Length);
                    internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();

                    residual = VectorOperations.VectorVectorSubtraction(BoundaryConditionsImposition.ReducedVector(internalForcesTotalVector, boundaryDof), incrementalExternalForcesVector);

                    residualNorm = VectorOperations.VectorNorm2(residual);

                    iteration = iteration + 1;
                }

                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(deltaU, boundaryDof));
                checking = incrementalExternalForcesVector;
            }
            //checking = incrementalExternalForcesVector;
        }

        public void SolveWithMethod(string method)
        {
            switch (method)
            {
                case "Newton-Raphson":
                    NewtonRaphson();
                    break;
                case "Load Controlled Newton-Raphson":
                    LoadControlledNewtonRaphson();
                    break;
            }
        }
    }
}
