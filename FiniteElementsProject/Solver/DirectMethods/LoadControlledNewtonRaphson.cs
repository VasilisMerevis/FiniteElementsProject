using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class LoadControlledNewtonRaphson : NonLinearSolution
    {
        int numberOfLoadSteps;
        int[] boundaryDof;
        Discretization2DFrame discretization;
        double lambda;
        double tolerance = 1e-5;
        int maxIterations = 1000;


        public LoadControlledNewtonRaphson(int[] boundaryDof, int numberOfLoadSteps, Discretization2DFrame discretization)
        {
            this.boundaryDof = boundaryDof;
            this.numberOfLoadSteps = numberOfLoadSteps;
            this.discretization = discretization;
            this.lambda = 1 / numberOfLoadSteps;
            
        }
        
        public void LoadControlledNR(double[] forceVector)
        {
            double[] incrementDf = VectorOperations.VectorScalarProduct(forceVector, lambda);
            double[] solutionVector = new double[forceVector.Length];
            double[] incrementalExternalForcesVector = new double[forceVector.Length];
            double[] tempSolutionVector = new double[solutionVector.Length];
            //double lambda;
            //double[] df = new double[forceVector.Length];
            double[] deltaU = new double[solutionVector.Length - boundaryDof.Length];
            //numberOfLoadSteps = 10;
            //double[] reducedTempSolutionVector = BoundaryConditionsImposition.ReducedVector(tempSolutionVector, boundaryDof);
            //lambda = 1 / Convert.ToDouble(numberOfLoadSteps);
            //df = VectorOperations.VectorScalarProduct(forceVector, lambda);
            double[] internalForcesTotalVector;
            double[] dU;
            double[] residual;
            double residualNorm;
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
                double[] checking = incrementalExternalForcesVector;
            }
            //checking = incrementalExternalForcesVector;
        }

        
    }
}
