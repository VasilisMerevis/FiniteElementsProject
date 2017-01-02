using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class LoadControlledNewtonRaphson : NonLinearSolution
    {
        public LoadControlledNewtonRaphson(IAssembly discretization, LinearSolution linearSolver)
        {
            this.discretization = discretization;
            lambda = 1.0 / numberOfLoadSteps;
            this.linearSolver = linearSolver;
        }
        
        private double[] LoadControlledNR(double[] forceVector)
        {
            double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVector, lambda);
            double[] solutionVector = new double[forceVector.Length];
            double[] incrementalExternalForcesVector = new double[forceVector.Length];
            double[] tempSolutionVector = new double[solutionVector.Length];
            
            double[] deltaU = new double[solutionVector.Length];
            
            double[] internalForcesTotalVector;
            double[] dU;
            double[] residual;
            double residualNorm;
            for (int i = 0; i < numberOfLoadSteps; i++)
            {
                incrementalExternalForcesVector = VectorOperations.VectorVectorAddition(incrementalExternalForcesVector, incrementDf);
                discretization.UpdateValues(solutionVector);
                internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();

                double[,] stiffnessMatrix = discretization.CreateTotalStiffnessMatrix();
                

                dU = linearSolver.Solve(stiffnessMatrix, incrementDf);
                
                
                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, dU);
                

                residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);

                residualNorm = VectorOperations.VectorNorm2(residual);

                int iteration = 0;
                Array.Clear(deltaU, 0, deltaU.Length);
                while (residualNorm > tolerance && iteration < maxIterations)
                {
                    stiffnessMatrix = discretization.CreateTotalStiffnessMatrix();
                    

                    deltaU = VectorOperations.VectorVectorSubtraction(deltaU, linearSolver.Solve(stiffnessMatrix, residual));

                    
                    tempSolutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                    

                    discretization.UpdateValues(tempSolutionVector);
                    internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();

                    residual = VectorOperations.VectorVectorSubtraction(internalForcesTotalVector, incrementalExternalForcesVector);

                    residualNorm = VectorOperations.VectorNorm2(residual);

                    iteration = iteration + 1;
                }

                solutionVector = VectorOperations.VectorVectorAddition(solutionVector, deltaU);
                double[] checking = incrementalExternalForcesVector;
            }
           
            return solutionVector;
        }

        public override double[] NLSolve(double[] forceVector)
        {
            double[] solution = LoadControlledNR(forceVector);
            return solution;
        }

        #region Backup 
        //public LoadControlledNewtonRaphson(Discretization2DFrame discretization, LinearSolution linearSolver)
        //{
        //    this.discretization = discretization;
        //    lambda = 1.0 / numberOfLoadSteps;
        //    this.linearSolver = linearSolver;
        //}

        //private double[] LoadControlledNR(double[] forceVector)
        //{
        //    double[] incrementDf = VectorOperations.VectorScalarProductNew(forceVector, lambda);
        //    double[] solutionVector = new double[forceVector.Length + boundaryDof.Length];
        //    double[] incrementalExternalForcesVector = new double[forceVector.Length];
        //    double[] tempSolutionVector = new double[solutionVector.Length];

        //    double[] deltaU = new double[solutionVector.Length - boundaryDof.Length];

        //    double[] internalForcesTotalVector;
        //    double[] dU;
        //    double[] residual;
        //    double residualNorm;
        //    for (int i = 0; i < numberOfLoadSteps; i++)
        //    {
        //        incrementalExternalForcesVector = VectorOperations.VectorVectorAddition(incrementalExternalForcesVector, incrementDf);
        //        discretization.UpdateValues(solutionVector);
        //        internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();

        //        discretization.CreateTotalStiffnessMatrix();
        //        double[,] reducedStiffnessMatrix = BoundaryConditionsImposition.ReducedTotalStiff(discretization.TotalStiffnessMatrix, boundaryDof);

        //        dU = linearSolver.Solve(reducedStiffnessMatrix, incrementDf);

        //        double[] reducedSolutionVector = BoundaryConditionsImposition.ReducedVector(solutionVector, boundaryDof);
        //        reducedSolutionVector = VectorOperations.VectorVectorAddition(reducedSolutionVector, dU);
        //        solutionVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(reducedSolutionVector, boundaryDof);

        //        residual = VectorOperations.VectorVectorSubtraction(BoundaryConditionsImposition.ReducedVector(internalForcesTotalVector, boundaryDof), incrementalExternalForcesVector);

        //        residualNorm = VectorOperations.VectorNorm2(residual);

        //        int iteration = 0;
        //        Array.Clear(deltaU, 0, deltaU.Length);
        //        while (residualNorm > tolerance && iteration < maxIterations)
        //        {
        //            discretization.CreateTotalStiffnessMatrix();
        //            reducedStiffnessMatrix = BoundaryConditionsImposition.ReducedTotalStiff(discretization.TotalStiffnessMatrix, boundaryDof);

        //            deltaU = VectorOperations.VectorVectorSubtraction(deltaU, linearSolver.Solve(reducedStiffnessMatrix, residual));

        //            double[] reducedTempSolutionVector = BoundaryConditionsImposition.ReducedVector(tempSolutionVector, boundaryDof);
        //            reducedTempSolutionVector = VectorOperations.VectorVectorAddition(reducedSolutionVector, deltaU);
        //            tempSolutionVector = BoundaryConditionsImposition.CreateFullVectorFromReducedVector(reducedTempSolutionVector, boundaryDof);

        //            discretization.UpdateValues(tempSolutionVector);
        //            internalForcesTotalVector = discretization.CreateTotalInternalForcesVector();

        //            residual = VectorOperations.VectorVectorSubtraction(BoundaryConditionsImposition.ReducedVector(internalForcesTotalVector, boundaryDof), incrementalExternalForcesVector);

        //            residualNorm = VectorOperations.VectorNorm2(residual);

        //            iteration = iteration + 1;
        //        }

        //        solutionVector = VectorOperations.VectorVectorAddition(solutionVector, BoundaryConditionsImposition.CreateFullVectorFromReducedVector(deltaU, boundaryDof));
        //        double[] checking = incrementalExternalForcesVector;
        //    }

        //    return solutionVector;
        //}

        //public override double[] NLSolve(double[] forceVector)
        //{
        //    double[] solution = LoadControlledNR(forceVector);
        //    return solution;
        //}
        #endregion
    }
}
