using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject.Solver
{
    public class CentralDifferencesSolver
    {
        public double totalTime;
        public int timeStepsNumber;
        public double time;
        public Dictionary<double,double[]> explicitSolution ;
        private double timeStep;
        private Discretization2DFrame disretization;
        int solutionLength;
        double[,] invertedMassMatrix;
        double[,] stiffenessMatrix;
        double[] externalForcesVector;

        public CentralDifferencesSolver(double totalTime, int timeStepsNumber, Discretization2DFrame discretization)
        {
            this.totalTime = totalTime;
            this.timeStepsNumber = timeStepsNumber;
            timeStep = totalTime / timeStepsNumber;
            invertedMassMatrix = MatrixOperations.InvertDiagMatrix(discretization.TotalMassMatrix);
            solutionLength = discretization.TotalStiffnessMatrix.GetLength(0);
        }

        public void SolveExplicit()
        {
            double[,] twoI = MatrixOperations.CreateDiagonalMatrix(solutionLength, 2);
            double[,] oneI = MatrixOperations.CreateDiagonalMatrix(solutionLength, 1);
            for (int i = 1; i < timeStepsNumber; i++)
            {
                time = i * timeStep;
                double[] solution = new double[solutionLength];
                double[] dt2MR = VectorOperations.VectorScalarProduct(
                                    VectorOperations.MatrixVectorProduct(invertedMassMatrix, externalForcesVector), timeStep * timeStep);
                double[,] dt2MK = MatrixOperations.ScalarMatrixProduct(timeStep * timeStep,
                                    MatrixOperations.MatrixProduct(invertedMassMatrix, stiffenessMatrix));
                double[] dt2MKminus2IbyU = VectorOperations.MatrixVectorProduct(
                                                MatrixOperations.MatrixSubtraction(dt2MK, twoI), solution);
                double[] IbyU = VectorOperations.MatrixVectorProduct(oneI, solution);
                double[] explicitSolutionVector = VectorOperations.VectorVectorSubtraction(dt2MR,
                                                VectorOperations.VectorVectorAddition(dt2MKminus2IbyU, IbyU));

                explicitSolution.Add(time, explicitSolutionVector);
            }
        }
    }
}
