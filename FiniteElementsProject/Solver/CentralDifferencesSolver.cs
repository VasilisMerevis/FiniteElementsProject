using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Media;

namespace FiniteElementsProject
{
    public class CentralDifferencesSolver
    {
        private double totalTime, timeStep;
        private int timeStepsNumber;
        private Dictionary<int,double[]> explicitSolution = new Dictionary<int, double[]>();
        int totalDOFs;
        double[,] massMatrix, dampingMatrix;
        double[,] stiffenessMatrix;
        double[] externalForcesVector = new double[12];
        double a0, a1, a2, a3;
        double[] initialDisplacementVector, initialVelocityVector, initialAccelerationVector;
        double initialTime;

        public CentralDifferencesSolver(double initialTime, double[] initialDisp, double[] initialVel, double[] initialAcc, double totalTime, int timeStepsNumber, double[,] stiffnessMatrix, double[,] massMatrix)
        {
            totalDOFs = stiffnessMatrix.GetLength(0);
            this.totalTime = totalTime;
            this.timeStepsNumber = timeStepsNumber;
            timeStep = totalTime / timeStepsNumber;
            this.massMatrix = massMatrix;
            this.stiffenessMatrix = stiffnessMatrix;
            dampingMatrix = new double[totalDOFs, totalDOFs];
            initialDisplacementVector = initialDisp;
            initialVelocityVector = initialVel;
            initialAccelerationVector = initialAcc;
            this.initialTime = initialTime;
            
            a0 = 1 / (timeStep * timeStep);
            a1 = 1 / (2 * timeStep);
            a2 = 2 * a0;
            a3 = 1 / a2;
        }

        private double[] CalculatePreviousDisplacementVector()
        {
            double[] previousDisp = VectorOperations.VectorVectorAddition(
                                    VectorOperations.VectorVectorSubtraction(initialDisplacementVector,
                                    VectorOperations.VectorScalarProduct(initialVelocityVector, timeStep)),
                                    VectorOperations.VectorScalarProduct(initialAccelerationVector, a3));
            return previousDisp;
        }

        private double[,] CalculateHatMMatrix()
        {
            double[,] a0M = MatrixOperations.ScalarMatrixProduct(a0, massMatrix);
            double[,] a1C = MatrixOperations.ScalarMatrixProduct(a1, dampingMatrix);
            double[,] hutM = MatrixOperations.MatrixAddition(a0M, a1C);
            return hutM;
        }

        private double[,] CalculateHatKMatrix()
        {
            double[,] hatK = MatrixOperations.MatrixSubtraction(stiffenessMatrix,
                                MatrixOperations.ScalarMatrixProduct(a2, massMatrix));
            return hatK;
        }

        private double[] CalculateHatRVector(int timeIncrement)
        {
            double[,] hatKMatrix = CalculateHatKMatrix();
            double[,] hatMMatrix = CalculateHatMMatrix();
            double[] hatCurrentU = VectorOperations.MatrixVectorProduct(hatKMatrix, explicitSolution[timeIncrement]);
            double[] hatPreviousU = VectorOperations.MatrixVectorProduct(hatMMatrix, explicitSolution[timeIncrement - 1]);

            double[] hatR = VectorOperations.VectorVectorSubtraction(externalForcesVector,
                            VectorOperations.VectorVectorAddition(hatCurrentU, hatPreviousU));
            return hatR;
        }

        
        public void SolveExplicit()
        {
            //double[,] twoI = MatrixOperations.CreateDiagonalMatrix(solutionLength, 2);
            //double[,] oneI = MatrixOperations.CreateDiagonalMatrix(solutionLength, 1);

            double[,] hatMassMatrix = CalculateHatMMatrix();
            double[] previousDisplacement = CalculatePreviousDisplacementVector();
            explicitSolution.Add(-1, previousDisplacement);
            explicitSolution.Add(0, initialDisplacementVector);
            for (int i = 0; i < timeStepsNumber; i++)
            {
                double time = i * timeStep + initialTime;
                double[] hatRVector = CalculateHatRVector(i);
                IterativeSolver linearSolver = new IterativeSolver(hatMassMatrix, hatRVector, 1000);
                linearSolver.SolveWithMethod("PCG");

                //DirectSolver linearSolver = new DirectSolver(hatMassMatrix, hatRVector);
                //linearSolver.SolveWithMethod("Cholesky");

                double[] nextSolution = linearSolver.GetSolutionVector;
                explicitSolution.Add(i+1, nextSolution);

                //double[] solution = new double[solutionLength];
                //double[] dt2MR = VectorOperations.VectorScalarProduct(
                //                    VectorOperations.MatrixVectorProduct(invertedMassMatrix, externalForcesVector), timeStep * timeStep);
                //double[,] dt2MK = MatrixOperations.ScalarMatrixProduct(timeStep * timeStep,
                //                    MatrixOperations.MatrixProduct(invertedMassMatrix, stiffenessMatrix));
                //double[] dt2MKminus2IbyU = VectorOperations.MatrixVectorProduct(
                //                                MatrixOperations.MatrixSubtraction(dt2MK, twoI), solution);
                //double[] IbyU = VectorOperations.MatrixVectorProduct(oneI, solution);
                //double[] explicitSolutionVector = VectorOperations.VectorVectorSubtraction(dt2MR,
                //                                VectorOperations.VectorVectorAddition(dt2MKminus2IbyU, IbyU));

                //explicitSolution.Add(time, explicitSolutionVector);
            }
        }
    }
}
