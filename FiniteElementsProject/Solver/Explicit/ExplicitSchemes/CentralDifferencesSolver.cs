﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Media;

namespace FiniteElementsProject
{
    public class CentralDifferencesSolver : ExplicitSolution
    {
        private double totalTime, timeStep;
        private int timeStepsNumber;
        private Dictionary<int,double[]> explicitSolution = new Dictionary<int, double[]>();
        int totalDOFs;
        double[,] massMatrix, dampingMatrix;
        double[,] stiffnessMatrix;
        double[] externalForcesVector;
        double a0, a1, a2, a3;
        double[] initialDisplacementVector, initialVelocityVector, initialAccelerationVector;
        double initialTime;

        public double[] GetExplicitSolution
        {
            get { return explicitSolution[timeStepsNumber]; }
        }

        public CentralDifferencesSolver(double initialTime, double[] initialDisp, double[] initialVel, double[] initialAcc, double totalTime, 
            int timeStepsNumber, double[,] stiffnessMatrix, double[,] massMatrix, double[] externalForcesVector)
        {
            totalDOFs = stiffnessMatrix.GetLength(0);
            this.totalTime = totalTime;
            this.timeStepsNumber = timeStepsNumber;
            timeStep = totalTime / timeStepsNumber;
            this.massMatrix = massMatrix;
            this.stiffnessMatrix = stiffnessMatrix;
            dampingMatrix = new double[totalDOFs, totalDOFs];
            this.externalForcesVector = externalForcesVector;
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
                                    VectorOperations.VectorScalarProductNew(initialVelocityVector, timeStep)),
                                    VectorOperations.VectorScalarProductNew(initialAccelerationVector, a3));
            return previousDisp;
        }

        private double[,] CalculateHatMMatrix()
        {
            double[,] a0M = MatrixOperations.ScalarMatrixProductNew(a0, massMatrix);
            double[,] a1C = MatrixOperations.ScalarMatrixProductNew(a1, dampingMatrix);
            double[,] hutM = MatrixOperations.MatrixAddition(a0M, a1C);
            return hutM;
        }

        private double[,] CalculateHatKMatrix()
        {
            double[,] hatK = MatrixOperations.MatrixSubtraction(stiffnessMatrix,
                                MatrixOperations.ScalarMatrixProductNew(a2, massMatrix));
            return hatK;
        }

        private double[] CalculateHatRVector(int i)
        {
            double[,] hatKMatrix = CalculateHatKMatrix();
            double[,] hatMMatrix = CalculateHatMMatrix();
            double[] hatCurrentU = VectorOperations.MatrixVectorProduct(hatKMatrix, explicitSolution[i-1]);
            double[] hatPreviousU = VectorOperations.MatrixVectorProduct(hatMMatrix, explicitSolution[i - 2]);

            double[] hatR = VectorOperations.VectorVectorSubtraction(externalForcesVector,
                            VectorOperations.VectorVectorAddition(hatCurrentU, hatPreviousU));
            return hatR;
        }

        
        public void SolveExplicit()
        {
            double[,] hatMassMatrix = CalculateHatMMatrix();
            explicitSolution.Add(-1, CalculatePreviousDisplacementVector());
            explicitSolution.Add(0, initialDisplacementVector);
            for (int i = 1; i < timeStepsNumber; i++)
            {
                double time = i * timeStep + initialTime;
                double[] hatRVector = CalculateHatRVector(i);
                DirectSolver linearSolver = new DirectSolver(hatMassMatrix, hatRVector);
                linearSolver.SolveWithMethod("Cholesky");
                double[] nextSolution = linearSolver.GetSolutionVector;
                explicitSolution.Add(i, nextSolution);
            }
        }

        public void PrintExplicitSolution()
        {
            foreach (KeyValuePair<int, double[]> element in explicitSolution)
            {
                int step = element.Key;
                double[] solutionInStep = element.Value;
                Console.WriteLine("Step is{0}", step);
                VectorOperations.PrintVector(solutionInStep);
            }
        }
    }
}
