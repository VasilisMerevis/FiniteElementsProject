using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    public class LinearSolution : ILinearSolution
    {
        protected double[] linearSolutionVector;
        //protected double[,] stiffnessMatrix;
        //protected double[] forceVector;
        public LinearSolution solutionMethod;

        public double[] GetSolutionVector
        {
            get { return this.linearSolutionVector; }
        }

        

        public void SetSolutionMethodToGauss()
        {
            solutionMethod = new GaussSolver();
        }

        public virtual void Solve(double[,] stiffnessMatrix, double[] forceVector)
        {
            solutionMethod.Solve(stiffnessMatrix, forceVector);
        }

        public virtual void PrintSolution()
        {
            solutionMethod.PrintSolution();
        }

        protected double[] BackSubstitution(double[,] upperTriangMatrix, double[] forceVector)
        {
            int rows = forceVector.GetLength(0);
            double[] solutionVector = new double[rows];
            double total;

            for (int i = rows - 1; i >= 0; i--)
            {
                total = forceVector[i];
                if (i < rows - 1)
                {
                    for (int j = i + 1; j < rows; j++)
                    {
                        total = total - upperTriangMatrix[i, j] * solutionVector[j];
                    }
                }
                solutionVector[i] = total / upperTriangMatrix[i, i];
            }
            return solutionVector;
        }

        protected double[] ForwardSubstitution(double[,] lowerTriangMatrix, double[] forceVector)
        {
            int rows = forceVector.GetLength(0);
            double[] solutionVector = new double[rows];
            double total;
            for (int i = 0; i < rows; i++)
            {
                total = forceVector[i];
                if (i > 0)
                {
                    for (int j = 0; j <= i - 1; j++)
                    {
                        total = total - lowerTriangMatrix[i, j] * solutionVector[j];
                    }
                }
                solutionVector[i] = total / lowerTriangMatrix[i, i];
            }
            return solutionVector;
        }
    }
}
