using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    public class DirectMethods : LinearSolver
    {
        //protected double[] solutionVector;
        //protected double[,] stiffnessMatrix;
        //protected double[] forceVector;

        

        override public void SolveWithMethod(string method)
        {

        }

        public virtual void PrintDirectSolution()
        {
            
        }

        

        public virtual void Solve(double[,] stiffnessMatrix, double[] forceVector)
        {
            
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
