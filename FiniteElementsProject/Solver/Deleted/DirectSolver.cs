using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementsProject
{
    public class DirectSolver : LinearSolver
    {
		public DirectSolver(double[,] stiffnessMatrix, double[] forceVector)
		{
			this.stiffnessMatrix = stiffnessMatrix;
			this.forceVector = forceVector;
		}

        override public void SolveWithMethod(string method)
		{
			switch(method)
			{
			case "Cholesky":
				double[,] lowerMatrix = Cholesky ();
				double[,] upperMatrix = MatrixOperations.Transpose (lowerMatrix);
				double[] intermediateVector = ForwardSubstitution (lowerMatrix, forceVector);
				double[] solutionuberVector = BackSubstitution (upperMatrix, intermediateVector);
				this.solutionVector = solutionuberVector;
				break;
			case "Gauss":
				GaussElimination (stiffnessMatrix, forceVector);
				this.solutionVector = BackSubstitution (stiffnessMatrix, forceVector);
				break;
			default:
				GaussElimination (stiffnessMatrix, forceVector);
				this.solutionVector = BackSubstitution (stiffnessMatrix, forceVector);
				break;
			}
		}

        private double [,] Cholesky()
        {
            int rows = stiffnessMatrix.GetLength(0);
            int cols = stiffnessMatrix.GetLength(1);

            double[,] lowerPart = new double[rows, cols];
            double sumr;
            double sumc;
            
            for (int j = 0; j < cols; j++)
            {
                sumr = 0;
                for (int k = 0; k <= j-1; k++)
                {
                    sumr = sumr + Math.Pow(lowerPart[j, k], 2);
                }
                if (stiffnessMatrix[j, j] - sumr < 0)
                {
                    throw new Exception("Cholesky: Negative number in square root");
                }
                lowerPart[j, j] = Math.Sqrt(stiffnessMatrix[j, j] - sumr);
                for (int i = j+1; i < rows; i++ )
                {
                    sumc = 0;
                    for (int k = 0; k <= j-1; k++)
                    {
                        sumc = sumc + lowerPart[i, k] * lowerPart[j, k];
                    }
                    lowerPart[i, j] = (stiffnessMatrix[i, j] - sumc) / lowerPart[j, j];
                }
            }
            return lowerPart;
        }

        private void GaussElimination(double[,] matrix, double[] vector)
        {
            for (int k = 0; k < vector.Length-1; k++)
            {
                for (int i = k+1; i < vector.Length; i++)
                {
                    for (int j = k+1; j < vector.Length; j++)
                    {
                        matrix[i, j] = matrix[i, j] - (matrix[i, k] / matrix[k, k]) * matrix[k, j];
                    }
                    vector[i] = vector[i] - (matrix[i, k] / matrix[k, k]) * vector[k];
                    matrix[i, k] = 0;
                }
            }
        }

        private double [] BackSubstitution(double [,] upperTriangMatrix, double [] forceVector)
        {
            int rows = forceVector.GetLength(0);
            double[] solutionVector = new double [rows];
            double total;

            for (int i = rows-1; i >= 0; i-- )
            {
                total = forceVector[i];
                if (i < rows-1)
                {
                    for (int j = i+1; j < rows; j++)
                    {
                        total = total - upperTriangMatrix[i, j] * solutionVector[j];
                    }
                }
                solutionVector[i] = total / upperTriangMatrix[i, i];
            }
            return solutionVector;
        }

        private double [] ForwardSubstitution(double [,] lowerTriangMatrix, double [] forceVector)
        {
            int rows = forceVector.GetLength(0);
            double[] solutionVector = new double[rows];
            double total;
            for (int i = 0; i < rows; i++)
            {
                total = forceVector[i];
                if (i > 0)
                {
                    for (int j = 0; j <= i-1; j++ )
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
