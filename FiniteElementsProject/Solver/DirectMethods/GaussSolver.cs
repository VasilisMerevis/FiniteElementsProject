﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class GaussSolver : LinearSolution
    {
        private void GaussElimination(double[,] matrix, double[] vector)
        {
            for (int k = 0; k < vector.Length - 1; k++)
            {
                for (int i = k + 1; i < vector.Length; i++)
                {
                    for (int j = k + 1; j < vector.Length; j++)
                    {
                        matrix[i, j] = matrix[i, j] - (matrix[i, k] / matrix[k, k]) * matrix[k, j];
                    }
                    vector[i] = vector[i] - (matrix[i, k] / matrix[k, k]) * vector[k];
                    matrix[i, k] = 0;
                }
            }
        }

        public override void Solve(double[,] stiffnessMatrix, double[] forceVector)
        {
            linearSolutionVector = new double[forceVector.Length];
            GaussElimination(stiffnessMatrix, forceVector);
            linearSolutionVector = BackSubstitution(stiffnessMatrix, forceVector);
        }

        public override void PrintSolution()
        {
            VectorOperations.PrintVector(linearSolutionVector);
        }



    }
}
