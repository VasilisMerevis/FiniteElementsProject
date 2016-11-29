using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class StaticSolver
    {
        private double[] staticSolutionVector;
        private LinearSolution solutionMethod;

        public void SetSolutionMethodToGauss()
        {
            solutionMethod = new GaussSolver();
        }

        public void SetSolutionMethodToCholesky()
        {
            solutionMethod = new CholeskyFactorization();
        }

        public void Solve(double[,] coefMatrix, double[] rhsVector)
        {
            staticSolutionVector = solutionMethod.Solve(coefMatrix, rhsVector);
        }

        public void PrintSolution()
        {
            VectorOperations.PrintVector(staticSolutionVector);
        }

    }
}
