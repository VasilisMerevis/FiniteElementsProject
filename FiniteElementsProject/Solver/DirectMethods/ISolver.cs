using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    interface ISolver
    {
        void SetSolutionMethodToGauss();
        void SetSolutionMethodToCholesky();
        void SetSolutionMethodToPCG();
        void SetNonLinearMethodToLoadControlledNewtonRaphson(int[] boundaryDof, int numberOfLoadSteps, Discretization2DFrame discretization);
        void Solve(double[,] coefMatrix, double[] rhsVector);
        void NLSolve(double[] rhsVector);
        void PrintSolution();
    }
}
