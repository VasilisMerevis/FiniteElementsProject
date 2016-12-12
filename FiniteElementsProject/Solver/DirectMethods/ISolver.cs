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
        void Solve(double[,] coefMatrix, double[] rhsVector);
        void PrintSolution();
    }
}
