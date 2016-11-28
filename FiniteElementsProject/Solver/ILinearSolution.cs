using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    interface ILinearSolution
    {
        void Solve(double[,] stiffnessMatrix, double[] forceVector);
        void SetSolutionMethodToGauss();
        void PrintSolution();
    }
}
