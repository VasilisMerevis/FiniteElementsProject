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
        void SetNonLinearMethodToLoadControlledNewtonRaphson();
        void Solve(double[,] coefMatrix, double[] rhsVector);
        void NLSolve(double[] rhsVector);
        void ReadBoundaryConditions(int[] boundaryCond);
        void PrintSolution();

        ILinearSolution LinearScheme { get; set; }
        IAssembly AssemblyData { get; set; }
        void SolveStatic(double[] rhsVector);
        bool ActivateNonLinearSolver { get; set; }

    }
}
