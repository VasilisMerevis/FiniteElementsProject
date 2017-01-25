using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class StaticSolver : ISolver
    {
        private double[] staticSolutionVector;
        
        private NonLinearSolution nonLinearSolutionMethod;
        public ILinearSolution LinearScheme { get; set; }
        public INonLinearSolution NonLinearScheme { get; set; }
        public bool ActivateNonLinearSolver { get;  set; }
        public IAssembly AssemblyData { get; set; }

        public StaticSolver()
        {
            ActivateNonLinearSolver = false;
        }
        
        //public void SetSolutionMethodToGauss()
        //{
        //    solutionMethod = new GaussSolver();
        //}

        //public void SetSolutionMethodToCholesky()
        //{
        //    solutionMethod = new CholeskyFactorization();
        //}

        //public void SetSolutionMethodToPCG()
        //{
        //    solutionMethod = new PCGSolver();
        //}

        //public void SetNonLinearMethodToLoadControlledNewtonRaphson()
        //{
        //    nonLinearSolutionMethod = new LoadControlledNewtonRaphson(AssemblyData, LinearScheme);
        //}

        public void ReadBoundaryConditions(int[] boundaryCond)
        {
            nonLinearSolutionMethod.DefineBoundaryConditions(boundaryCond);
        }

        //public void Solve(double[,] coefMatrix, double[] rhsVector)
        //{
        //    staticSolutionVector = solutionMethod.Solve(coefMatrix, rhsVector);
        //}

        //public void NLSolve(double[] rhsVector)
        //{
        //    staticSolutionVector = nonLinearSolutionMethod.NLSolve(rhsVector);
        //}

        public void SolveStatic(double[] rhsVector)
        {
            if (ActivateNonLinearSolver == true)
            {
                //staticSolutionVector = nonLinearSolutionMethod.NLSolve(AssemblyData, LinearScheme, rhsVector);
                staticSolutionVector = NonLinearScheme.NLSolve(AssemblyData, LinearScheme, rhsVector);
            }
            else
            {
                double[,] coefMatrix = AssemblyData.CreateTotalStiffnessMatrix();
                staticSolutionVector = LinearScheme.Solve(coefMatrix, rhsVector);
            }
        }

        public void PrintSolution()
        {
            VectorOperations.PrintVector(staticSolutionVector);
        }

    }
}
