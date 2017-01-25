using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    public class NonLinearSolution : INonLinearSolution
    {
        protected int numberOfLoadSteps = 10;
        protected int[] boundaryDof;
        protected IAssembly discretization;
        protected double lambda;
        protected double tolerance = 1e-5;
        protected int maxIterations = 1000;
        protected ILinearSolution linearSolver;
        //public ILinearSolution LinearScheme { get; set; }
        

        public void DefineBoundaryConditions(int[] boundaryConditionsVector)
        {
            boundaryDof = boundaryConditionsVector;
        }

        public virtual double[] Solve(double[,] stiffnessMatrix, double[] forceVector)
        {
            throw new Exception("LinearSolution.Solve not implemented");
        }

        public virtual double[] NLSolve(IAssembly assembly, ILinearSolution linearScheme, double[] forceVector)
        {
            throw new Exception("LinearSolution.Solve not implemented");
        }


    }
}
