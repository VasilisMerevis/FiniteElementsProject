using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    public abstract class NonLinearSolution
    {
        protected int numberOfLoadSteps = 10;
        protected int[] boundaryDof;
        protected Discretization2DFrame discretization;
        protected double lambda;
        protected double tolerance = 1e-5;
        protected int maxIterations = 1000;
        protected LinearSolution linearSolver;

        public void DefineBoundaryConditions(int[] boundaryConditionsVector)
        {
            boundaryDof = boundaryConditionsVector;
        }

        public virtual double[] Solve(double[,] stiffnessMatrix, double[] forceVector)
        {
            throw new Exception("LinearSolution.Solve not implemented");
        }

        public virtual double[] NLSolve(double[] forceVector)
        {
            throw new Exception("LinearSolution.Solve not implemented");
        }


    }
}
