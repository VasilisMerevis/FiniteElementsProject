using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    public abstract class LinearSolver
    {
        protected double[] solutionVector;
        protected double[,] stiffnessMatrix;
        protected double[] forceVector;
        public DirectMethods directMethod;    

        public double[] GetSolutionVector
        {
            get { return this.solutionVector; }
        }

        abstract public void SolveWithMethod(string method);

        public void SetSolutionMethodToGauss()
        {
            directMethod = new Gauss();
        }

        public void Solve(double[,] stiffnessMatrix, double[] forceVector)
        {
            directMethod.Solve(stiffnessMatrix, forceVector);
        }

        public void PrintSolution()
        {
            directMethod.PrintDirectSolution();
        }
    }
}
