using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    public abstract class Solver
    {
        protected double[] solutionVector;
        protected double[,] stiffnessMatrix;
        protected double[] forceVector;
       

        public double[] GetSolutionVector
        {
            get { return this.solutionVector; }
        }

        abstract public void SolveWithMethod(string method);


    }
}
