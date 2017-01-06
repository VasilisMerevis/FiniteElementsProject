using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject.Solver.DirectMethods
{
    public abstract class ExplicitSolution
    {
        public virtual double[] Solve(double[,] stiffnessMatrix, double[] forceVector)
        {
            throw new Exception("Explicit.Solve not implemented");
        }
    }
}
