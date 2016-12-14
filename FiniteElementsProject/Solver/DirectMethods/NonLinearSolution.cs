using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    public abstract class NonLinearSolution
    {
        public virtual double[] Solve(double[,] stiffnessMatrix, double[] forceVector)
        {
            throw new Exception("LinearSolution.Solve not implemented");
        }

        
    }
}
