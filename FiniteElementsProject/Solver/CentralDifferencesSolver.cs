using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject.Solver
{
    public class CentralDifferencesSolver
    {
        public double totalTime;
        public int timeStepsNumber;
        public double time;
        public Dictionary<double,double[]> explicitSolution;
        private double timeStep;

        public CentralDifferencesSolver(double totalTime, int timeStepsNumber)
        {
            this.totalTime = totalTime;
            this.timeStepsNumber = timeStepsNumber;
            timeStep = totalTime / timeStepsNumber;
        }

        public void SolveExplicit()
        {
            for (int i = 0; i < timeStepsNumber; i++)
            {

            }
        }
    }
}
