using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FiniteElementsProject
{
    class ContactNTN2D : Element1D
    {
        private double PenaltyFactor { get; set; }
        public ContactNTN2D(double E, double A, double[] nodesX, double[] nodesY)
            : base(E, A, nodesX, nodesY)
        {
            lambdaMatrix = new double[4, 4];
            localStiffnessMatrix = new double[6, 6];
            globalStiffnessMatrix = new double[4, 4];
        }

        private double[] CalculateNormalUnitVector()
        {

        }

        public override double[,] CreateGlobalStiffnessMatrix()
        {
            
        }
    }
}
