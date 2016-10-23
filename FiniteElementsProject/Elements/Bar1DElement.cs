using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementsProject
{
    class Bar1DElement : Element1D
    {
        
        public double density = 8000;

        public Bar1DElement(double E, double A, double[] nodesX, double[] nodesY)
            :base(E, A, nodesX, nodesY)
        {
            lambdaMatrix = new double[6, 6];
            localStiffnessMatrix = new double[6, 6];
            globalStiffnessMatrix = new double[6, 6];
        }

        public override double[,] CreateLocalStiffnessMatrix()
        {
            double length = VectorOperations.CalculateVectorLengthFromEndPoints (nodesX [0], nodesX [1], nodesY [0], nodesY [1]);
            localStiffnessMatrix = new[,]
            { { E * A / length, 0, 0, -E * A / length, 0, 0 },
                { 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0},
                { -E * A / length, 0, 0, E* A/ length, 0, 0 },
                { 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0}};

            return localStiffnessMatrix;
        }

        public override double[,] CreateLambdaMatrix()
        {
            double length = Math.Sqrt (Math.Pow ((nodesY [1] - nodesY [0]), 2) + Math.Pow ((nodesX [1] - nodesX [0]), 2));
            double cosine = (nodesX [1] - nodesX [0]) / length;
            double sine = (nodesY [1] - nodesY [0]) / length;
            lambdaMatrix = new[,]
            { { cosine, sine, 0, 0, 0, 0 },
                { -sine, cosine, 0, 0, 0, 0 },
                { 0, 0, 1, 0, 0, 0 },
                { 0, 0, 0, cosine, sine, 0 },
                { 0, 0, 0, -sine, cosine, 0 },
                { 0, 0, 0, 0, 0, 1 }};

            return lambdaMatrix;
        }

        public override double[,] CreateMassMatrix()
        {
            double length = Math.Sqrt(Math.Pow((nodesY[1] - nodesY[0]), 2) + Math.Pow((nodesX[1] - nodesX[0]), 2));
            double elementMass = density * A * length;
            massMatrix = new double[,]
            {
                {elementMass/2, 0, 0, 0, 0, 0},
                {0, elementMass/2, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, elementMass/2, 0, 0},
                {0, 0, 0, 0, elementMass/2, 0},
                {0, 0, 0, 0, 0, 0}
            };
            
            return massMatrix;
        }
    }
}
